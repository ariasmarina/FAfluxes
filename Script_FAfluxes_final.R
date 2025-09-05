# R Script for: Arias et al. 
# Title: Multiple stressors across ecosystem boundaries: do invasive species and 
# light pollution change the quality of aquatic subsidies to terrestrial 
# predators?

# Authors: Marina Arias, Gemma Burgazzi, Sebastian Pietz, et al. 

# The code has been written by: 
# Dr. Marina Arias
# Institute for Environmental Sciences - RPTU. Landau, Germany &
# Instituto de Limnologia Dr. Raul A. Ringuelet. La Plata, Argentina
# Email: arias@ilpla.edu.ar

# Revised by: Dr. Gemma Burgazzi, University of Turin, Italy.

# R version 4.1.2 (2021) -- "Bird Hippie"

# In this script you can find the analysis and plots done for the manuscript.

# This study addresses the individual and joint effect of ALAN and signal 
# crayfish on the FA fluxes from aquatic to terrestrial ecosystems in a fully 
# crossed experimental design employed at the Riparian Stream Mesocosm (RSM) 
# facility in Landau, Germany. 

# Code sections:
# 1. Load libraries
# 2. Emergent aquatic insects
## 2.1. Data import and manipulation
## 2.2. Models
## 2.3. Summary
## 2.4 Plots
# 3. FA concentrations of chironomids
## 3.1. Data import and manipulation
## 3.2. Models
## 3.3. Summary
## 3.4 Plots
# 4. FA concentrations of spiders 
## 4.1. Data import and manipulation
## 4.2. Models
## 4.3. Summary
## 4.4 Plots
# 5. Supplementary Information
## 5.1. Signal Crayfish activity under ALAN
## 5.2. FA content in chironomids
## 5.3. FA content in spiders

#1. Load libraries ----
setwd("C:/R FA Analysis")

library(car) #version 3.1-1, checking linear models assumptions
library(cowplot) #version 1.1.3, arrange figures
library(gridExtra) #version 2.3, organise figures
library(DHARMa) #version 0.4.7. generalised model diagnostics
library(dplyr) #version 1.1.3, data manipulation
library(emmeans) #version 1.8.6 posthoc test
library(ggplot2) #version 3.4.3, prepare figures
library(lme4) #version 1.1-27.1, mixed models
library(performance) #version 0.13.0, checking mixed models assumptions
library(see) #version 0.9.0, visualization of diagnosis graphs 
library(tidyr) #version 1.3.0, data manipulation
library(tidyverse) #version 2.0.0, data manipulation 

# 2. Emergent aquatic insects ----

## 2.1. Data import and manipulation ----
taxa <- read.csv("emergence_biomass.csv", h=T, sep=";")

#set factors
taxa$treat<-factor(taxa$treat, levels = c("Control", "ALAN", "Crayfish", "ALAN+Crayfish"))
taxa$trap <- factor(taxa$trap, levels = c("up", "mid", "down"))
taxa$week <- factor(taxa$week, levels = c("Week1", "Week4", "Week6"))

# Calculations of biomass rate (mean mg/m2/d)
taxa_red <- taxa[,-3] #delete the grouping factor "day"

#obtain total mg per week
taxa_bio <- taxa_red %>% group_by(treat, week, flume, trap) %>% summarise(
  bio.week = sum(chiro_bio))

#calculate biomass per square meter
taxa_bio2 <- taxa_bio #copy of the dataset
taxa_bio2$bio.week[taxa_bio$trap == "up"] <- taxa_bio2$bio.week[taxa_bio$trap == "up"]/0.7  
taxa_bio2$bio.week[taxa_bio$trap == "mid"] <- taxa_bio2$bio.week[taxa_bio$trap == "mid"]/0.75
taxa_bio2$bio.week[taxa_bio$trap == "down"] <- taxa_bio2$bio.week[taxa_bio$trap == "down"]/0.8

#calculate the rate per day (mg/m2/d)
taxa_bio3 <- taxa_bio2
taxa_bio3$bio.week.m2.day <- taxa_bio3$bio.week/7  

#average the three traps per flume
taxa_bio4 <- taxa_bio3[,-4]
taxa_bio4 <- taxa_bio3 %>% group_by(treat, week, flume) %>% summarise(
  mean.bio = mean(bio.week.m2.day))

#create final data frame
bio_mean <- taxa_bio4
bio_mean$mean.bio <- taxa_bio4$mean.bio 

write.csv(bio_mean, "biomass_mean.csv") #print

chiro <- read.csv("biomass_mean.csv", h=T, sep=",") #visualization 

## 2.2. Models ----

#set dummy variables for stressors
chiro$ALAN <- chiro$treat
chiro$Crayfish <- chiro$treat

chiro$ALAN <- plyr::mapvalues(chiro$ALAN, from=c("ALAN+Crayfish", "Crayfish", "ALAN", "Control"), to=c(1,0,1,0))
chiro$Crayfish <- plyr::mapvalues(chiro$Crayfish, from=c("ALAN+Crayfish", "Crayfish", "ALAN", "Control"), to=c(1,1,0,0))

#exploration of variable distribution
hist(chiro$mean.bio)
hist(sqrt(chiro$mean.bio)) 

#fit model
lme_bio<- lmer(sqrt(mean.bio) ~ ALAN * Crayfish * week + (1|flume), data = chiro)

#model check
check_overdispersion(lme_bio) 
text(qqnorm(resid(lme_bio)))
qqline(resid(lme_bio))
check_collinearity(lme_bio)
check_model(lme_bio) 

# model summary and anova
summary(lme_bio)
Anova(lme_bio, type = 2)

# R2
r2_nakagawa(lme_bio)

## 2.3 Summary ----

# Mean, standard deviation, 95% confidence intervals and effect sizes for each treatment
summary_chiro <- chiro %>%
  group_by(treat) %>%
  summarise(
    mean_value = mean(mean.bio),
    sd_value = sd(mean.bio),
    count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lower_ci = mean_value - qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    upper_ci = mean_value + qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    effect = case_when(
      treat == "ALAN" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "ALAN+Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      TRUE ~ NA_real_
    ))

#Summary per week
#separation per time point
chiro_w1 <- subset(chiro, chiro$week=="Week1")
chiro_w4 <- subset(chiro, chiro$week=="Week4")
chiro_w6 <- subset(chiro, chiro$week=="Week6")

#CHIRO W1
summary_chiro_w1 <- chiro_w1 %>%
  group_by(treat) %>%
  summarise(
    mean_value = mean(mean.bio),
    sd_value = sd(mean.bio),
    count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lower_ci = mean_value - qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    upper_ci = mean_value + qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    effect = case_when(
      treat == "ALAN" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "ALAN+Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      TRUE ~ NA_real_
    ))

#CHIRO W4
summary_chiro_w4 <- chiro_w4 %>%
  group_by(treat) %>%
  summarise(
    mean_value = mean(mean.bio),
    sd_value = sd(mean.bio),
    count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lower_ci = mean_value - qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    upper_ci = mean_value + qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    effect = case_when(
      treat == "ALAN" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "ALAN+Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      TRUE ~ NA_real_
    ))

#CHIRO W6
summary_chiro_w6 <- chiro_w6 %>%
  group_by(treat) %>%
  summarise(
    mean_value = mean(mean.bio),
    sd_value = sd(mean.bio),
    count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lower_ci = mean_value - qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    upper_ci = mean_value + qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    effect = case_when(
      treat == "ALAN" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "ALAN+Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      TRUE ~ NA_real_
    ))

## 2.4. Plots ----
#Figure 2

#set factors
chiro$treat<-factor(chiro$treat, levels = c("Control", "ALAN", "Crayfish", "ALAN+Crayfish"))
chiro$week <- factor(chiro$week, levels = c("Week1", "Week4", "Week6"))

#Figure2
#to export the figure:
jpeg('Figure2.jpeg', width = 10, height = 6, units = 'in', res=300)
Figure2 <- ggplot(chiro, aes(x = treat, y = mean.bio)) +
  geom_jitter(aes(color = treat), width = 0.2, alpha = 0.6, size = 3) +
  stat_summary(fun = median, geom = "crossbar", width = 0.4, fatten = 3, color = "black") + 
  stat_summary(fun.min = min, fun.max = max, fun = median,
               geom = "errorbar", width = 0.2, color = "black") +
  scale_color_manual(values = c("#303030", "#FFB90F","#009ACD", "#5D478B")) +
  ylab(bquote("Emergence biomass " ~ (mg * ~ m^{-2} ~ d^{-1}))) +
  facet_wrap(~week) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size = 14, face= "bold"),
        legend.position = "none") + 
  scale_x_discrete(labels = c("Control", "ALAN", "Crayfish", "ALAN+Crayfish"))
Figure2
dev.off()

# 3. FA composition of chironomids ----
## 3.1. Data import and manipulation ----
FAgroups <- read.table ("FA_emergence.csv",h=T,sep=";" ,fill = TRUE)

#tidy dataset deleting extra columns
FAgroups <-FAgroups[,-27]
FAgroups <-FAgroups[,-32]

#set factors
FAgroups$treat<- factor(FAgroups$treat, levels = c("Control", "ALAN", "Crayfish", "ALAN+Crayfish"))
FAgroups$week <- factor(FAgroups$week, levels = c("Week1", "Week4", "Week6"))

#reorganize dataset for calculations by grouping FAs in categories
FAgroups$SFA <- apply(FAgroups[,4:15], 1, sum) #saturated FA
FAgroups$MUFA <- apply(FAgroups[,16:23], 1, sum) #monounsaturated FA
FAgroups$PUFA <- apply(FAgroups[24:33], 1, sum) #polyunsaturated FA (PUFA)
FAgroups$totalFA <- apply(FAgroups[,34:36], 1, sum) #total FA

#copy for Supplementary Material 4
sm4 <- FAgroups [,-c(4:29, 32, 33)]

#Calculation FA fluxes
#security copy
chiro2<- chiro  #mean biomass
chiro2<- chiro2[,-1] #tidy extra columns
chiro2$treat<- factor(chiro2$treat, levels = c("Control", "ALAN", "Crayfish", "ALAN+Crayfish"))
chiro2$week <- factor(chiro2$week, levels = c("Week1", "Week4", "Week6"))
chiro2 <- chiro2[, c(2, 1, 3, 4)] #reorganize columns

fluxes <- FAgroups[,-(4:33)]

#merge data
fluxes_fa<- merge(fluxes, chiro2, by = c("week", "treat", "flume"))

#add physiologically important PUFA
fluxes_phy <- FAgroups %>%
  select(week, treat, flume, C20.4_ARA, C20.5_EPA)

#merge data
fluxes_fa<- merge(fluxes_phy, fluxes_fa, by = c("week", "treat", "flume"))

#fluxes calculations
fluxes_fa$SFAflux <- fluxes_fa$SFA * fluxes_fa$mean.bio
fluxes_fa$MUFAflux <- fluxes_fa$MUFA * fluxes_fa$mean.bio
fluxes_fa$PUFAflux <- fluxes_fa$PUFA * fluxes_fa$mean.bio
fluxes_fa$totalFAflux <- fluxes_fa$totalFA * fluxes_fa$mean.bio
fluxes_fa$ARAflux <- fluxes_fa$C20.4_ARA * fluxes_fa$mean.bio
fluxes_fa$EPAflux <- fluxes_fa$C20.5_EPA * fluxes_fa$mean.bio

flux <- fluxes_fa[,-c(4:10)]

#export new dataset
write.csv(flux, "fluxes.csv")

## 3.2. Models ----
#set dummy variables
flux$ALAN <- flux$treat 
flux$Crayfish <- flux$treat

flux$ALAN <- plyr::mapvalues(flux$ALAN, from=c("ALAN+Crayfish", "Crayfish", "ALAN", "Control"), to=c(1,0,1,0))
flux$Crayfish <- plyr::mapvalues(flux$Crayfish, from=c("ALAN+Crayfish", "Crayfish", "ALAN", "Control"), to=c(1,1,0,0))

## totalFA flux  ----

#exploration of variable distribution
hist(flux$totalFAflux)
hist(log(flux$totalFAflux)+1) 

#mixed model total fa
lme_fa<- lmer(log((totalFAflux)+1) ~ ALAN * Crayfish * week + (1|flume), data=flux)

#model check
check_overdispersion(lme_fa) 
text(qqnorm(resid(lme_fa)))
qqline(resid(lme_fa))
check_collinearity(lme_fa)
check_model(lme_fa) 

#model summary and anova
summary(lme_fa)
Anova(lme_fa, type=2) #Week p<0.001, Crayfish:Week p=0.016

# R2
r2_nakagawa(lme_fa)

##SFA flux ----

#exploration of variable distribution
hist(flux$SFAflux)
hist(log(flux$SFAflux)+1) 

#mixed model SFA
lme_sfa<- lmer(log((SFAflux)+1)~ ALAN * Crayfish * week + (1|flume), data=flux)

#model check
check_overdispersion(lme_sfa) 
text(qqnorm(resid(lme_sfa)))
qqline(resid(lme_sfa))
check_collinearity(lme_sfa)
check_model(lme_sfa) 

#model summary and anova
summary(lme_sfa)
Anova(lme_sfa, type=2) #Week p<0.001, Crayfish:Week p=0.007

# R2
r2_nakagawa(lme_sfa)

##MUFA flux -----

#exploration of variable distribution
hist(flux$MUFAflux)
hist(log(flux$MUFAflux)+1)

#mixed model MUFA
lme_mufa<- lmer(log((MUFAflux)+1) ~ ALAN * Crayfish * week + (1|flume), data=flux)

#model check
check_overdispersion(lme_mufa) 
text(qqnorm(resid(lme_mufa)))
qqline(resid(lme_mufa))
check_collinearity(lme_mufa)
check_model(lme_mufa) 

#model summary and anova
summary(lme_mufa)
Anova(lme_mufa, type=2) #Week p<0.001, Crayfish:Week p=0.02

# R2
r2_nakagawa(lme_mufa)

##PUFA flux ----

#exploration of variable distribution
hist(flux$PUFAflux)
hist(log(flux$PUFAflux)+1)

#mixed model PUFA fa
lme_pufa<- lmer(log((PUFAflux)+1) ~ ALAN * Crayfish * week + (1|flume), data=flux)

#model check
check_overdispersion(lme_pufa) 
text(qqnorm(resid(lme_pufa)))
qqline(resid(lme_pufa))
check_collinearity(lme_pufa)
check_model(lme_pufa) 

#model summary and anova
summary(lme_pufa)
Anova(lme_pufa, type=2) #Week p<0.001, Crayfish:Week p=0.036

# R2
r2_nakagawa(lme_pufa)

##Physiologically important PUFA fluxes----

#ARA
#exploration of variable distribution
hist(flux$ARAflux)
hist(log(flux$ARAflux)+1)

#mixed model PUFA fa
lme_ara<- lmer(log((ARAflux)+1) ~ ALAN * Crayfish * week + (1|flume), data=flux)

#model check
check_overdispersion(lme_ara) 
text(qqnorm(resid(lme_ara)))
qqline(resid(lme_ara))
check_collinearity(lme_ara)
check_model(lme_ara) 

#model summary and anova
summary(lme_ara)
Anova(lme_ara, type=2) #Week p<0.001

# R2
r2_nakagawa(lme_ara)

#EPA

#exploration of variable distribution
hist(flux$EPAflux)
hist(log(flux$EPAflux)+1)

#fit model
lme_epa<- lmer(log((EPAflux)+1) ~ ALAN * Crayfish * week + (1|flume), data=flux)

#model check
check_overdispersion(lme_epa) 
text(qqnorm(resid(lme_epa)))
qqline(resid(lme_epa))
check_collinearity(lme_epa)
check_model(lme_epa) 

#model summary and anova
summary(lme_epa)
Anova(lme_epa, type=2) #Week p<0.001; Crayfish:Week p=0.047

# R2
r2_nakagawa(lme_epa_mix)

#For all FA groups, ARA and EPA, factor "week" was significant (changes over time) 
#and, except for ARA, the interaction "Crayfish:Week" (presence of crayfish over time) 

## 3.3. Summary ----

# Mean, standard deviation and effect sizes for each treatment

#total FA flux
summary_faflux <- flux %>%
  group_by(treat) %>%
  summarise(
    mean_value = mean(totalFAflux),
    sd_value = sd(totalFAflux),
    count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lower_ci = mean_value - qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    upper_ci = mean_value + qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    effect = case_when(
      treat == "ALAN" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "ALAN+Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      TRUE ~ NA_real_
    ))

#PUFA flux
summary_pufaflux <- flux %>%
  group_by(treat) %>%
  summarise(
    mean_value = mean(PUFAflux),
    sd_value = sd(PUFAflux),
    count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lower_ci = mean_value - qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    upper_ci = mean_value + qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    effect = case_when(
      treat == "ALAN" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "ALAN+Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      TRUE ~ NA_real_
    ))

#Temporal pattern in PUFA among treatments 
#subset pero week
flux_w1 <- subset(flux, flux$week=="Week1")
flux_w4 <- subset(flux, flux$week=="Week4")
flux_w6 <- subset(flux, flux$week=="Week6")

#PUFA W1
summary_pufa_w1 <- flux_w1 %>%
  group_by(treat) %>%
  summarise(
    mean_value = mean(PUFAflux),
    sd_value = sd(PUFAflux),
    count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lower_ci = mean_value - qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    upper_ci = mean_value + qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    effect = case_when(
      treat == "ALAN" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "ALAN+Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      TRUE ~ NA_real_
    ))

#PUFA W4
summary_pufa_w4 <- flux_w4 %>%
  group_by(treat) %>%
  summarise(
    mean_value = mean(PUFAflux),
    sd_value = sd(PUFAflux),
    count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lower_ci = mean_value - qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    upper_ci = mean_value + qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    effect = case_when(
      treat == "ALAN" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "ALAN+Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      TRUE ~ NA_real_
    ))

#PUFA W6
summary_pufa_w6 <- flux_w6 %>%
  group_by(treat) %>%
  summarise(
    mean_value = mean(PUFAflux),
    sd_value = sd(PUFAflux),
    count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lower_ci = mean_value - qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    upper_ci = mean_value + qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    effect = case_when(
      treat == "ALAN" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "ALAN+Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      TRUE ~ NA_real_
    ))

#ARA flux w1
summary_ara_w1 <- flux_w1 %>%
  group_by(treat) %>%
  summarise(
    mean_value = mean(ARAflux),
    sd_value = sd(ARAflux),
    count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lower_ci = mean_value - qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    upper_ci = mean_value + qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    effect = case_when(
      treat == "ALAN" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "ALAN+Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      TRUE ~ NA_real_
    ))

#ARA flux w4
summary_ara_w4 <- flux_w4 %>%
  group_by(treat) %>%
  summarise(
    mean_value = mean(ARAflux),
    sd_value = sd(ARAflux),
    count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lower_ci = mean_value - qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    upper_ci = mean_value + qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    effect = case_when(
      treat == "ALAN" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "ALAN+Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      TRUE ~ NA_real_
    ))

#ARA flux w6
summary_ara_w6 <- flux_w6 %>%
  group_by(treat) %>%
  summarise(
    mean_value = mean(ARAflux),
    sd_value = sd(ARAflux),
    count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lower_ci = mean_value - qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    upper_ci = mean_value + qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    effect = case_when(
      treat == "ALAN" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "ALAN+Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      TRUE ~ NA_real_
    ))

#ARA flux OVERALL
summary_araflux <- flux %>%
  group_by(treat) %>%
  summarise(
    mean_value = mean(ARAflux),
    sd_value = sd(ARAflux),
    count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lower_ci = mean_value - qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    upper_ci = mean_value + qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    effect = case_when(
      treat == "ALAN" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "ALAN+Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      TRUE ~ NA_real_
    ))

#EPA flux
summary_epaflux <- flux %>%
  group_by(treat) %>%
  summarise(
    mean_value = mean(EPAflux),
    sd_value = sd(EPAflux),
    count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lower_ci = mean_value - qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    upper_ci = mean_value + qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    effect = case_when(
      treat == "ALAN" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "ALAN+Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      TRUE ~ NA_real_
    ))

#FA fluxes per week despite treatment
#total FA flux
allweeks_faflux <- flux %>%
  group_by(week) %>%
  summarise(
    mean_value = mean(totalFAflux),
    sd_value = sd(totalFAflux),
    count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lower_ci = mean_value - qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    upper_ci = mean_value + qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    effect = case_when(
      week == "Week4" ~ ((mean_value[which(week == "Week1")] - mean_value) / mean_value[which(week == "Week1")]) * 100,
      week == "Week6" ~ ((mean_value[which(week == "Week1")] - mean_value) / mean_value[which(week == "Week1")]) * 100,
      TRUE ~ NA_real_
    ))

#PUFA flux
allweeks_pufaflux <- flux %>%
  group_by(week) %>%
  summarise(
    mean_value = mean(PUFAflux),
    sd_value = sd(PUFAflux),
    count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lower_ci = mean_value - qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    upper_ci = mean_value + qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    effect = case_when(
      week == "Week4" ~ ((mean_value[which(week == "Week1")] - mean_value) / mean_value[which(week == "Week1")]) * 100,
      week == "Week6" ~ ((mean_value[which(week == "Week1")] - mean_value) / mean_value[which(week == "Week1")]) * 100,
      TRUE ~ NA_real_
    ))

## 3.4 Plots ----
## FA fluxes 

#set factors
chiro$treat<-factor(chiro$treat, levels = c("Control", "ALAN", "Crayfish", "ALAN+Crayfish"))
chiro$week <- factor(chiro$week, levels = c("Week1", "Week4", "Week6"))

plot_fa <-ggplot(flux, aes(x = treat, y = totalFAflux))+
  geom_jitter(aes(color = treat), width = 0.2, alpha = 0.6, size = 3)+
  stat_summary(fun = median, geom = "crossbar", width = 0.4, fatten = 3, color = "black") +
  stat_summary(fun.min = min, fun.max = max, fun = median,
               geom = "errorbar", width = 0.2, color = "black")+
  scale_color_manual(values = c("#303030", "#FFB90F","#009ACD", "#5D478B"))+
  ylab(bquote("total FA flux" ~ (mu * g ~ m^{-2} ~ d^{-1}))) + 
  facet_wrap(~week) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size = 14, face= "bold"),
        legend.position = "none")  
plot_fa

plot_sfa <-ggplot(flux, aes(x = treat, y = SFAflux))+
  geom_jitter(aes(color = treat), width = 0.2, alpha = 0.6, size = 3)+
  stat_summary(fun = median, geom = "crossbar", width = 0.4, fatten = 3, color = "black") +
  stat_summary(fun.min = min, fun.max = max, fun = median,
               geom = "errorbar", width = 0.2, color = "black")+
  scale_color_manual(values = c("#303030", "#FFB90F","#009ACD", "#5D478B"))+
  ylab(bquote("SFA flux" ~ (mu * g ~ m^{-2} ~ d^{-1}))) + 
  facet_wrap(~week) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14),
        strip.text = element_blank(),
        legend.position = "none")  
plot_sfa

plot_mufa <-ggplot(flux, aes(x = treat, y = MUFAflux))+
  geom_jitter(aes(color = treat), width = 0.2, alpha = 0.6, size = 3)+
  stat_summary(fun = median, geom = "crossbar", width = 0.4, fatten = 3, color = "black") +
  stat_summary(fun.min = min, fun.max = max, fun = median,
               geom = "errorbar", width = 0.2, color = "black")+
  scale_color_manual(values = c("#303030", "#FFB90F","#009ACD", "#5D478B"))+
  ylab(bquote("MUFA flux" ~ (mu * g ~ m^{-2} ~ d^{-1}))) + 
  facet_wrap(~week) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14),
        strip.text = element_blank(),
        legend.position = "none")  
plot_mufa

plot_pufa <-ggplot(flux, aes(x = treat, y = PUFAflux))+
  geom_jitter(aes(color = treat), width = 0.2, alpha = 0.6, size = 3)+
  stat_summary(fun = median, geom = "crossbar", width = 0.4, fatten = 3, color = "black") +
  stat_summary(fun.min = min, fun.max = max, fun = median,
               geom = "errorbar", width = 0.2, color = "black")+
  scale_color_manual(values = c("#303030", "#FFB90F","#009ACD", "#5D478B"))+
  ylab(bquote("PUFA flux" ~ (mu * g ~ m^{-2} ~ d^{-1}))) + 
  facet_wrap(~week) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14),
        strip.text = element_blank(),
        legend.position = "none")  
plot_pufa

plot_ara <-ggplot(flux, aes(x = treat, y = ARAflux))+
  geom_jitter(aes(color = treat), width = 0.2, alpha = 0.6, size = 3)+
  stat_summary(fun = median, geom = "crossbar", width = 0.4, fatten = 3, color = "black") +
  stat_summary(fun.min = min, fun.max = max, fun = median,
               geom = "errorbar", width = 0.2, color = "black")+
  scale_color_manual(values = c("#303030", "#FFB90F","#009ACD", "#5D478B"))+
  ylab(bquote("ARA flux" ~ (mu * g ~ m^{-2} ~ d^{-1}))) + 
  facet_wrap(~week) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 14),
        strip.text = element_blank(),
        legend.position = "none")  
plot_ara

plot_epa <-ggplot(flux, aes(x = treat, y = EPAflux))+
  geom_jitter(aes(color = treat), width = 0.2, alpha = 0.6, size = 3)+
  scale_y_log10()+
  stat_summary(fun = median, geom = "crossbar", width = 0.4, fatten = 3, color = "black") +
  stat_summary(fun.min = min, fun.max = max, fun = median,
               geom = "errorbar", width = 0.2, color = "black")+
  scale_color_manual(values = c("#303030", "#FFB90F","#009ACD", "#5D478B"))+
  ylab(bquote("log EPA flux" ~ (mu * g ~ m^{-2} ~ d^{-1}))) + 
  facet_wrap(~week) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 14),
        strip.text = element_blank(),
        legend.position = "none") + 
  scale_x_discrete(labels = c("Control", "ALAN", "Crayfish", "ALAN+Crayfish"))  
plot_epa

#Prepare Figure 3
jpeg('Figure3b.jpeg', width =8, height = 12, units = 'in', res=300)
figure3b <- grid.arrange(
  plot_fa, plot_pufa, plot_ara, plot_epa, ncol = 1)
figure3b
dev.off()

# 4. FA concentrations of spiders ----
## 4.1. Data import and manipulation ----
spider_all <- read.csv("FA_spiders.csv", h=T, sep=",")

#remove negative values (due to blanks)
spider_all[spider_all < 0] <-0 

spider <- spider_all[,-1] #delete first column

#set factors
spider$treat<- factor(spider_all$treat, levels = c("Control", "ALAN", "Crayfish", "ALAN+Crayfish"))

#reorganize dataset for calculations by grouping FAs in categories
spider$SFA <- apply(spider[,4:15], 1, sum) #saturated FA
spider$MUFA <- apply(spider[,16:23], 1, sum) #monoinsaturated FA
spider$PUFA <- apply(spider[24:33], 1, sum) #poliunsaturated FA
spider$totalFA <- apply(spider[,34:36], 1, sum) #total FA

## 4.2. Models ----
#set dummy variables
spider$ALAN <- spider$treat
spider$Crayfish <- spider$treat

spider$ALAN <- plyr::mapvalues(spider$ALAN, from=c("ALAN+Crayfish", "Crayfish", "ALAN", "Control"), to=c(1,0,1,0))
spider$Crayfish <- plyr::mapvalues(spider$Crayfish, from=c("ALAN+Crayfish", "Crayfish", "ALAN", "Control"), to=c(1,1,0,0))

#anova type=2 for unbalaced design (Langsrud, 2003)

#totalFA
#explore distribution  
hist(spider$totalFA)

#homogeneity of variance
leveneTest(totalFA ~ treat, data = spider)

#linear model
lm_fa<- lm(totalFA ~ALAN*Crayfish, data=spider)

#check model parameters
res <- residuals(lm_fa, Uso="normalized")
qqnorm(res); qqline(res)                              
plot(fitted(lm_fa), res); abline(h = 0, lty = 2)       
vif(lm_fa, type="predictor") #check variance inflation factor

#ANOVA
summary(lm_fa)
Anova (lm_fa, type=2) #alan=0.09

#SFA
#explore distribution  
hist(spider$SFA)

#homogeneity of variance
leveneTest(SFA ~ treat, data = spider)

#linear model
lm_sa<- lm(SFA~ALAN*Crayfish, data=spider)

#check model parameters
res <- residuals(lm_sa, Uso="normalized")
qqnorm(res); qqline(res)                              
plot(fitted(lm_fa), res); abline(h = 0, lty = 2)
vif(lm_sa, type="predictor")

#ANOVA
summary(lm_sa)
Anova(lm_sa, type=2) 

#MUFA
#explore distribution 
hist(spider$MUFA)

#homogeneity of variance
leveneTest(MUFA ~ treat, data = spider)

#linear model
lm_mu<- lm(MUFA~ALAN*Crayfish, data=spider)

#check model parameters
Res <- residuals(lm_mu, Uso="normalized")
qqnorm(res); qqline(res)                              
plot(fitted(lm_mu), res); abline(h = 0, lty = 2)
vif(lm_sa, type="predictor")

#ANOVA
summary(lm_mu)
Anova(lm_mu, type=2) #alan=0.081

# PUFA
#explore distribution 
hist(spider$PUFA)

#homogeneity of variance
leveneTest(PUFA ~ treat, data = spider)

#linear model
lm_pu<- lm(PUFA~ALAN*Crayfish, data=spider)

#check model parameters
res <- residuals(lm_pu, Uso="normalized")
qqnorm(res); qqline(res)                              
plot(fitted(lm_pu), res); abline(h = 0, lty = 2)
vif(lm_pu, type="predictor")

#ANOVA
summary(lm_pu)
Anova(lm_pu, type=2) #alan=0.078

#Physiologically important PUFA of spiders ----

#ARA
#explore distribution 
hist(spider$C20.4_ARA)
hist(sqrt(spider$C20.4_ARA))

#homogeneity of variance
leveneTest(C20.4_ARA ~ treat, data = spider)
leveneTest(sqrt(C20.4_ARA) ~ treat, data = spider)

#linear model
lm_ara<- lm((sqrt(C20.4_ARA)) ~ ALAN*Crayfish, data=spider)

#check model parameters
res <- residuals(lm_ara, Uso="normalized")
qqnorm(res); qqline(res)                              
plot(fitted(lm_ara), res); abline(h = 0, lty = 2)
vif(lm_ara, type="predictor")

#ANOVA
summary(lm_ara)
Anova(lm_ara, type=2) #alan=0.058

#EPA
#explore distribution 
hist(spider$C20.5_EPA)

#homogeneity of variance
leveneTest(C20.5_EPA ~ treat, data = spider)

#linear model
lm_epa<- lm(spider$C20.5_EPA~ALAN*Crayfish, data=spider)

#check model parameters
res <- residuals(lm_epa, Uso="normalized")
qqnorm(res); qqline(res)                              
plot(fitted(lm_epa), res); abline(h = 0, lty = 2)
vif(lm_epa, type="predictor")

#ANOVA
summary(lm_epa)
Anova(lm_epa, type=2) #alan p=0.094

## 4.3. Summary ----
#FA content
sum_fa <- spider %>%
  group_by(treat) %>%
  summarise(
    mean_value = mean(totalFA),
    sd_value = sd(totalFA),
    count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lower_ci = mean_value - qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    upper_ci = mean_value + qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    effect = case_when(
      treat == "ALAN" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "ALAN+Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      TRUE ~ NA_real_
    ))

#labels for including effect sizes in plots
sum_fa$efsize <- ifelse(sum_fa$effect > 0,
                        paste0("+", round(sum_fa$effect, 0), "%"),
                        paste0(round(sum_fa$effect, 0), "%"))

#PUFA content
sum_pufa <- spider %>%
  group_by(treat) %>%
  summarise(
    mean_value = mean(PUFA),
    sd_value = sd(PUFA),
    count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lower_ci = mean_value - qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    upper_ci = mean_value + qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    effect = case_when(
      treat == "ALAN" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "ALAN+Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      TRUE ~ NA_real_
    ))

#labels for including effect sizes in plots
sum_pufa$efsize <- ifelse(sum_pufa$effect > 0,
                        paste0("+", round(sum_pufa$effect, 0), "%"),
                        paste0(round(sum_pufa$effect, 0), "%"))

#ARA content
sum_ara <- spider %>%
  group_by(treat) %>%
  summarise(
    mean_value = mean(C20.4_ARA),
    sd_value = sd(C20.4_ARA),
    count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lower_ci = mean_value - qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    upper_ci = mean_value + qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    effect = case_when(
      treat == "ALAN" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "ALAN+Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      TRUE ~ NA_real_
    ))

#labels for including effect sizes in plots
sum_ara$efsize <- ifelse(sum_ara$effect > 0,
                          paste0("+", round(sum_ara$effect, 0), "%"),
                          paste0(round(sum_ara$effect, 0), "%"))

#EPA content
sum_epa <- spider %>%
  group_by(treat) %>%
  summarise(
    mean_value = mean(C20.5_EPA),
    sd_value = sd(C20.5_EPA),
    count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    lower_ci = mean_value - qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    upper_ci = mean_value + qt(0.975, df = count - 1) * (sd_value / sqrt(count)),
    effect = case_when(
      treat == "ALAN" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      treat == "ALAN+Crayfish" ~ ((mean_value - mean_value[which(treat == "Control")]) / mean_value[which(treat == "Control")]) * 100,
      TRUE ~ NA_real_
    ))

#labels for including effect sizes in plots
sum_epa$efsize <- ifelse(sum_epa$effect > 0,
                         paste0("+", round(sum_epa$effect, 0), "%"),
                         paste0(round(sum_epa$effect, 0), "%"))

#summary for addition of effect sizes to figures
sum_fa2<- sum_fa[,-c(2:7)]
sum_pufa2<- sum_pufa[,-c(2:7)]
sum_ara2 <- sum_ara[,-c(2:7)]
sum_epa2 <- sum_epa[,-c(2:7)]

## 4.4 Plots ----
spider$treat<- factor(spider_all$treat, levels = c("Control", "ALAN", "Crayfish", "ALAN+Crayfish"))

spiderplot1 <-ggplot(spider, aes(x = treat, y = totalFA))+
  geom_jitter(aes(color = treat), width = 0.2, alpha = 0.6, size = 3)+
  stat_summary(fun = median, geom = "crossbar", width = 0.4, fatten = 3, color = "black") +
  stat_summary(fun.min = min, fun.max = max, fun = median,
               geom = "errorbar", width = 0.2, color = "black")+
  scale_color_manual(values = c("#303030", "#FFB90F","#009ACD", "#5D478B"))+
  ylab(bquote("total FA" ~ (mu * g ~ mg^{-1} ))) + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(size = 12, face= "bold"),
        legend.position = "none")  +
  geom_text(data = sum_fa2, aes(x = treat, y = max(spider$totalFA) * 1.05, label = efsize), 
            size = 4, fontface = "bold")
spiderplot1
#missing values correspond to the effect size value from the control. 

spiderplot2 <-ggplot(spider, aes(x = treat, y = PUFA))+
  geom_jitter(aes(color = treat), width = 0.2, alpha = 0.6, size = 3)+
  stat_summary(fun = median, geom = "crossbar", width = 0.4, fatten = 3, color = "black") +
  stat_summary(fun.min = min, fun.max = max, fun = median,
               geom = "errorbar", width = 0.2, color = "black")+
  scale_color_manual(values = c("#303030", "#FFB90F","#009ACD", "#5D478B"))+
  ylab(bquote("PUFA" ~ (mu * g ~ mg^{-1} ))) + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(size = 12, face= "bold"),
        legend.position = "none")  +
  geom_text(data = sum_pufa2, aes(x = treat, y = max(spider$PUFA) * 1.05, label = efsize), 
            size = 4, fontface = "bold")
spiderplot2

spiderplot3 <-ggplot(spider, aes(x = treat, y = C20.4_ARA))+
  geom_jitter(aes(color = treat), width = 0.2, alpha = 0.6, size = 3)+
  stat_summary(fun = median, geom = "crossbar", width = 0.4, fatten = 3, color = "black") +
  stat_summary(fun.min = min, fun.max = max, fun = median,
               geom = "errorbar", width = 0.2, color = "black")+
  scale_color_manual(values = c("#303030", "#FFB90F","#009ACD", "#5D478B"))+
  ylab(bquote("ARA" ~ (mu * g ~ mg^{-1} ))) + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(size = 12, face= "bold"),
        legend.position = "none")  +
  geom_text(data = sum_ara2, aes(x = treat, y = max(spider$C20.4_ARA) * 1.05, label = efsize), 
            size = 4, fontface = "bold")
spiderplot3

spiderplot4 <-ggplot(spider, aes(x = treat, y = C20.5_EPA))+
  geom_jitter(aes(color = treat), width = 0.2, alpha = 0.6, size = 3)+
  stat_summary(fun = median, geom = "crossbar", width = 0.4, fatten = 3, color = "black") +
  stat_summary(fun.min = min, fun.max = max, fun = median,
               geom = "errorbar", width = 0.2, color = "black")+
  scale_color_manual(values = c("#303030", "#FFB90F","#009ACD", "#5D478B"))+
  ylab(bquote("EPA" ~ (mu * g ~ mg^{-1} ))) + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(size = 12, face= "bold"),
        legend.position = "none")  +
  geom_text(data = sum_epa2, aes(x = treat, y = max(spider$C20.5_EPA) * 1.05, label = efsize), 
            size = 4, fontface = "bold")
spiderplot4

#spiders PUFA plot
jpeg('Figure4.jpeg', width = 6, height = 9, units = 'in', res = 300)
spiderplot <- grid.arrange(
  spiderplot1, spiderplot2, spiderplot3,spiderplot4, ncol=1)
dev.off()

# 5. Supplementary Material ----

## 5.1. Signal crayfish activity under ALAN ----

#data
cray <- read.csv("crayfish.csv", h=T, sep=";")

str(cray)

# set factors
cray$Week <- as.numeric(cray$Week)
cray$Treat <- as.factor(cray$Treat)
cray$Treat <- factor(cray$Treat, levels = c("NoALAN", "ALAN"))
#cray$Flume <- as.factor(cray$Flume)

#Mixed Model with Binomial distribution
lme_cray <- glmer((cbind(Active, Hidden)) ~ Treat + Week + (1|Flume), 
                  family=binomial, data = cray)

summary(lme_cray)

#model check
text(qqnorm(resid(lme_cray)))
qqline(resid(lme_cray))
check_overdispersion(lme_cray) 
check_collinearity(lme_cray)

check_singularity(lme_cray) #TRUE
summary(lme_cray)$varcor #variance of "flume" = 0 
#random effect structure is not supported by the data
#change to generalised linear model

# Generalised Linear Model
glm_cray <- glm(cbind(Active, Hidden) ~ Treat + Week,
           family = binomial, data = cray)

glm_cray0 <- glm(cbind(Active, Hidden) ~ 1,
                family = binomial, data = cray)

AIC(lme_cray,glm_cray, glm_cray0) #glm_cray lowest AIC

#Model results
summary(glm_cray)

# ANOVA
Anova(glm_cray, type = "II", test.statistic = "LR") #Treat: p=0.035

#model diagnostics
sim<- simulateResiduals(fittedModel = glm_cray, plot = F)
plot(sim)

#Summary
cray_sum <- cray %>%
  group_by(Week, Treat) %>%
  summarise(
    Active = sum(Active),
    Hidden = sum(Hidden),
    .groups = "drop"
  ) %>%
  mutate(
    Total = Active + Hidden,
    Active_prop = (Active / Total) * 100,
    Hidden_prop = (Hidden / Total) * 100
  ) %>%
  pivot_longer(
    cols = c(Active_prop, Hidden_prop),
    names_to = "Status",
    values_to = "Percentage"
  ) %>%
  group_by(Treat, Status) %>%
  summarise(
    mean = mean(Percentage),
    se = sd(Percentage) / sqrt(n()),
    .groups = "drop"
  )

#Plot
jpeg('Crayfish.jpeg', width = 10, height = 6, units = 'in', res=300)
cray_plot <- ggplot(cray_sum, aes(x = Status, y = mean, fill = Treat)) +
  geom_col(position = position_dodge(0.9), width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(0.9), width = 0.2) + 
  scale_fill_manual(values = c("#009ACD", "#FFB90F")) +
  scale_x_discrete(labels = c("Active", "Hidden")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) + 
  labs(
    x = NULL,
    y = "Percentage (%)",
    fill = "Treatment"
  ) +
  theme_minimal(base_size = 14) + 
  theme(
    axis.text.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15)
  ) +
  annotate("text", x = 1.5, y = 105, label = "GLM: p = 0.035", size = 5) +
  coord_cartesian(ylim = c(0, 110)) 
cray_plot
dev.off()

## 5.2. FA content in chironomids ----

#Calculation of mean concentrations and standard deviation of FA groups 
#content in chironomids. The obtained data are presented in Table S5. 

##security copy from L.269

#subsets per time point
sm4_w1 <- subset(sm4, sm4$week=="Week1")
sm4_w4 <- subset(sm4, sm4$week=="Week4")
sm4_w6 <- subset(sm4, sm4$week=="Week6")

data_list <- list(w1 = sm4_w1, w4 = sm4_w4, w6 = sm4_w6)

response_vars <- c("totalFA", "SFA", "MUFA", "PUFA", "C20.4_ARA", "C20.5_EPA")

#function for summary calculations
summary_wide <- function(data, response_var, time_name) {
  data %>%
    group_by(treat) %>%
    summarise(
      mean = mean(.data[[response_var]], na.rm = TRUE),
      sd = sd(.data[[response_var]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(mean, sd), names_to = "statistic", values_to = "value") %>%
    pivot_wider(names_from = treat, values_from = value) %>%
    mutate(time = time_name, .before = 1) 
}

results <- list()

for (time_name in names(data_list)) {
  for (var in response_vars) {
    df <- data_list[[time_name]]
    result_name <- paste0("sum_", var, "_", time_name)
    results[[result_name]] <- summary_wide(df, var, time_name)  
  }
}

final_summary_fa <-bind_rows(results, .id = "source")

## 5.3. FA content in spiders ----

#Calculation of mean concentrations and standard deviations of FA groups 
#content in spiders. The obtained data is presented as Supplementary Table 7. 

#total FA 
# mean and standard deviation for each treatment
summary_fa_sp <- spider %>%
  group_by(treat) %>%
  summarise(mean_value = mean(totalFA),
            sd_value = sd(totalFA),
            count = n())

#SAFA 
# mean and standard deviation for each treatment
summary_safa_sp <- spider %>%
  group_by(treat) %>%
  summarise(mean_value = mean(SFA),
            sd_value = sd(SFA),
            count = n())

#MUFA 
# mean and standard deviation for each treatment
summary_mufa_sp <- spider %>%
  group_by(treat) %>%
  summarise(mean_value = mean(MUFA),
            sd_value = sd(MUFA),
            count = n())

#PUFA 
# mean and standard deviation for each treatment
summary_pufa_sp <- spider %>%
  group_by(treat) %>%
  summarise(mean_value = mean(PUFA),
            sd_value = sd(PUFA),
            count = n())

#ARA
# mean and standard deviation for each treatment
summary_ara_sp <- spider %>%
  group_by(treat) %>%
  summarise(mean_value = mean(C20.4_ARA),
            sd_value = sd(C20.4_ARA),
            count = n())

#EPA
# mean and standard deviation for each treatment
summary_epa_sp <- spider %>%
  group_by(treat) %>%
  summarise(mean_value = mean(C20.5_EPA),
            sd_value = sd(C20.5_EPA),
            count = n())

###

