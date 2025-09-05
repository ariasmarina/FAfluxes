# FA fluxes - Overview
This is the accompanying research data and code for the paper "Multiple stressors across ecosystem boundaries: do invasive species and light pollution change the quality of aquatic subsidies to terrestrial predators?" by Arias et al. (submitted to ----)

This study addresses the individual and joint effect of Artificial Light at Night (ALAN) and the invasive signal crayfish (Pacifastacus leniusculus) on the fatty acid (FA) fluxes from aquatic to terrestrial ecosystems in a fully crossed experimental design employed at the Riparian Stream Mesocosm (RSM) facility in Landau, Germany. 

The FA profile of emergent aquatic insects, represented by chironomids, and riparian predators, represented by the web-building spider Tetragnatha extensa L., was assessed, as shifts in the FA fluxes from stream to land are likely mirrored in the spiders that depend strongly on aquatic resources.

Four datasets are provided, and two new datasets are created.

#1. Biomass of emergent aquatic insects 
### emergence_biomass.csv
This dataset contains the biomass of each sample in mg (chiro_bio) collected from three emergence traps in 16 flumes, twice a week, after one, four, and six weeks of stressor exposure (N=288).

treat: treatments; "Control" (no stressors); "Crayfish" (presence of signal crayfish, absence of ALAN); "ALAN" (presence of ALAN, absence of signal crayfish); "ALAN+Crayfish" (presence of both ALAN and signal crayfish).
week: time after initiation of stressors exposure measured in weeks= one week "Week1", four weeks "Week4", and six weeks "Week6".
day: sampling was done twice a week, day 1 (d1) and day 2 (d2).
flume: number of each artificial stream of the experimental system (1 to 16). 
trap: emergence traps (n=3) located in the upper  ("up"), middle ("mid") and downstream ("down") reaches.
chiro_bio: biomass (mg) of emerged non-tanypodinae chironomids in each trap

The code describes the calculation of emergence biomass rate, which will be used later for the calculation of FA fluxes. 

First, emergent biomass was summed among the two sampling days to obtain the total emergence biomass over a week.
Then, emergent biomass was standardised to 1m2 considering the width of the flume in each section 
i.e. upper reach (trap "up"): 0.7cm
i.e. middle reach (trap "mid"): 0.75 cm 
i.e. downstream reach (trap "down") 0.8 cm

The total biomass per week was then converted to rate (mg/m2/d). Finally, data was averaged among the three traps, obtaining the mean emergence biomass (mg/m2/d) (n=48). 
A new dataset is created, entitled "mean_biomass.csv".

### mean_biomass.csv
x: number of observation
treat: treatments (same as above)
week: time after initiation of stressor exposure (in weeks, same as above)
flume: number of each artificial stream of the experimental system (1 to 16). 
mean.bio: mean emergence biomass (mg/m2/d)

#2. Fatty acids composition of chironomids
### FA_emergence.csv

This dataset contains the FA composition of chironomid samples (n=10) of the 16 streams of the RSM (four per treatment) after 1, 4 and 6 weeks of exposure to stressors.

week: time after initiation of stressors exposure measured in weeks= one week "Week1", four weeks "Week4", and six weeks "Week6". 
treat: treatments; "Control" (no stressors); "Crayfish" (presence of signal crayfish, absence of ALAN); "ALAN" (presence of ALAN, absence of signal crayfish); "ALAN+Crayfish" (presence of both ALAN and signal crayfish).
flume: number of each artificial stream of the experimental system (1 to 16). 
C12.0: Lauric acid (µg mg-1)
C13.0: Tridecanoic acid (µg mg-1)
C14.0: Myristic acid (µg mg-1)
C15.0: Pentadecanoic acid (µg mg-1)
C16.0: Palmitic acid (µg mg-1)
C17.0: Heptadecanoic acid (µg mg-1)
C18.0: Stearic acid (µg mg-1)
C20.0: Arachidic acid (µg mg-1)
C21.0: Heneicosanoid acid (µg mg-1)
C22.0: Behenic acid (µg mg-1)
C23.0: Tricosanoic acid (µg mg-1)
C24.0: Lignoceric acid (µg mg-1)
C14.1: Myristoleic acid (µg mg-1)
C15.1: Pentadecenoic acid (µg mg-1)
C16.1: Palmitoleic acid (µg mg-1)
C17.1: Heptadecenoic acid (µg mg-1)
C18.1: Oleic acid (µg mg-1)
C20.1: Gadoleic acid (µg mg-1)
C22.1: Erucic acid (µg mg-1)
C24.1: Nervonic acid (µg mg-1)
C18.2: Linoleic acid (µg mg-1)
C18.3_GLA: y-linolenic acid (µg mg-1)
C18.3_ALA: a-linolenic acid (µg mg-1)
C20.2: Eicosadienoic acid (µg mg-1)
C20.3a Eicosatrienoic acid (cis-11,14,17) (µg mg-1)
C20.3b: Eicosatrienoic acid (cis-8,11,14) (µg mg-1)
C20.4_ARA: Arachidonic acid (µg mg-1)
C20.5_EPA: Eicosapentaenoic acid (µg mg-1)
C22.2: Docosadienoic acid (µg mg-1)
C22.6_DHA: Docosahexaenoic acid (µg mg-1)

The code describes the calculation main FA groups concentration, which will be used later for the calculation of FA fluxes. Main FA groups concentration was multiplied to the obtained mean emergence biomass (mg m-2 d-1) for FA fluxes calculation, study object of the present study. 
A new dataset is created, entitled "fluxes.csv".

### fluxes.csv
x: number of observation
week: time after initiation of stressors exposure measured in weeks (same as above) 
treat: treatments "Control" (no stressors); "Crayfish" (presence of signal crayfish, absence of ALAN); "ALAN" (presence of ALAN, absence of signal crayfish); "ALAN+Crayfish" (presence of both ALAN and signal crayfish).
flume: number of each artificial stream of the experimental system (1 to 16).
SFAflux: saturated FA flux (µg m-2 d-1)
MUFAflux: monounsaturated FA flux (µg m-2 d-1)
PUFAflux: polyunsaturated FA flux (µg m-2 d-1)
totalFAflux: total FA flux (µg m-2 d-1)
ARAflux: ARA flux (µg m-2 d-1)
EPAflux: EPA flux (µg m-2 d-1)

#3. Fatty acid composition of Riparian Spiders
This dataset contains the FA composition of samples of Tetragnata extensa (n=1) of the 16 streams of the RSM (four per treatment) after 6 weeks of exposure to stressors.

### FA_spiders.csv
x: number of observation
treat: "Control" (no stressors); "Crayfish" (presence of signal crayfish, absence of ALAN); "ALAN" (presence of ALAN, absence of signal crayfish); "ALAN+Crayfish" (presence of both ALAN and signal crayfish).
flume: number of each artificial stream of the experimental system (1 to 16).
sample: sample name
C12.0: Lauric acid (µg mg-1)
C13.0: Tridecanoic acid (µg mg-1)
C14.0: Myristic acid (µg mg-1)
C15.0: Pentadecanoic acid (µg mg-1)
C16.0: Palmitic acid (µg mg-1)
C17.0: Heptadecanoic acid (µg mg-1)
C18.0: Stearic acid (µg mg-1)
C20.0: Arachidic acid (µg mg-1)
C21.0: Heneicosanoid acid (µg mg-1)
C22.0: Behenic acid (µg mg-1)
C23.0: Tricosanoic acid (µg mg-1)
C24.0: Lignoceric acid (µg mg-1 )
C14.1: Myristoleic acid (µg mg-1)
C15.1: Pentadecenoic acid (µg mg-1)
C16.1: Palmitoleic acid (µg mg-1)
C17.1: Heptadecenoic acid (µg mg-1)
C18.1: Oleic acid (µg mg-1)
C20.1: Gadoleic acid (µg mg-1)
C22.1: Erucic acid (µg mg-1)
C24.1: Nervonic acid (µg mg-1)
C18.2: Linoleic acid (µg mg-1)
C18.3_GLA: y-linolenic acid (µg mg-1)
C18.3_ALA: a-linolenic acid (µg mg-1)
C20.2: Eicosadienoic acid (µg mg-1)
C20.3a Eicosatrienoic acid (cis-11,14,17) (µg mg-1)
C20.3b: Eicosatrienoic acid (cis-8,11,14) (µg mg-1)
C20.4_ARA: Arachidonic acid (µg mg-1)
C20.5_EPA: Eicosapentaenoic acid (µg mg-1)
C22.2: Docosadienoic acid (µg mg-1)
C22.6_DHA: Docosahexaenoic acid (µg mg-1)

#4. Crayfish activity

### crayfish.csv
This dataset contains the weekly records of signal crayfish activity in flumes exposed and non-exposed to ALAN. Individuals were considered “Hidden” when they remain inside the hiding places or in burrows in the sediment, while those visible outside hiding places or walking were considered “Active”. 

Week: time after initiation of stressors exposure measured in weeks (1, 2, 3, 4, 5, 6)
Date: date of night-check
Treat: presence (ALAN) or absence (NoALAN) or artificial light at night. 
Flume: number of each artificial stream of the experimental system (1 to 16). 
Hidden: number of individuals of signal crayfish inside the hiding places or in burrows in the sediment
Active: number of individuals f signal crayfish visible outside hiding places or walking
Total: Active + Hidden
