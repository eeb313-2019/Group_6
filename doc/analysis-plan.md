# Possible analyses
### I. How to answer research questions 
community composition: identify dominant genus, abundance/importance (if have over 10 genus)

#### 1. community composition with changing seasons
community composition 
season
sites

years 
linear mixed model - year random effect 

#### 2. community composition will be the same across sites with the same water-body type, and different acorss sites with different water-body types
water body type
community composition 
size of species 
sites

PCA 

#### 3. Changing pH levels will affect species community composition and abundances.
pH
community composition 
abundances/richness 
similar are species compositions in different treatments, you will need ordination methods (RDA, ANOSIM, NMDS)


#### II. How do we formulate the data

how you anticipate getting from the raw data to whatever summary data you will use to generate a given plot 
(i.e. explain a data **cleaning/transformation plan**)

Merge data from each year and season (years)
individual count to relative abundance (across sites)
retain only genus information 

For question 1: 
Group data by year and season (year - random effect)

For question 2: 
Group by water body type for sites 

For question 3:
Group by pH (site - random effect or spatial autocorrelated)

#### III. what kinds of statistical approaches you anticipate employing
principal component analysis 
PCA - accross season over average year 
compare 3 PCA 
multivariant model, 3 season accross sites 
spatial autocorrelation test 
ANOVA 
PERMANOVA 


# Possilbe results tables
Anova and linear mixed model 


# Possible results figures 
variogram 
make a map, what sites for conservation, what is healthy, rare species vulnerable 





