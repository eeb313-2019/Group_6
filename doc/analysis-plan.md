# Possible analyses
### I. How to answer research questions 
For variable community composition, we plan to plot them with PCA and spatial autocorrelation, identify dominant genus, abundance/importance (if have over 10 genus).
Including species richness as another variable related to community composition. 
Other variables includes: land use, pH, season, sites, year, water-body type. 
We predict that there is a direct relationship between land use, pH, and season to community composition, while land use and season also influence pH. Our model selection and PCA will be based on these predictions. 

#### Q1
	community composition/species richness
	season
PCA for community composition 
linear mixed model and model selection - year random effect 

#### Q2
	community composition/species richness
	land use (use proxy for land use)
Spatial correlation between community composition and land use
PCA 
Model selection with water-body type and sites as random effect

#### Q3
	pH
	community composition/species richness
Spatial correlation between community composition and pH
PCA 
Model selection with water-body type and sites as random effect


#### II. How do we formulate the data

**Cleaning/transformation plan: **

Merge data from each year and season (all years, all sites)
individual count to relative abundance (across sites)
retain only genus information 

For question 1: 
Group data by year and season (year - random effect)

For question 2: 
Group by land use proxy of sites 

For question 3:
Group by pH (site - random effect)

#### III. what kinds of statistical approaches you anticipate employing
Principal component analysis (PCA)
Linear mixed model and model selection  
spatial autocorrelation test 
ANOVA 
MANOVA/PERMANOVA


# Possilbe results tables
Anova/MANOVA/PERMANOVA
Linear mixed model 
Model selection 


# Possible results figures 
PCA
Variogram 
Make a map to show sites and factors (legend). 





