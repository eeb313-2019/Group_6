# Analyses plan
 
1)	Measuring pH variability across sites and seasons
-	The effect of site, year, and season on pH:
-	Simple boxplots by sites and seasons: 3 for each site 
-	**somewhere in lecture** there is code that tests for 4 assumptions of normality – do that before running linear model 
-	Run multiple models with different structures: run AIC for the best structure 
-	Account for spatial autocorrelation 
-	Identify which variable to include as random effect (whatever the fixed effect is) 

2)	Measuring conductance variability across sites and seasons
-	The effect of site, year, and season on conductance: 
-	Simple boxplots by sites and seasons: 3 for each site 
-	**somewhere in lecture** there is code that tests for 4 assumptions of normality – do that before running linear model 
-	Run multiple models with different structures: run AIC for the best structure 
-	Account for spatial autocorrelation 
-	Identify which variable to include as random effect (whatever the fixed effect is) 

3)	Create dataset for % Riparian Cover 
-	Change Y/N for landuse to 0 and 1 and get mean = call this column “% landuse” 
-	Only keep: mean riparian cover for each site and this new column

4)	Measure correlation between %Riparian Cover and other variables 
-	Measure correlation between riparian cover, land use, seasonality and site 
-	Account for spatial autocorrelation 
-	Identify which variable to include as random effect (whatever the fixed effect is) 

5)	Measuring Simpson’s Index across sites 
-	Convert species data to a matrix 
-	Spatial autocorrelation 
-	Use Simpson’s function to quantify diversity

6)	Either linear model or mixed model to measure effect of pH, conductance, riperian cover, and land use on Simpson's index
-	Need to account for spatial autocorrelation from step 3 
-	Model structure will depend on results from step 1 

8)	Measuring Shannon’s Index across sites and 
-	Convert species data to a matrix 
-	Spatial autocorrelation 
-	Use Simpson’s function to quantify richness 

9)	Either linear model or mixed model to measure effect of pH, conductance, riperian cover, and land use on Species evenness  
-	Need to account for spatial autocorrelation from step 3 
-	Model structure will depend on results from step 1 




