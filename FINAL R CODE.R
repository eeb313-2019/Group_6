## packages needed

```{r}
library(dplyr)
library(ggplot)
library (readr) 
library(vegan)
library(ggplot2)
library(tidyverse)
library(sp)
library(raster)
library(ape)
library(lawstat) 
library(rcompanion)
library(lme4)
library(MuMIn)
library(lawstat)
library(nlme)
library(lavaan)
library(PerformanceAnalytics)
```

## files needed
```{r}
invert_genera <- read_csv("invert.genera.season.csv")
FINAL.PH <- read_csv("FINAL.PH.GEO.csv") 
riparian_data <- read_csv("FINAL.RIP.GEO.csv", col_names = TRUE)
merged.final.file <- read_csv("merged.final.fle.csv")

```

##map of the sites
```{r}
map.US <- map_data("world", region = "USA")
invert_genera %>%
  ggplot() + 
  geom_polygon(data = map.US, aes(x = long, y = lat, group = group), fill = NA, colour = "gray") + 
  coord_fixed() + 
  geom_point(data = invert_genera, aes(x = decimalLongitude.x, y = decimalLatitude.x)) + 
  xlim(-200, -50) + 
  xlab("Longitude") + ylab("Latitude")
```

##inverse distance weight matrix for testing spatial autocorrelation

```{r}
coord_invert <- SpatialPoints(cbind(invert_evenness_coords$decimalLongitude, invert_evenness_coords$decimalLatitude), 
proj4string=CRS("+proj=longlat +ellps=WGS84"))
UTM.US <- spTransform(coord_invert, CRS("+init=epsg:2163")) # US National Atlas Equal Area
dist.matrix <- as.matrix(dist(data.frame(UTM.US)))
inv.dist<- 1/dist.matrix
diag(inv.dist) <- 0
# getting rid of infinite values
inv.dist[is.infinite(inv.dist)] <- 0
```

# 1)	Measuring pH variability across sites and seasons
-	The effect of site, year, and season on pH:
o	Simple boxplots by sites and seasons: 3 for each site 
o	there is code that tests for 4 assumptions of normality – do that before running linear model 
o	Run multiple models with different structures: run AIC for the best structure 
o	Account for spatial autocorrelation 
o	Identify which variable to include as random effect (whatever the fixed effect is) 

```{r}
#removing outliers that don't fit within the normal range of pH values (0-14) & removing TECR site
cleaned.final.ph <- FINAL.PH %>% 
  filter(pH < 14 & pH > 0) %>% 
  filter(!siteID == "TECR") %>% 
  dplyr::select(siteID, year, month, Season, pH, conductance, latitude, longitude) 

cleaned.final.ph
```

##boxplots to visualize variation in pH
```{r}
#by site
cleaned.final.ph %>% 
  ggplot(aes(x= siteID, y= pH))+ 
  geom_boxplot() + 
  labs(x="Site ID", y="pH") 

#by season
cleaned.final.ph %>% 
  ggplot(aes(x= Season, y= pH))+ 
  geom_boxplot() + 
  labs(x="Season", y="pH") 

#by seasons and sites
cleaned.final.ph %>% 
  ggplot(aes(x= Season, y= pH))+ 
  geom_boxplot() + 
  facet_wrap(~siteID) 
```

##assumptions
``` {r}
# Model that includes all fixed effects 
pH.lm<- lm(pH~Season*siteID, data= cleaned.final.ph)
par(mfrow=c(2,2))
plot(pH.lm)
#looking at panels 1 and 2, the residuals are homogenous and normally distributed 

# Independence of residuals 
lawstat::runs.test(pH.lm$residuals)
#p value is 0.08, i.e no autocorrelation 
```

# AIC vs. AICc
```{r}
cleaned.final.ph
n<- 74 
k<- 4 
AIC_mod <- n/k 
AIC_mod
#USING AICc because n/k is 18.5 (i.e. less than 40)
```

##linear models
```{r}
#Saturated Model 
pH.lm.sat.f <-lmer(pH~Season*siteID +(1|year/month), data= cleaned.final.ph)

#best random effect structure
pH.lm.sat.f<-lmer(pH~Season*siteID + (1|year/month), data= cleaned.final.ph, REML= TRUE)
pH.lm.sat<-lmer(pH~Season*siteID +(1|year), data= cleaned.final.ph, REML= TRUE)

MuMIn::AICc(pH.lm.sat, pH.lm.sat.f)
#pH.lm.sat.f has the lowest AICc value

#best fixed effect structure 
pH.lm.sat.f<-lmer(pH~Season*siteID + (1|year/month), data= cleaned.final.ph, REML= FALSE)
pH.lm.sat.season<-lmer(pH~Season +(1|year/month), data= cleaned.final.ph, REML= FALSE)
pH.lm.sat.siteID<-lmer(pH~siteID +(1|year/month), data= cleaned.final.ph, REML= FALSE)
pH.lm.sat.no.int<-lmer(pH~Season +siteID +(1|year/month) , data= cleaned.final.ph, REML= FALSE)
pH.lm.sat.ind<-lmer(pH~1  +(1|year/month), data= cleaned.final.ph, REML= FALSE)

MuMIn::AICc(pH.lm.sat.f, pH.lm.sat.season, pH.lm.sat.siteID, pH.lm.sat.no.int, pH.lm.sat.ind)
#pH.lm.sat.f has the lowest AICc value; best overall model

summary(pH.lm.sat.f)
r.squaredGLMM(pH.lm.sat.f)
```


##spatial autocorrelation 
```{r}
coord_invert <- SpatialPoints(cbind(cleaned.final.ph$longitude, cleaned.final.ph$latitude), 
proj4string=CRS("+proj=longlat +ellps=WGS84"))
UTM.US <- spTransform(coord_invert, CRS("+init=epsg:2163")) # US National Atlas Equal Area
dist.matrix <- as.matrix(dist(data.frame(UTM.US)))
inv.dist<- 1/dist.matrix
diag(inv.dist) <- 0
# getting rid of infinite values
inv.dist[is.infinite(inv.dist)] <- 0
Moran.I(cleaned.final.ph$pH, inv.dist, alternative = "two.sided") #p value is 0.59 so pH is not spatially autocorrelated 

```

# 2)	Measuring conductance variability across sites and seasons
-	The effect of site, year, and season on conductance: 
o	Simple boxplots by sites and seasons: 3 for each site 
o	test assumptions – do that before running linear model 
o	Run multiple models with different structures: run AIC for the best structure 
o	Account for spatial autocorrelation 
o	Identify which variable to include as random effect (whatever the fixed effect is) 

```{r}
#removed the TECR site since it didn't have any taxonomic data 
FINAL.PH <- read_csv("FINAL.PH.GEO.csv") 
df_ph <- FINAL.PH %>% 
  dplyr::select(siteID, year, month, Season, pH, conductance, latitude, longitude) %>% 
  filter(!siteID == "TECR")
```

##boxplots for each site & season
```{r}
#by season
df_ph_logged %>%
  ggplot(aes(x=Season, y=conductance)) + 
  geom_boxplot() + 
  labs(x="Season", y="Conductance (mS/cm)") 

#by site 
df_ph_logged %>%
  ggplot(aes(x=siteID, y=conductance)) + 
  geom_boxplot() +
  labs(x="Site ID", y="Conductance (mS/cm)")

#grouped by site, with each season plotted per site
df_ph_logged %>%
  ggplot(aes(x=Season, y=conductance)) + 
  geom_boxplot() + 
  facet_wrap(~siteID)
```

summary:
- there is not a big difference in conductance between seasons across all sites 
- there are big difference in the conductance between sites across all seasons 
- there doesn't seem to be a big difference in conductance between seasons within sites 
  + caution: not very many data points exist for each season for each individual site 

##assumptions:
```{r}
df_ph %>% 
  ggplot(aes(x=siteID, y=conductance)) + 
  geom_point() 

#remove two outliers due to sampling or inputting error 
df_ph_clean<- df_ph %>% 
  filter(conductance < 10000)
df_ph_clean

#ASSESSING NORMALITY:
##here is the plot for normality & homogeneity: 
basic.conductance.lm <- lmer(conductance ~ siteID*Season, data=df_ph_clean)
par(mfrow=c(2,2))
plot(basic.conductance.lm)

qqnorm(residuals(basic.conductance.lm))
plotNormalHistogram(basic.conductance.lm$residuals)

df_ph_logged <- df_ph_clean %>% 
  mutate(logged.conductance = log(conductance))

#ASSESSING HOMOGENEITY OF VARIANCE: 
basic.conductance.lm <- lm(logged.conductance ~ Season*siteID, data=df_ph_logged)
par(mfrow=c(2,2))
plot(basic.conductance.lm)
#panel 1 (no funnel shape) & 3 (approx horizontal line)
#although not great, the variance is much more equal than before 

#ASSESSING AUTOCORRELATION (IE. INDEPENDANCE OF RESIDUALS): 
lawstat::runs.test(basic.conductance.lm$residuals)
#since pvalue is not significant, that means no autocorrelation :) 

```
summary:
- removed two rows from the data that contained outliers
- plotted histograms for conductance values (all, by site, by season)
  + plotted log transformed conductance values (all, by season)
- created a new dataframe where conductance was log transformed 
- checked autocorrelationa & homogeneity of variance 


# AIC vs. AICc
```{r}
df_ph_logged
n<- 76
k<- 4 
AIC_mod <- n/k 
AIC_mod
#USING AICc because n/k is 19 (i.e. less than 40)
```

##linear models
```{r}
#saturated model
conductance.lm.sat <- lm(logged.conductance ~ Season*siteID + (1|year/month), data=df_ph_logged)

#best random effect structure 
conductance.lm.sat <- lmer(logged.conductance ~ Season*siteID + (1|year/month), data=df_ph_logged, REML=TRUE)
conductance.lm.year <- lmer(logged.conductance ~ Season*siteID + (1|year) , data=df_ph_logged, REML=TRUE)

AICc(conductance.lm.sat, conductance.lm.year)
#use conductance.lm.year since it has the lower AICc score

#best fixed effected structure 
conductance.lmer1 <- lmer(logged.conductance ~ Season*siteID + (1|year), data=df_ph_logged, REML=FALSE)
conductance.lmer2 <- lmer(logged.conductance ~ Season + (1|year), data=df_ph_logged, REML=FALSE)
conductance.lmer3 <- lmer(logged.conductance ~ siteID + (1|year), data=df_ph_logged,REML=FALSE)
conductance.lmer4 <- lmer(logged.conductance ~ 1 + (1|year), data=df_ph_logged,REML=FALSE)

AICc(conductance.lmer1, conductance.lmer2, conductance.lmer3, conductance.lmer4)
#best structure is: conductance.lmer3 

#checking normality of residuals 
qqnorm(residuals(conductance.lmer3))

summary(conductance.lmer3)
r.squaredGLMM(conductance.lmer3)
```
summary: 
- conductance.lmer3 model is the best fit, since lowest AICc score 
- only site ID predicts conductance
  + conductance value varies based on site 
  + site becomes random effect in final models predicting richness, evenness, etc.

##looking at spatial autocorrelation
```{r}
coord_invert <- SpatialPoints(cbind(df_ph_logged$longitude, df_ph_logged$latitude), proj4string = CRS("+proj=longlat +ellps=WGS84"))

UTM.US <- spTransform(coord_invert, CRS("+init=epsg:2163"))

dist_matrix <- as.matrix(dist(data.frame(UTM.US)))
inv_dist <- 1/dist_matrix
diag(inv_dist) <- 0

#get rid of infinite values 
inv_dist[is.infinite(inv_dist)] <- 0

#moran's I
Moran.I(df_ph_logged$logged.conductance, inv_dist, alternative = "two.sided") 
#pvalue is 2.09x10-6, so conductance is spatially autocorrelated 

#visualizing the null model using a variogram 
df_ph_logged <- bind_cols(df_ph_logged, as.data.frame(UTM.US))
df_ph_logged
conductance.null <- gls(logged.conductance ~ siteID, data=df_ph_logged, method="ML")
plot(Variogram(conductance.null, form = ~coords.x1+coords.x2, resType = "normalized"), main="Semivariogram including CARI site")
##it appears that one point is an outlier - the comparison is of locations that are very far from each other 
##we could try excluding the alaska site (cari) to see if it resolves the issue of spatial autocorrelation

```

##removoing site in alaska ("CARI") & rechecking spatial autocorrelation 
```{r}
#excluding mayf site in alaska 
conductance.cari.out <- df_ph_logged %>% 
  filter(!siteID == "CARI")
conductance.null.nocari <- gls(logged.conductance ~ siteID, data=conductance.cari.out, method="ML")
plot(Variogram(conductance.null.nocari, form = ~coords.x1+coords.x2, resType = "normalized"), main="Semivariogram excluding CARI site")
##the semivariogram now has a decreasing slope

#calculating moran's i without "cari" site
coord_invert.cari.out <- SpatialPoints(cbind(conductance.cari.out$longitude, conductance.cari.out$latitude), proj4string = CRS("+proj=longlat +ellps=WGS84"))
UTM.US <- spTransform(coord_invert.cari.out, CRS("+init=epsg:2163"))

dist_matrix.cari.out <- as.matrix(dist(data.frame(UTM.US)))
inv_dist.cari.out <- 1/dist_matrix.cari.out
diag(inv_dist.cari.out) <- 0
inv_dist.cari.out[is.infinite(inv_dist.cari.out)] <- 0

Moran.I(conductance.cari.out$logged.conductance, inv_dist.cari.out, alternative = "two.sided") 
##pvalue << 0.05, so there is still signficant spatial autocorrelation - in this case it's oddly negative 
```
summary: 
- the pvalue for moran's i with all the sites included was << 0.05, indicating that there was spatial autocorrelation in the data
  + some points were more similar depending on distance
  + we noticed in the variogram that there was one outlier point (high semivariance, far distance), which might have been skewing the data 
- reran moran's i & a semivariogram, and found that spatial autocorrelation was still significant 
  + however, it was negatively significant - increasing distance between sites made variance decrease
- for the next step, we will use cari, but note this site skews the direction of the autocorrelation

##creating a model corrected for spatial autocorrelation
```{r}
#jittering coords
##because the same location was sampled over successive seasons/years, we have to jitter coordinates slightly
set.seed(777)

coord_invert_jitter <- SpatialPoints(cbind(jitter(df_ph_logged$longitude), jitter(df_ph_logged$latitude)), proj4string = CRS("+proj=longlat +ellps=WGS84"))

UTM.US.jitter <- spTransform(coord_invert_jitter, CRS("+init=epsg:2163"))
UTM.US.jitter

#binding together jittered location & conductance data set 
df_ph_logged_jitter <- df_ph_logged %>% 
  dplyr::select(-coords.x1, -coords.x2) %>% 
  bind_cols(as.data.frame(UTM.US.jitter))

#all model types 
conductance.exp <- gls(logged.conductance~siteID, 
               data=df_ph_logged_jitter, method="ML", 
               corr=corSpatial(form=~coords.x1+coords.x2, type ="exponential"))
conductance.gau <- gls(logged.conductance~siteID, 
               data=df_ph_logged_jitter, method="ML", 
               corr=corSpatial(form=~coords.x1+coords.x2, type ="gaussian"))
conductance.lin <- gls(logged.conductance~siteID, 
               data=df_ph_logged_jitter, method="ML", 
               corr=corSpatial(form=~coords.x1+coords.x2, type ="linear"))
conductance.rat <- gls(logged.conductance~siteID, 
               data=df_ph_logged_jitter, method="ML", 
               corr=corSpatial(form=~coords.x1+coords.x2, type ="rational"))
conductance.sph <- gls(logged.conductance~siteID, 
               data=df_ph_logged_jitter, method="ML", 
               corr=corSpatial(form=~coords.x1+coords.x2, type ="spherical"))

#testing for the best model
AICc(conductance.null, conductance.exp, conductance.gau, conductance.lin, conductance.rat, conductance.sph)

##since the gaussian model has the lowest AICc score, it has the best fit to the data

#best model to predict conductance is: 
conductance.gau <- gls(logged.conductance~siteID, 
               data=df_ph_logged_jitter, method="ML", 
               corr=corSpatial(form=~coords.x1+coords.x2, type ="gaussian"))
```
summary: 
- set a seed so that the long/lat points jitter the same way each time 
  + jittering the long/lat points allows us to run spatial correction models where the same location is sampled over months/years 
  + the coordinates are changed very slightly, but not in a way that affects our analyses
- next, we ran multiple "corrected" models, where the spatial autocorrelation is accounted for in the regression
  + the best model was the gaussian model 

# 3)	Measuring correlation between %Riparian Cover and other variables 
o	Measure correlation between riparian cover, land use, seasonality and site 
o	Account for spatial autocorrelation 
o	Identify which variable to include as random effect (whatever the fixed effect is) 

```{r}
riparian_data <- read.csv("FINAL.RIP.GEO.csv")
```

##visualizing spread of mean canopy observations
```{r}
##over months
riparian_data %>%
  mutate(month = as.factor(month)) %>%
  ggplot(aes(x = month, y = mean_canopy)) +
  geom_boxplot()
#March and April have lower canopy cover than other months

##over seasons
riparian_data_season %>%
  ggplot(aes(x = season, y = mean_canopy)) +
  geom_boxplot()
#not a noticable difference between seasons 

##over sites
riparian_data %>%
  ggplot(aes(x = siteID, y = mean_canopy))+
  geom_boxplot()
# sites differ in canopy cover 

##over years
riparian_data %>%
  mutate(year = as.factor(year)) %>%
  ggplot(aes(x = year, y = mean_canopy))+
  geom_boxplot()
#each year is similar 
```

##visualizing correlation between canopy cover and land use 
```{r}
# buildings
riparian_data %>%
  ggplot(aes(x = fraction_buildings, y = mean_canopy)) +
  geom_point() +
  labs(x="Fraction of buildings", y="Mean canopy cover (%)")

# roads
riparian_data %>%
  ggplot(aes(x = fraction_roads, y = mean_canopy)) +
  geom_point() +
  labs(x="Fraction of roads", y="Mean canopy cover (%)")

# parks/lawns
riparian_data %>%
  ggplot(aes(x = fraction_parklawn, y = mean_canopy)) +
  geom_point() +
  labs(x="Fraction of parks/lawns", y="Mean canopy cover (%)")

# agriculture
riparian_data %>%
  ggplot(aes(x = fraction_agriculture, y = mean_canopy)) +
  geom_point() +
  labs(x="Fraction of agriculture", y="Mean canopy cover (%)")

# industry
riparian_data %>%
  ggplot(aes(x = fraction_industry, y = mean_canopy)) +
  geom_point() +
  labs(x="Fraction of industry", y="Mean canopy cover (%)")
```

##spatial autocorrelation between riparian cover and sites
```{r}
Moran.I(riparian_data$mean_canopy, inv.dist, alternative = "two.sided")
#p.value = 0.603 so its not spatially correlated. 
```

##calculating correlation between canopy cover & land use types 
```{r}
#buildings
cor(riparian_data$fraction_buildings, riparian_data$mean_canopy, use="complete.obs")

#agriculture 
cor(riparian_data$fraction_agriculture, riparian_data$mean_canopy, use="complete.obs")

#industry
cor(riparian_data$fraction_industry, riparian_data$mean_canopy, use="complete.obs")

#lawns/parks
cor(riparian_data$fraction_parklawn, riparian_data$mean_canopy, use="complete.obs")

#roads
cor(riparian_data$fraction_roads, riparian_data$mean_canopy, use="complete.obs")
```

# 4)	Measuring Simpson’s Index across sites 
o	Convert species data to a matrix 
o	Spatial autocorrelation 
o	Use Simpson’s function to quantify richness 

## Simpson's index of diversity
```{r}
ig_all <- read_csv("invert.genera.season.csv", col_names = TRUE)
ig_simpson <- ig_all %>%
  group_by(siteID, collectDate) %>%
  summarise(simpson = diversity(estimatedTotalCount, index = "simpson"),
            richness = specnumber(estimatedTotalCount)) 
ig_all_simpson <- inner_join(ig_all, ig_simpson)
```
## Is the simpson's index spatially autocorrelated?
### richness 
```{r}
Moran.I(ig_all_simpson$richness, inv.dist, alternative = "two.sided") # p < 0.05, therefore it is spatially autocorrelated.
```
### simpson's index
```{r}
Moran.I(ig_all_simpson$simpson, inv.dist, alternative = "two.sided") 
# p = 0.01546975 < 0.05, therefore it is spatially autocorrelated.
```

# 5)	Measuring Shannon’s Index across sites 
o	Convert species data to a matrix 
o	Spatial autocorrelation 
o	Use Simpson’s function to quantify richness 

## Evenness
```{r}
invert_evenness <- invert_genera %>%
  group_by(siteID, collectDate) %>%
  summarise(shannon = diversity(estimatedTotalCount, index = "shannon"),
            evenness = shannon/log(specnumber(estimatedTotalCount))) 
```

## Is evenness spatially autocorrelated?
```{r}
Moran.I(invert_evenness_coords$evenness, inv.dist, alternative = "two.sided") ## p = 0.432 so it is not.
```

# 4)	running SEM models 
```{r}
##modifying data 
merged.final.file<- read.csv('merged.final.file.csv')

merged.final.file <- merged.final.file %>%  
  distinct()

merged.final.file <- merged.final.file %>% 
  mutate(as.factor(siteID)) %>% 
  mutate(numericsiteID = as.numeric(siteID))
merged.final.file %>% 
  mutate(numericsiteID = as.factor(siteID)) %>% 
  mutate(numericsiteID = as.numeric(siteID))

merged.final.file <- merged.final.file %>% 
  mutate(as.factor(season)) %>% 
  mutate(numericseason = as.numeric(season))

##looking at correlation between variables 
merged.final.file %>% 
  dplyr::select(c(richness, pH, conductance, mean_canopy)) %>%
  chart.Correlation(.) 
merged.final.file %>% 
  dplyr::select(c(evenness, pH, conductance, mean_canopy)) %>%
  chart.Correlation(.) # Not very normal

##SEM for richness
sem_richness <- '
#regression 
richness ~ a*pH + b*mean_canopy + c*conductance
conductance ~ h*numericsiteID
pH ~ g*numericsiteID
pH ~ i*numericseason 

#correlations 
pH ~~ f*conductance
mean_canopy ~~ e*conductance 
mean_canopy ~~ d*pH 

#defined parameters 
total.pH := a + (d*b) + (f*c) + (d*e*c) + (f*e*b)
total.canopy := b + (c*e) + (d*a)
total.conductance := c + (e*b) + (f*a) + (e*d*a) + (f*e*b)
total.siteID := g*a + g*d*b + g*d*e*c + g*f*c + g*f*e*b + h*c + h*e*b + h*e*d*a + h*f*a + h*f*d*b
total.season := i*f*c + i*f*e*b + i*a + i*d*b + i*d*e*c
'

set.seed(77)
sem_richness_fit <- sem(sem_richness, data=merged.final.file, se="boot", bootstrap=100, verbose=TRUE)
sem_richness_est <- parameterEstimates(sem_richness_fit, boot.ci.type="bca.simple", standardized = TRUE)
sem_richness_est

write.csv(sem_richness_est, "sem.richness.csv")

##SEM for evenness 
sem_evenness <- '
#regression 
evenness ~ a*pH + b*mean_canopy + c*conductance
conductance ~ h*numericsiteID
pH ~ g*numericsiteID
pH ~ i*numericseason 

#correlations 
pH ~~ f*conductance
mean_canopy ~~ e*conductance 
mean_canopy ~~ d*pH 

#defined parameters 
total.pH := a + (d*b) + (f*c) + (d*e*c) + (f*e*b)
total.canopy := b + (c*e) + (d*a)
total.conductance := c + (e*b) + (f*a) + (e*d*a) + (f*e*b)
total.siteID := g*a + g*d*b + g*d*e*c + g*f*c + g*f*e*b + h*c + h*e*b + h*e*d*a + h*f*a + h*f*d*b
total.season := i*f*c + i*f*e*b + i*a + i*d*b + i*d*e*c
'

set.seed(77)
sem_evenness_fit <- sem(sem_evenness, data=merged.final.file, se="boot", bootstrap=100, verbose=TRUE)
sem_evenness_est <- parameterEstimates(sem_evenness_fit, boot.ci.type="bca.simple", standardized = TRUE)
sem_evenness_est

getwd()
write_csv(sem_evenness_est, "sem.evenness.csv")
```

```{r}
# plots 
merged.final.file %>% 
  ggplot(aes(x= pH, y= richness)) + 
  geom_point() + 
  geom_smooth(method= 'lm')


merged.final.file %>% 
  ggplot(aes(x= conductance, y= richness)) + 
  geom_point() + 
  geom_smooth(method= 'lm')

merged.final.file %>% 
  ggplot(aes(x= mean_canopy, y= richness)) + 
  geom_point() + 
  geom_smooth(method= 'lm')


merged.final.file %>% 
  ggplot(aes(x= pH, y= evenness)) + 
  geom_point() + 
  geom_smooth(method= 'lm')


merged.final.file %>% 
  ggplot(aes(x= conductance, y= evenness)) + 
  geom_point() + 
  geom_smooth(method= 'lm')

merged.final.file %>% 
  ggplot(aes(x= mean_canopy, y= evenness)) + 
  geom_point() + 
  geom_smooth(method= 'lm')
```
