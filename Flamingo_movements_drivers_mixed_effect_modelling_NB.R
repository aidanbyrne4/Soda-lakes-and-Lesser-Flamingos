### R script for generalised linear mixed modelling to determine the drivers influencing Lesser Flamingo abundance at East African soda lakes
### Data availability is described in the STAR Methods section 
### Author: Aidan Byrne - aidan.byrne@kcl.ac.uk 
### Last modified: 22/01/2024

############################
###### load packages #######
############################

library(GLMMadaptive) ## for generalised linear mixed models with negative binomial family 
library(AICcmodavg) ## for AIC model comparison 
library(MuMIn) ## for model averaging 
library(cAIC4) ## for model averaging 

#########################################################################
############# read data, remove NA values and pre-process ###############
#########################################################################

data <- read.csv('Flamingo_count_and_drivers_data.csv', header = TRUE) ## cannot share data due to count data ownership rights - see Data availability section of STAR methods

completedata <- na.omit(data)

head(completedata)

str(completedata)

# attach dataset

attach(completedata)

## rescale predictors for model fitting - predictors are on different scales 

completedata_scaled <- completedata

completedata_scaled[, c("Total_lake_area", "Rainfall_on_lake", "Water_temp_Celsius", "Windspeed", "Chla_water_Tebbs_method", "Weighted_mean_highbiomass_elsewhere", "Season", "High_biomass_area", "Percent_shoreline_vegetation", "Weighted_sum_highbiomass_elsewhere", "Mean_highbiomass_elsewhere", "Sum_highbiomass_elsewhere", "Microphytobenthos_area")] <- scale(completedata[, c("Total_lake_area", "Rainfall_on_lake", "Water_temp_Celsius", "Windspeed", "Chla_water_Tebbs_method", "Weighted_mean_highbiomass_elsewhere", "Season", "High_biomass_area", "Percent_shoreline_vegetation", "Weighted_sum_highbiomass_elsewhere", "Mean_highbiomass_elsewhere", "Sum_highbiomass_elsewhere", "Microphytobenthos_area")])

# remove scientific notation

options(scipen=10000)

##################################################
################# modelling ######################
##################################################

#### Negative binomial Generalised linear mixed models

# Model 1 
formula <- Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation + Rainfall_on_lake + Water_temp_Celsius + Windspeed

model1_nb <- mixed_model(
  fixed = formula,
  data = completedata_scaled,
  random = ~1|Site,
  family = negative.binomial()
)

# Model 2
formula2 <- Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation + Rainfall_on_lake + Water_temp_Celsius + Windspeed + High_biomass_area

model2_nb <- mixed_model(
  fixed = formula2,
  data = completedata_scaled,
  random = ~1|Site,
  family = negative.binomial()
)

# Model 3
formula3 <- Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation + Rainfall_on_lake + Water_temp_Celsius + Windspeed + Microphytobenthos_area

model3_nb <- mixed_model(
  fixed = formula3,
  data = completedata_scaled,
  random = ~1|Site,
  family = negative.binomial()
)

# Model 4
formula4 <- Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation + Rainfall_on_lake + Water_temp_Celsius + Windspeed + Chla_water_Tebbs_method

model4_nb <- mixed_model(
  fixed = formula4,
  data = completedata_scaled, ####### removed scale for plot
  random = ~1|Site,
  family = negative.binomial()
)

# Model 5
formula5 <- Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation +
  Rainfall_on_lake + Water_temp_Celsius + Windspeed + High_biomass_area +
  Sum_highbiomass_elsewhere

model5_nb <- mixed_model(
  fixed = formula5,
  data = completedata_scaled,
  random = ~1|Site,
  family = negative.binomial()
)

# Model 6
formula6 <- Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation +
  Rainfall_on_lake + Water_temp_Celsius + Windspeed + Microphytobenthos_area +
  Sum_highbiomass_elsewhere

model6_nb <- mixed_model(
  fixed = formula6,
  data = completedata_scaled,
  random = ~1|Site,
  family = negative.binomial()
)

# Model 7
formula7 <- Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation +
  Rainfall_on_lake + Water_temp_Celsius + Windspeed + Chla_water_Tebbs_method +
  Sum_highbiomass_elsewhere

model7_nb <- mixed_model(
  fixed = formula7,
  data = completedata_scaled,
  random = ~1|Site,
  family = negative.binomial()
)

# Model 8
formula8 <- Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation +
  Rainfall_on_lake + Water_temp_Celsius + Windspeed + High_biomass_area +
  Weighted_sum_highbiomass_elsewhere

model8_nb <- mixed_model(
  fixed = formula8,
  data = completedata_scaled,
  random = ~1|Site,
  family = negative.binomial()
)

# Model 9
formula9 <- Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation +
  Rainfall_on_lake + Water_temp_Celsius + Windspeed + Microphytobenthos_area +
  Weighted_sum_highbiomass_elsewhere

model9_nb <- mixed_model(
  fixed = formula9,
  data = completedata_scaled,
  random = ~1|Site,
  family = negative.binomial()
)

# Model 10
formula10 <- Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation +
  Rainfall_on_lake + Water_temp_Celsius + Windspeed + Chla_water_Tebbs_method +
  Weighted_sum_highbiomass_elsewhere

model10_nb <- mixed_model(
  fixed = formula10,
  data = completedata_scaled,
  random = ~1|Site,
  family = negative.binomial()
)

# Model 11
formula11 <- Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation +
  Rainfall_on_lake + Water_temp_Celsius + Windspeed + High_biomass_area +
  Mean_highbiomass_elsewhere

model11_nb <- mixed_model(
  fixed = formula11,
  data = completedata_scaled,
  random = ~1|Site,
  family = negative.binomial()
)

# Model 12
formula12 <- Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation +
  Rainfall_on_lake + Water_temp_Celsius + Windspeed + Microphytobenthos_area +
  Mean_highbiomass_elsewhere

model12_nb <- mixed_model(
  fixed = formula12,
  data = completedata_scaled,
  random = ~1|Site,
  family = negative.binomial()
)

# Model 13
formula13 <- Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation +
  Rainfall_on_lake + Water_temp_Celsius + Windspeed + Chla_water_Tebbs_method +
  Mean_highbiomass_elsewhere

model13_nb <- mixed_model(
  fixed = formula13,
  data = completedata_scaled,
  random = ~1|Site,
  family = negative.binomial()
)

# Model 14
formula14 <- Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation +
  Rainfall_on_lake + Water_temp_Celsius + Windspeed + High_biomass_area +
  Weighted_mean_highbiomass_elsewhere

model14_nb <- mixed_model(
  fixed = formula14,
  data = completedata_scaled,
  random = ~1|Site,
  family = negative.binomial()
)

# Model 15
formula15 <- Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation +
  Rainfall_on_lake + Water_temp_Celsius + Windspeed + Microphytobenthos_area +
  Weighted_mean_highbiomass_elsewhere

model15_nb <- mixed_model(
  fixed = formula15,
  data = completedata_scaled,
  random = ~1|Site,
  family = negative.binomial()
)

# Model 16
formula16 <- Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation +
  Rainfall_on_lake + Water_temp_Celsius + Windspeed + Chla_water_Tebbs_method +
  Weighted_mean_highbiomass_elsewhere

model16_nb <- mixed_model(
  fixed = formula16,
  data = completedata_scaled,
  random = ~1|Site,
  family = negative.binomial()
)

############################
###### Model outputs #######
############################

## print outputs

summary(model4_nb) ##print output - model4 was the best model 

AIC(model4_nb) ## print AIC

##AIC comparison

modellist <- list(model1_nb, model2_nb, model3_nb, model4_nb, model5_nb, model6_nb,
                  model7_nb, model8_nb, model9_nb, model10_nb, model11_nb, model12_nb,
                  model13_nb, model14_nb, model15_nb, model16_nb)

aictab(modellist, modnames = NULL, second.ord = TRUE, nobs = NULL,
       sort = TRUE)

## model averaging

bestmodels <- list(model4_nb, model7_nb, model10_nb, model13_nb, model16_nb)

summary(model.avg(bestmodels)) ##muin package

summary(model.avg(modellist, subset = delta < 2)) ## only models with AIC within 2 of top model 
