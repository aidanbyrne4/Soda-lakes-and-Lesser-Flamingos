### R script for linear mixed-effect modelling to determine the drivers influencing Lesser Flamingo abundance at East African soda lakes
### Data availability is described in the STAR Methods section 
### Author: Aidan Byrne - aidan.byrne@kcl.ac.uk 
### Last modified: 03/11/2023

############################
###### load packages #######
############################

library(lme4) ## for linear mixed-effect modelling 
library(lmerTest) ## gives p values for lmer output when using summary()
library(olsrr) ## for testing normality of error distribution in log transformed models
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

## log transformed counts 

transformed <- completedata

transformed$Flamingo_count_IWC <- log(transformed$Flamingo_count_IWC)

## check distribution of log transformed data 

hist(transformed$Flamingo_count_IWC)

# remove scientific notation

options(scipen=10000)

##################################################
################# modelling ######################
##################################################

#### models for transformed data


#### Model 1 

model1 <-lmer(Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation
              + Rainfall_on_lake + Water_temp_Celsius + Windspeed + (1|Site),
              data = transformed, REML = TRUE, verbose = FALSE)

#### Model 2

model2 <- lmer(Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation
               + Rainfall_on_lake + Water_temp_Celsius + Windspeed + High_biomass_area + (1|Site),
               data = transformed, REML = TRUE, verbose = FALSE)

#### Model 3 

model3 <-lmer(Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation
              + Rainfall_on_lake + Water_temp_Celsius + Windspeed + Microphytobenthos_area + (1|Site),
              data = transformed, REML = TRUE, verbose = FALSE)

#### Model 4 

model4 <-lmer(Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation
              + Rainfall_on_lake + Water_temp_Celsius + Windspeed + Chla_water_Tebbs_method + (1|Site),
              data = transformed, REML = TRUE, verbose = FALSE)

#### Model 5 

model5 <-lmer(Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation
              + Rainfall_on_lake + Water_temp_Celsius + Windspeed + High_biomass_area + Sum_highbiomass_elsewhere + (1|Site),
              data = transformed, REML = TRUE, verbose = FALSE)

#### Model 6 

model6 <-lmer(Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation
              + Rainfall_on_lake + Water_temp_Celsius + Windspeed + Microphytobenthos_area + Sum_highbiomass_elsewhere + (1|Site),
              data = transformed, REML = TRUE, verbose = FALSE)

#### Model 7 

model7 <-lmer(Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation
              + Rainfall_on_lake + Water_temp_Celsius + Windspeed + Chla_water_Tebbs_method + Sum_highbiomass_elsewhere + (1|Site),
              data = transformed, REML = TRUE, verbose = FALSE)

#### Model 8 

model8 <-lmer(Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation
              + Rainfall_on_lake + Water_temp_Celsius + Windspeed + High_biomass_area + Weighted_sum_highbiomass_elsewhere + (1|Site),
              data = transformed, REML = TRUE, verbose = FALSE)


#### Model 9 

model9 <-lmer(Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation
              + Rainfall_on_lake + Water_temp_Celsius + Windspeed + Microphytobenthos_area + Weighted_sum_highbiomass_elsewhere + (1|Site),
              data = transformed, REML = TRUE, verbose = FALSE)

#### Model 10

model10 <-lmer(Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation
               + Rainfall_on_lake + Water_temp_Celsius + Windspeed + Chla_water_Tebbs_method + Weighted_sum_highbiomass_elsewhere + (1|Site),
               data = transformed, REML = TRUE, verbose = FALSE)

#### Model 11

model11 <-lmer(Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation
               + Rainfall_on_lake + Water_temp_Celsius + Windspeed + High_biomass_area + Mean_highbiomass_elsewhere + (1|Site),
               data = transformed, REML = TRUE, verbose = FALSE)

#### Model 12

model12 <-lmer(Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation
               + Rainfall_on_lake + Water_temp_Celsius + Windspeed + Microphytobenthos_area + Mean_highbiomass_elsewhere + (1|Site),
               data = transformed, REML = TRUE, verbose = FALSE)

#### Model 13

model13 <-lmer(Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation
               + Rainfall_on_lake + Water_temp_Celsius + Windspeed + Chla_water_Tebbs_method + Mean_highbiomass_elsewhere + (1|Site),
               data = transformed, REML = TRUE, verbose = FALSE)

#### Model 14

model14 <-lmer(Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation
               + Rainfall_on_lake + Water_temp_Celsius + Windspeed + High_biomass_area + Weighted_mean_highbiomass_elsewhere + (1|Site),
               data = transformed, REML = TRUE, verbose = FALSE)

#### Model 15

model15 <-lmer(Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation
               + Rainfall_on_lake + Water_temp_Celsius + Windspeed + Microphytobenthos_area + Weighted_mean_highbiomass_elsewhere + (1|Site),
               data = transformed, REML = TRUE, verbose = FALSE)

#### Model 16

model16 <-lmer(Flamingo_count_IWC ~ Total_lake_area + Season + Percent_shoreline_vegetation
               + Rainfall_on_lake + Water_temp_Celsius + Windspeed + Chla_water_Tebbs_method + Weighted_mean_highbiomass_elsewhere + (1|Site),
               data = transformed, REML = TRUE, verbose = FALSE)

############################
###### Model outputs #######
############################

## print outputs

summary(model16) ##print output - model16 was the best model 

AIC(model16) ## print AIC

##AIC comparison

modellist <- list(model1, model2, model3, model4, model5, model6,
                  model7, model8, model9, model10, model11, model12,
                  model13, model14, model15, model16)

aictab(modellist, modnames = NULL, second.ord = TRUE, nobs = NULL,
       sort = TRUE)

## model averaging

bestmodels <- list(model16, model14, model4, model8, model10, model2, model1)

summary(model.avg(bestmodels)) ##muin package

summary(model.avg(modellist, subset = delta < 2)) ## only models with AIC within 2 of top model 

##################################
###### Model checking ############
##################################

### test for normality assumption of error distribution 

qqnorm(resid(model16)) ##qqplot for residuals
qqline(resid(model16))

hist(resid(model16)) ##histogram for residuals

plot(model16)###fitted vs residuals plot

plot(cooks.distance(model16)) #cooks distance

lev<-hat(model.matrix(model16))#calculate leverage
plot(resid(model16,type="pearson")~lev,las=1,ylab="Standardised residuals",xlab="Leverage") #Plot leverage against standardised residuals

plot(model16, ##scale location plot
     sqrt(abs(resid(.)))~fitted(.),
     type=c("p","smooth"), col.line=1)
