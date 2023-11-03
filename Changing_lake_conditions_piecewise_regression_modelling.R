### R script for segmented regression modelling of trends in lake surface areas and chlorophyll-a concentrations for East African soda lakes
### Data availability: Data availability is described in the STAR methods section 
### Author: Aidan Byrne - aidan.byrne@kcl.ac.uk 
### Last modified: 03/11/2023

##############################
####### load packages ########
##############################

library(dplyr) ## for data manipulation 
library(lme4) ## for linear models 
library(lmerTest) ## gives p values for lmer output when using summary()
library(olsrr) ## for testing normality of error distribution in log transformed models
library(AICcmodavg) ## for AIC model comparison 
library(MuMIn) ## for model averaging 
library(cAIC4) ## for model averaging 
library(segmented) ## for piecewise regression models 

#######################################################
############# read data and pre-process ###############
#######################################################

data <- read.csv("Soda_lakes_Chla_and_surface_areas_1999_2022.csv", header=TRUE)

head(data)

## convert date column to date class and add Year column

data$Date <- as.Date(paste(data$Date,"-01",sep=""))

data$Date <- as.Date(data$Date,
                         format = "%Y/%m/%d")

data$Year <- format(data$Date,
                        format = "%Y")

## convrt year to numeric for models

data$Year <- as.numeric(data$Year)

# remove scientific notation

options(scipen=10000)

####################################################
#### Subset data for individual lake modelling #####
####################################################

## Define site name for individual lake models

Bogorianame = "Bogoria"

Bogoria <- data[data$Site == Bogorianame, ]

head(Bogoria)

##

Elmenteitaname = "Elmenteita"

Elmenteita <- data[data$Site == Elmenteitaname, ]

head(Elmenteita)

##

Magadiname = "Magadi"

Magadi <- data[data$Site == Magadiname, ]

head(Magadi)

##

Magadi_Ngorongoroname = "Magadi - Ngorongoro"

Magadi_Ngorongoro <- data[data$Site == Magadi_Ngorongoroname, ]

head(Magadi_Ngorongoro)

##

Nakuruname = "Nakuru"

Nakuru <- data[data$Site == Nakuruname, ]

head(Nakuru)

##

Oloidienname = "Oloidien"

Oloidien <- data[data$Site == Oloidienname, ]

head(Oloidien)

##

Sonachiname = "Sonachi"

Sonachi <- data[data$Site == Sonachiname, ]

head(Sonachi)

## 

Solainame = "Solai"

Solai <- data[data$Site == Solainame, ]

head(Solai)

##

Simbi_Cratername = "Simbi Crater"

Simbi_Crater <- data[data$Site == Simbi_Cratername, ]

head(Simbi_Crater)

##

Natronname = "Natron"

Natron <- data[data$Site == Natronname, ]

head(Natron)

##

Manyaraname = "Manyara"

Manyara <- data[data$Site == Manyaraname, ]

head(Manyara)

##

Eyasiname = "Eyasi"

Eyasi <- data[data$Site == Eyasiname, ]

head(Eyasi)

##

Burunginame = "Burungi"

Burungi <- data[data$Site == Burunginame, ]

head(Burungi)

##

Green_lakename = "Green lake"

Green_lake <- data[data$Site == Green_lakename, ]

head(Green_lake)

##

Abijataname = "Abijata"

Abijata <- data[data$Site == Abijataname, ]

head(Abijata)

##

Shallaname = "Shalla"

Shalla <- data[data$Site == Shallaname, ]

head(Shalla)

##

Chituname = "Chitu"

Chitu <- data[data$Site == Chituname, ]

head(Chitu)

##

Chew_Bahirname = "Chew Bahir"

Chew_Bahir <- data[data$Site == Chew_Bahirname, ]

head(Chew_Bahir)

##

Logipiname = "Logipi"

Logipi <- data[data$Site == Logipiname, ]

head(Logipi)

##

Empakainame = "Empakai"

Empakai <- data[data$Site == Empakainame, ]

head(Empakai)

##

Big_Momellaname = "Big Momella"

Big_Momella <- data[data$Site == Big_Momellaname, ]

head(Big_Momella)

##

Small_Momellaname = "Small Momella"

Small_Momella <- data[data$Site == Small_Momellaname, ]

head(Small_Momella)


##################################################
################# modelling ######################
##################################################

## develop 8 models for each of the 22 soda lakes - 4 chl-a and 4 surface area. 
## Each set of 4 contains 1 linear model, 1 breakpoint segmented regression, 2 breakpoint segmented regression, 3 breakpoint segmented regression.
## Note: Model 2 for each lake was originally a polynomial regression, however was removed during analyses

###############  bog

#### Model 1 

bog1 <-lm(Chla_water_Tebbs_method ~ Year,
            data = Bogoria)

#### bog 3 

bog3 <- segmented(bog1, seg.Z = ~Year, npsi= 1)

#### bog 4 

bog4 <- segmented(bog1, seg.Z = ~Year, npsi= 2)

#### bog 5

bog5 <- segmented(bog1, seg.Z = ~Year, npsi= 3)

#### bog 6 

bog6 <-lm(Surface.area ~ Year,
            data = Bogoria)

#### bog 8 

bog8 <- segmented(bog6, seg.Z = ~Year, npsi= 1)

#### bog 9 

bog9 <- segmented(bog6, seg.Z = ~Year, npsi= 2)

#### bog 10

bog10 <- segmented(bog6, seg.Z = ~Year, npsi= 3)

############### elm

#### elm 1 

elm1 <-lm(Chla_water_Tebbs_method ~ Year,
            data = Elmenteita)

#### elm 3 

elm3 <- segmented(elm1, seg.Z = ~Year, npsi= 1)

#### elm 4 

elm4 <- segmented(elm1, seg.Z = ~Year, npsi= 2)

#### elm 5

elm5 <- segmented(elm1, seg.Z = ~Year, npsi= 3)

#### elm 6 

elm6 <-lm(Surface.area ~ Year,
            data = Elmenteita)

#### elm 8 

elm8 <- segmented(elm6, seg.Z = ~Year, npsi= 1)

#### elm 9 

elm9 <- segmented(elm6, seg.Z = ~Year, npsi= 2)

#### elm 10

elm10 <- segmented(elm6, seg.Z = ~Year, npsi= 3)

################ mag

#### mag 1 

mag1 <-lm(Chla_water_Tebbs_method ~ Year,
            data = Magadi)

#### mag 3 

mag3 <- segmented(mag1, seg.Z = ~Year, npsi= 1)

#### mag 4 

mag4 <- segmented(mag1, seg.Z = ~Year, npsi= 2)

#### mag 5

mag5 <- segmented(mag1, seg.Z = ~Year, npsi= 3)

#### mag 6 

mag6 <-lm(Surface.area ~ Year,
            data = Magadi)

#### mag 8 

mag8 <- segmented(mag6, seg.Z = ~Year, npsi= 1)

#### mag 9 

mag9 <- segmented(mag6, seg.Z = ~Year, npsi= 2)

#### mag 10

mag10 <- segmented(mag6, seg.Z = ~Year, npsi= 3)

############### mag ngo

#### magngo 1 

magngo1 <-lm(Chla_water_Tebbs_method ~ Year,
            data = Magadi_Ngorongoro)

#### magngo 3 

magngo3 <- segmented(magngo1, seg.Z = ~Year, npsi= 1)

#### magngo 4 

magngo4 <- segmented(magngo1, seg.Z = ~Year, npsi= 2)

#### magngo 5

magngo5 <- segmented(magngo1, seg.Z = ~Year, npsi= 3)

#### magngo 6 

magngo6 <-lm(Surface.area ~ Year,
            data = Magadi_Ngorongoro)

#### magngo 8 

magngo8 <- segmented(magngo6, seg.Z = ~Year, npsi= 1)

#### magngo 9 

magngo9 <- segmented(magngo6, seg.Z = ~Year, npsi= 2)

#### magngo 10

magngo10 <- segmented(magngo6, seg.Z = ~Year, npsi= 3)

############## nak

#### nak 1 

nak1 <-lm(Chla_water_Tebbs_method ~ Year,
            data = Nakuru)

#### nak 3 

nak3 <- segmented(nak1, seg.Z = ~Year, npsi= 1)

#### nak 4 

nak4 <- segmented(nak1, seg.Z = ~Year, npsi= 2)

#### nak 5

nak5 <- segmented(nak1, seg.Z = ~Year, npsi= 3)

#### nak 6 

nak6 <-lm(Surface.area ~ Year,
            data = Nakuru)

#### nak 8 

nak8 <- segmented(nak6, seg.Z = ~Year, npsi= 1)

#### nak 9 

nak9 <- segmented(nak6, seg.Z = ~Year, npsi= 2)

#### nak 10

nak10 <- segmented(nak6, seg.Z = ~Year, npsi= 3)

############## oloi

#### olo 1 

olo1 <-lm(Chla_water_Tebbs_method ~ Year,
            data = Oloidien)

#### olo 3 

olo3 <- segmented(olo1, seg.Z = ~Year, npsi= 1)

#### olo 4 

olo4 <- segmented(olo1, seg.Z = ~Year, npsi= 2)

#### olo 5

olo5 <- segmented(olo1, seg.Z = ~Year, npsi= 3)

#### olo 6 

olo6 <-lm(Surface.area ~ Year,
            data = Oloidien)

#### olo 8 

olo8 <- segmented(olo6, seg.Z = ~Year, npsi= 1)

#### olo 9 

olo9 <- segmented(olo6, seg.Z = ~Year, npsi= 2)

#### olo 10

olo10 <- segmented(olo6, seg.Z = ~Year, npsi= 3)

############### son

#### son 1 

son1 <-lm(Chla_water_Tebbs_method ~ Year,
            data = Sonachi)

#### son 3 

son3 <- segmented(son1, seg.Z = ~Year, npsi= 1)

#### son 4 

son4 <- segmented(son1, seg.Z = ~Year, npsi= 2)

#### son 5

son5 <- segmented(son1, seg.Z = ~Year, npsi= 3)

#### son 6 

son6 <-lm(Surface.area ~ Year,
            data = Sonachi)

#### son 8 

son8 <- segmented(son6, seg.Z = ~Year, npsi= 1)

#### son 9 

son9 <- segmented(son6, seg.Z = ~Year, npsi= 2)

#### son 10

son10 <- segmented(son6, seg.Z = ~Year, npsi= 3)

############# new plots #######
#######

############### solai

#### sol 1 

sol1 <-lm(Chla_water_Tebbs_method ~ Year,
          data = Solai)

#### sol 3 

sol3 <- segmented(sol1, seg.Z = ~Year, npsi= 1)

#### sol 4 

sol4 <- segmented(sol1, seg.Z = ~Year, npsi= 2)

#### sol 5

sol5 <- segmented(sol1, seg.Z = ~Year, npsi= 3)

#### sol 6 

sol6 <-lm(Surface.area ~ Year,
          data = Solai)

#### sol 8 

sol8 <- segmented(sol6, seg.Z = ~Year, npsi= 1)

#### sol 9 

sol9 <- segmented(sol6, seg.Z = ~Year, npsi= 2)

#### sol 10

sol10 <- segmented(sol6, seg.Z = ~Year, npsi= 3)

############### sim

#### sim 1 

sim1 <-lm(Chla_water_Tebbs_method ~ Year,
          data = Simbi_Crater)

#### sim 3 

sim3 <- segmented(sim1, seg.Z = ~Year, npsi= 1)

#### sim 4 

sim4 <- segmented(sim1, seg.Z = ~Year, npsi= 2)

#### sim 5

sim5 <- segmented(sim1, seg.Z = ~Year, npsi= 3)

#### sim 6 

sim6 <-lm(Surface.area ~ Year,
          data = Simbi_Crater)

#### sim 8 

sim8 <- segmented(sim6, seg.Z = ~Year, npsi= 1)

#### sim 9 

sim9 <- segmented(sim6, seg.Z = ~Year, npsi= 2)

#### sim 10

sim10 <- segmented(sim6, seg.Z = ~Year, npsi= 3)

############### nat

#### nat 1 

nat1 <-lm(Chla_water_Tebbs_method ~ Year,
          data = Natron)

#### nat 3 

nat3 <- segmented(nat1, seg.Z = ~Year, npsi= 1)

#### nat 4 

nat4 <- segmented(nat1, seg.Z = ~Year, npsi= 2)

#### nat 5

nat5 <- segmented(nat1, seg.Z = ~Year, npsi= 3)

#### nat 6 

nat6 <-lm(Surface.area ~ Year,
          data = Natron)

#### nat 8 

nat8 <- segmented(nat6, seg.Z = ~Year, npsi= 1)

#### nat 9 

nat9 <- segmented(nat6, seg.Z = ~Year, npsi= 2)

#### nat 10

nat10 <- segmented(nat6, seg.Z = ~Year, npsi= 3)

############### man

#### man 1 

man1 <-lm(Chla_water_Tebbs_method ~ Year,
          data = Manyara)

#### man 3 

man3 <- segmented(man1, seg.Z = ~Year, npsi= 1)

#### man 4 

man4 <- segmented(man1, seg.Z = ~Year, npsi= 2)

#### man 5

man5 <- segmented(man1, seg.Z = ~Year, npsi= 3)

#### man 6 

man6 <-lm(Surface.area ~ Year,
          data = Manyara)

#### man 8 

man8 <- segmented(man6, seg.Z = ~Year, npsi= 1)

#### man 9 

man9 <- segmented(man6, seg.Z = ~Year, npsi= 2)

#### man 10

man10 <- segmented(man6, seg.Z = ~Year, npsi= 3)

############### eya

#### eya 1 

eya1 <-lm(Chla_water_Tebbs_method ~ Year,
          data = Eyasi)

#### eya 3 

eya3 <- segmented(eya1, seg.Z = ~Year, npsi= 1)

#### eya 4 

eya4 <- segmented(eya1, seg.Z = ~Year, npsi= 2)

#### eya 5

eya5 <- segmented(eya1, seg.Z = ~Year, npsi= 3)

#### eya 6 

eya6 <-lm(Surface.area ~ Year,
          data = Eyasi)

#### eya 8 

eya8 <- segmented(eya6, seg.Z = ~Year, npsi= 1)

#### eya 9 

eya9 <- segmented(eya6, seg.Z = ~Year, npsi= 2)

#### eya 10

eya10 <- segmented(eya6, seg.Z = ~Year, npsi= 3)

############### bur

#### bur 1 

bur1 <-lm(Chla_water_Tebbs_method ~ Year,
          data = Burungi)

#### bur 3 

bur3 <- segmented(bur1, seg.Z = ~Year, npsi= 1)

#### bur 4 

bur4 <- segmented(bur1, seg.Z = ~Year, npsi= 2)

#### bur 5

bur5 <- segmented(bur1, seg.Z = ~Year, npsi= 3)

#### bur 6 

bur6 <-lm(Surface.area ~ Year,
          data = Burungi)

#### bur 8 

bur8 <- segmented(bur6, seg.Z = ~Year, npsi= 1)

#### bur 9 

bur9 <- segmented(bur6, seg.Z = ~Year, npsi= 2)

#### bur 10

bur10 <- segmented(bur6, seg.Z = ~Year, npsi= 3)

####
#### newer lakes 
####

############### gre

#### gre 1 

gre1 <-lm(Chla_water_Tebbs_method ~ Year,
          data = Green_lake)

#### gre 3 

gre3 <- segmented(gre1, seg.Z = ~Year, npsi= 1)

#### gre 4 

gre4 <- segmented(gre1, seg.Z = ~Year, npsi= 2)

#### gre 5

gre5 <- segmented(gre1, seg.Z = ~Year, npsi= 3)

#### gre 6 

gre6 <-lm(Surface.area ~ Year,
          data = Green_lake)

#### gre 8 

gre8 <- segmented(gre6, seg.Z = ~Year, npsi= 1)

#### gre 9 

gre9 <- segmented(gre6, seg.Z = ~Year, npsi= 2)

#### gre 10

gre10 <- segmented(gre6, seg.Z = ~Year, npsi= 3)

############### abi

#### abi 1 

abi1 <-lm(Chla_water_Tebbs_method ~ Year,
          data = Abijata)

#### abi 3 

abi3 <- segmented(abi1, seg.Z = ~Year, npsi= 1)

#### abi 4 

abi4 <- segmented(abi1, seg.Z = ~Year, npsi= 2)

#### abi 5

abi5 <- segmented(abi1, seg.Z = ~Year, npsi= 3)

#### abi 6 

abi6 <-lm(Surface.area ~ Year,
          data = Abijata)

#### abi 8 

abi8 <- segmented(abi6, seg.Z = ~Year, npsi= 1)

#### abi 9 

abi9 <- segmented(abi6, seg.Z = ~Year, npsi= 2)

#### abi 10

abi10 <- segmented(abi6, seg.Z = ~Year, npsi= 3)

############### sha

#### sha 1 

sha1 <-lm(Chla_water_Tebbs_method ~ Year,
          data = Shalla)

#### sha 3 

sha3 <- segmented(sha1, seg.Z = ~Year, npsi= 1)

#### sha 4 

sha4 <- segmented(sha1, seg.Z = ~Year, npsi= 2)

#### sha 5

sha5 <- segmented(sha1, seg.Z = ~Year, npsi= 3)

#### sha 6 

sha6 <-lm(Surface.area ~ Year,
          data = Shalla)

#### sha 8 

sha8 <- segmented(sha6, seg.Z = ~Year, npsi= 1)

#### sha 9 

sha9 <- segmented(sha6, seg.Z = ~Year, npsi= 2)

#### sha 10

sha10 <- segmented(sha6, seg.Z = ~Year, npsi= 3)

############### chi

#### chi 1 

chi1 <-lm(Chla_water_Tebbs_method ~ Year,
          data = Chitu)

#### chi 3 

chi3 <- segmented(chi1, seg.Z = ~Year, npsi= 1)

#### chi 4 

chi4 <- segmented(chi1, seg.Z = ~Year, npsi= 2)

#### chi 5

chi5 <- segmented(chi1, seg.Z = ~Year, npsi= 3)

#### chi 6 

chi6 <-lm(Surface.area ~ Year,
          data = Chitu)

#### chi 8 

chi8 <- segmented(chi6, seg.Z = ~Year, npsi= 1)

#### chi 9 

chi9 <- segmented(chi6, seg.Z = ~Year, npsi= 2)

#### chi 10

chi10 <- segmented(chi6, seg.Z = ~Year, npsi= 3)

############### che

#### che 1 

che1 <-lm(Chla_water_Tebbs_method ~ Year,
          data = Chew_Bahir)

#### che 3 

che3 <- segmented(che1, seg.Z = ~Year, npsi= 1)

#### che 4 

che4 <- segmented(che1, seg.Z = ~Year, npsi= 2)

#### che 5

che5 <- segmented(che1, seg.Z = ~Year, npsi= 3)

#### che 6 

che6 <-lm(Surface.area ~ Year,
          data = Chew_Bahir)

#### che 8 

che8 <- segmented(che6, seg.Z = ~Year, npsi= 1)

#### che 9 

che9 <- segmented(che6, seg.Z = ~Year, npsi= 2)

#### che 10

che10 <- segmented(che6, seg.Z = ~Year, npsi= 3)

############### log

#### log 1 

log1 <-lm(Chla_water_Tebbs_method ~ Year,
          data = Logipi)

#### log 3 

log3 <- segmented(log1, seg.Z = ~Year, npsi= 1)

#### log 4 

log4 <- segmented(log1, seg.Z = ~Year, npsi= 2)

#### log 5

log5 <- segmented(log1, seg.Z = ~Year, npsi= 3)

#### log 6 

log6 <-lm(Surface.area ~ Year,
          data = Logipi)

#### log 8 

log8 <- segmented(log6, seg.Z = ~Year, npsi= 1)

#### log 9 

log9 <- segmented(log6, seg.Z = ~Year, npsi= 2)

#### log 10

log10 <- segmented(log6, seg.Z = ~Year, npsi= 3)

############### emp

#### emp 1 

emp1 <-lm(Chla_water_Tebbs_method ~ Year,
          data = Empakai)

#### emp 3 

emp3 <- segmented(emp1, seg.Z = ~Year, npsi= 1)

#### emp 4 

emp4 <- segmented(emp1, seg.Z = ~Year, npsi= 2)

#### emp 5

emp5 <- segmented(emp1, seg.Z = ~Year, npsi= 3)

#### emp 6 

emp6 <-lm(Surface.area ~ Year,
          data = Empakai)

#### emp 8 

emp8 <- segmented(emp6, seg.Z = ~Year, npsi= 1)

#### emp 9 

emp9 <- segmented(emp6, seg.Z = ~Year, npsi= 2)

#### emp 10

emp10 <- segmented(emp6, seg.Z = ~Year, npsi= 3)

############### big

#### big 1 

big1 <-lm(Chla_water_Tebbs_method ~ Year,
          data = Big_Momella)

#### big 3 

big3 <- segmented(big1, seg.Z = ~Year, npsi= 1)

#### big 4 

big4 <- segmented(big1, seg.Z = ~Year, npsi= 2)

#### big 5

big5 <- segmented(big1, seg.Z = ~Year, npsi= 3)

#### big 6 

big6 <-lm(Surface.area ~ Year,
          data = Big_Momella)

#### big 8 

big8 <- segmented(big6, seg.Z = ~Year, npsi= 1)

#### big 9 

big9 <- segmented(big6, seg.Z = ~Year, npsi= 2)

#### big 10

big10 <- segmented(big6, seg.Z = ~Year, npsi= 3)

############### sma

#### sma 1 

sma1 <-lm(Chla_water_Tebbs_method ~ Year,
          data = Small_Momella)

#### sma 3 

sma3 <- segmented(sma1, seg.Z = ~Year, npsi= 1)

#### sma 4 

sma4 <- segmented(sma1, seg.Z = ~Year, npsi= 2)

#### sma 5

sma5 <- segmented(sma1, seg.Z = ~Year, npsi= 3)

#### sma 6 

sma6 <-lm(Surface.area ~ Year,
          data = Small_Momella)

#### sma 8 

sma8 <- segmented(sma6, seg.Z = ~Year, npsi= 1)

#### sma 9 

sma9 <- segmented(sma6, seg.Z = ~Year, npsi= 2)

#### sma 10

sma10 <- segmented(sma6, seg.Z = ~Year, npsi= 3)

#############################
###### Print outputs ########
#############################

## print outputs 

summary(sma10)
AIC(sma10)

##AIC comparison

modellistfood <- list(sma1, sma3, sma4, sma5)

modellistwater <- list(sma6, sma8, sma9, sma10)

aictab(modellistfood, modnames = NULL, second.ord = TRUE, nobs = NULL,
       sort = TRUE)

###################################
######### Model checking ##########
###################################

qqnorm(resid(sma5)) ##qqplot for residuals
qqline(resid(sma5))

hist(resid(sma5)) ##histogram for residuals

plot(model5)###fitted vs residuals plot

plot(cooks.distance(sma5)) #cooks distance

lev<-hat(model.matrix(sma5))#calculate leverage
plot(resid(sma5,type="pearson")~lev,las=1,ylab="Standardised residuals",xlab="Leverage") #Plot leverage against standardised residuals

plot(sma5, ##scale location plot
     sqrt(abs(resid(.)))~fitted(.),
     type=c("p","smooth"), col.line=1)
