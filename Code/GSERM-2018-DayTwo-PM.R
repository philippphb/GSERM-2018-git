#########################################################
# GSERM - St. Gallen (2018)
#
# Binary & event count panel data models...
#
########################################################
# Load packages (install as needed), set options:

library(RCurl)
library(psych)
library(lme4)
library(plm)
library(gtools)
library(boot)
library(plyr)
library(texreg)
library(statmod)
library(prais)
library(nlme)
library(pcse)
library(plm)
library(glmmML)
library(censReg)

options(scipen = 6) # bias against scientific notation
options(digits = 2) # show fewer decimal places

#######################
# Segal data...

URL<-"https://raw.githubusercontent.com/PrisonRodeo/GSERM-2018-git/master/Data/SegalVotes.csv"
temp<-getURL(URL)
Segal<-read.csv(textConnection(temp))
rm(temp,URL)

summary(Segal)

# Fixed Effects model:

SegalFE<-glmmboot(vote~warrant+house+person+business+car+us+
         except,data=Segal,family="binomial",
         cluster=justid)
summary(SegalFE)

# Random-Effects:

SegalRE<-glmmML(vote~warrant+house+person+business+car+us+
                except+justideo,data=Segal,family="binomial",
                cluster=justid)
summary(SegalRE)

######################################
# Event Counts:
#
# Get SFTF data:

URL<-"https://raw.githubusercontent.com/PrisonRodeo/GSERM-2018-git/master/Data/SFTF.csv"
temp<-getURL(URL)
SFTF<-read.csv(textConnection(temp))
rm(temp,URL)

summary(SFTF)
pdim(SFTF)

# Panel Tobit: 

# library(plm)
SFTF.panel<-pdata.frame(SFTF,i="countryid")
# library(censReg)
Tobit.panel<-censReg(SumEvents~POLITY+unuurbpc+poldurab+year,
                     data=SFTF.panel,method="BHHH")
summary(Tobit.panel)

# Random effects Poisson

# library(lme4)
Poisson.RE<-glmer(ciob~POLITY+unuurbpc+poldurab+I(year-1900)+
                    (1|countryid),data=SFTF,family="poisson")
summary(Poisson.RE)

# Alternative RE Poisson, using glmmML:

# library(glmmML)
Poisson.RE.alt<-glmmML(ciob~POLITY+unuurbpc+poldurab+I(year-1900),
                       data=SFTF,cluster=countryid,
                       family="poisson")
summary(Poisson.RE.alt)

# Fixed effects Poisson:

Poisson.FE<-glm(ciob~POLITY+unuurbpc+poldurab+I(year-1900)+
            as.factor(countryid),data=SFTF,family="poisson")
summary(Poisson.FE)
