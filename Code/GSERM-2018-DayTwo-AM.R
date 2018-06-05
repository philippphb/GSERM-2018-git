#########################################################
# GSERM - St. Gallen (2018)
#
# Heirarchical models...
#
########################################################
# Load packages (install as needed), set options:

library(RCurl)
library(RColorBrewer)
library(colorspace)
library(foreign)
library(psych)
library(lme4)
library(plm)
library(gtools)
library(boot)
library(plyr)
library(texreg)
library(statmod)

options(scipen = 6) # bias against scientific notation
options(digits = 2) # show fewer decimal places

############################################
# "Robust" illustration:

id<-seq(1,100,1) # 100 observations
set.seed(7222009)
x<-rnorm(100) # N(0,1) noise
y<-1+1*x+rnorm(100)*abs(x)
library(rms)
fit<-ols(y~x,x=TRUE,y=TRUE)
fit

RVCV<-robcov(fit)
RVCV

bigID<-rep(id,16)
bigX<-rep(x,16)
bigY<-rep(y,16)
bigdata<-as.data.frame(cbind(bigID,bigY,bigX))
bigOLS<-ols(bigY~bigX,data=bigdata,x=TRUE,y=TRUE)
bigOLS

BigRVCV<-robcov(bigOLS,bigdata$bigID)
BigRVCV

# Plot:

xsim=c(-3,-2,-1,0,1,2,3)
hats100<-predict(fit,xsim,se.fit=TRUE)
hats1600<-predict(bigOLS,xsim,se.fit=TRUE)
hatsRVCV<-predict(BigRVCV,xsim,se.fit=TRUE)
ub100<-hats100$linear.predictors + (1.96*hats100$se.fit)
lb100<-hats100$linear.predictors - (1.96*hats100$se.fit)
ub1600<-hats1600$linear.predictors + (1.96*hats1600$se.fit)
lb1600<-hats1600$linear.predictors - (1.96*hats1600$se.fit)
ubRVCV<-hatsRVCV$linear.predictors + (1.96*hatsRVCV$se.fit)
lbRVCV<-hatsRVCV$linear.predictors - (1.96*hatsRVCV$se.fit)

pdf("GSERM-Robustness.pdf",6,5)
par(mar=c(4,4,2,2))
plot(x,y,pch=16,xlab="X",ylab="Y")
abline(lm(y~x),lwd=3)
lines(xsim,ub100,lty=2,lwd=2)
lines(xsim,lb100,lty=2,lwd=2)
lines(xsim,ub1600,lty=2,lwd=2,col="red")
lines(xsim,lb1600,lty=2,lwd=2,col="red")
lines(xsim,ubRVCV,lty=3,lwd=2,col="green")
lines(xsim,lbRVCV,lty=3,lwd=2,col="green")
legend("topleft",legend=c("Black Dashed Line is N=100",
                           "Red Dashed Line is Naive N=1600",
                           "Green Dashed Line is Robust N=1600"),
       bty="n")
dev.off()

############################################
# HLMs
############################################
# Data:

temp<-getURL("https://raw.githubusercontent.com/PrisonRodeo/GSERM-2018-git/master/Data/HIVDeaths.csv")
HIV<-read.csv(text=temp, header=TRUE)
HIV<-HIV[ which(is.na(HIV$HIVDeathRate)==FALSE), ]
HIV$LnDeathPM <- log(HIV$HIVDeathRate*1000)
summary(HIV)

# Plot:
  
pdf("HIVDeaths.pdf",6,5)
par(mar=c(4,4,2,2))
with(HIV, plot(density(HIVDeathRate,na.rm=TRUE),
               xlim=c(0,8.5),lwd=3,main="",
               xlab="Value"))
with(HIV, lines(density(LnDeathPM,na.rm=TRUE),
                lwd=3,lty=2,col="red"))
legend("topright",bty="n",lwd=3,col=c("black","red"),
       legend=c("Deaths Per Thousand","ln(Deaths Per Million)"))
dev.off()

###############################################
# Models...
#
# OLS:

OLSfit<-with(HIV, lm(LnDeathPM~GDPLagK+GDPGrowthLag+
                       OPENLag+POLITYLag+POLITYSQLag+CivilWarDummy+
                       InterstateWarLag+BatDeaths1000Lag))
summary(OLSfit)

# Fixed Effects:

FEfit<-plm(LnDeathPM~GDPLagK+GDPGrowthLag+
             OPENLag+POLITYLag+POLITYSQLag+CivilWarDummy+
             InterstateWarLag+BatDeaths1000Lag,
           data=HIV,effect="individual", model="within",
           index=c("ISO3","year"))
summary(FEfit)

# Random effects via lmer:

REfit<-lmer(LnDeathPM~GDPLagK+GDPGrowthLag+OPENLag+
              POLITYLag+POLITYSQLag+CivilWarDummy+
              InterstateWarLag+BatDeaths1000Lag+(1|ISO3),
            data=HIV,REML=FALSE)
summary(REfit)

# HLMs / varying coefficient models:
#
# Varying (random) intercept and slope on GDP

HLMfit1<-lmer(LnDeathPM~GDPLagK+(GDPLagK|ISO3)+GDPGrowthLag+
                OPENLag+POLITYLag+POLITYSQLag+CivilWarDummy+
                InterstateWarLag+BatDeaths1000Lag,
              data=HIV,REML=FALSE)
summary(HLMfit1)

# Testing:

anova(REfit,HLMfit1)
VarCorr(HLMfit1)

# Coefficients:

Bs<-data.frame(coef(HLMfit1)[1])

head(Bs)
mean(Bs$ISO3.GDPLagK)
var(Bs$ISO3.GDPLagK)

pdf("GDPRandomBetas.pdf",6,5)
par(mar=c(4,4,2,2))
with(Bs, plot(density(ISO3.GDPLagK),lwd=3,
              main="",xlab="Beta Values"))
abline(v=mean(Bs$ISO3.GDPLagK),lty=2)
dev.off()

# Varying intercept and slope on trade openness:

HLMfit2<-lmer(LnDeathPM~GDPLagK+GDPGrowthLag+OPENLag+
                (GDPGrowthLag|ISO3)+POLITYLag+POLITYSQLag+CivilWarDummy+
                InterstateWarLag+BatDeaths1000Lag,
              data=HIV,REML=TRUE)
summary(HLMfit2)
