#########################################################
# GSERM - St. Gallen (2018)
#
# Generalized Estimating Equations.
#
########################################################
# Load packages (install as needed), set options:

library(RCurl)
library(lme4)
library(plm)
library(gtools)
library(plyr)
library(texreg)
library(statmod)
library(prais)
library(nlme)
library(pcse)
library(plm)
library(glmmML)
library(geepack)

options(scipen = 6) # bias against scientific notation
options(digits = 2) # show fewer decimal places

#######################
# Bush approval data:

url <- getURL("https://raw.githubusercontent.com/PrisonRodeo/GSERM-2018-git/master/Data/BushApproval.csv")
Bush <- read.csv(text = url) 
rm(url)

summary(Bush)
pdim(Bush)

# Independence:

GEE.IND<-geeglm(approval~partyid+perfin+nateco+age+educ+class+
                nonwhite+female,data=Bush,id=idno,family=gaussian,
                corstr="independence")
summary(GEE.IND)

# Compare to GLM:

GLM <- glm(approval~partyid+perfin+nateco+age+educ+class+
          nonwhite+female,data=Bush,family=gaussian)

# Coefficients:

cbind(GEE.IND$coefficients,GLM$coefficients)

# Standard Errors:
cbind(sqrt(diag(GEE.IND$geese$vbeta.naiv)),sqrt(diag(vcov(GLM))))

# Exchangeable correlation:

GEE.EXC<-geeglm(approval~partyid+perfin+nateco+age+educ+class+
                nonwhite+female,data=Bush,id=idno,family=gaussian,
                corstr="exchangeable")
summary(GEE.EXC)

# AR(a) correlation:

GEE.AR1<-geeglm(approval~partyid+perfin+nateco+age+educ+class+
                nonwhite+female,data=Bush,id=idno,family=gaussian,
                corstr="ar1")
summary(GEE.AR1)

# Unstructured correlation:

GEE.UNSTR<-geeglm(approval~partyid+perfin+nateco+age+educ+class+
                  nonwhite+female,data=Bush,id=idno,family=gaussian,
                  corstr="unstructured")
summary(GEE.UNSTR)

# Plot the betas / SEs:

library(car)
betas<-cbind(GEE.IND$coefficients,GEE.EXC$coefficients,
             GEE.AR1$coefficients,
             GEE.UNSTR$coefficients)

pdf("GEEBetasR.pdf",7,7)
par(mar=c(4,4,2,2))
scatterplotMatrix(betas[-1,],smooth=FALSE,
                  var.labels=c("Ind","Exch","AR1","Unstr"),
                  diagonal="none")
dev.off()

# SEs:

ses<-cbind(sqrt(diag(GEE.IND$geese$vbeta)),
             sqrt(diag(GEE.EXC$geese$vbeta)),
             sqrt(diag(GEE.AR1$geese$vbeta)),
             sqrt(diag(GEE.UNSTR$geese$vbeta)))

pdf("GEE-SEsR.pdf",7,7)
par(mar=c(4,4,2,2))
scatterplotMatrix(ses[-1,],smooth=FALSE,
                  var.labels=c("Ind","Exch","AR1","Unstr"),
                  diagonal="none")
dev.off()
