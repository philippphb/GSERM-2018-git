#########################################################
# GSERM - St. Gallen (2018)
#
# GLS-ARMA Models + Dynamic Panel Data models...
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
library(prais)
library(nlme)
library(pcse)

options(scipen = 6) # bias against scientific notation
options(digits = 2) # show fewer decimal places

############################################
# Hall & Franzese (1998) Data:

temp<-getURL("https://raw.githubusercontent.com/PrisonRodeo/GSERM-2018-git/master/Data/HF1998.csv")
HF<-read.csv(text=temp, header=TRUE)
rm(temp)

# Summary statistics:

summary(HF)
HF$GDPPC <- HF$GDP_PC

# Nifty range plot by year:

Ann <- ddply(HF, .(year), summarise,
             MeanUE = mean(ue,na.rm=TRUE),
             MinUE = min(ue),
             MaxUE = max(ue))

pdf("HFMeanUE-R.pdf",6,5)
par(mar=c(4,4,2,2))
with(Ann, plot(year, MeanUE, t="l",lwd=2,ylim=c(0,18),
               xlab="Year",ylab="Unemployment Rate"))
with(Ann, points(year, MeanUE, pch=19))
with(Ann, segments(year,MinUE,year,MaxUE,lwd=2))
legend("topleft",bty="n",legend=c("Points Are OECD Means;",
                                  "lines indicate min-max ranges"))
dev.off()

# Models: OLS

HF.OLS<-plm(ue~GDPPC+open+uden+lcab+cbi+cwagebrg+wagexcbi,data=HF,model="pooling")
summary(HF.OLS)

library(lmtest)
coeftest(HF.OLS,vcov=vcovBK)

# Prais-Winsten:

HF.prais <- prais.winsten(ue~GDPPC+open+uden+lcab+cbi+cwagebrg+wagexcbi,
                          data=HF,iter=100)
HF.prais

# GLS with homoscedasticity & AR(1) autocorrelation:

HF.GLS <- gls(ue~GDPPC+open+uden+lcab+cbi+cwagebrg+wagexcbi,
                          HF,correlation=corAR1(form=~1|country))
summary(HF.GLS)

# GLS with unit-wise heteroscedasticity:

HF.GLS2 <- gls(ue~GDPPC+open+uden+lcab+cbi+cwagebrg+wagexcbi,
              HF,correlation=corAR1(form=~1|country),
              weights = varIdent(form = ~1|country))
summary(HF.GLS2)

# PCSEs
HF.lm<-lm(ue~GDPPC+open+uden+lcab+cbi+cwagebrg+wagexcbi,data=HF)
HF.pcse<-pcse(HF.lm,groupN = HF$country, groupT = HF$year)
summary(HF.pcse)

########################################
# Dynamic Models:

# Plot of long-run multiplier...

Phi <- seq(0,0.97,by=0.01)
Beta <- 1
BetaLR <- (Beta) / (1-Phi)

pdf("LDVLongRunPlotR.pdf",6,5)
par(mar=c(4,4,2,2))
plot(Phi,BetaLR,t="l",lwd=3,xlab=expression(phi),
     ylab="Long-Run Effect")
abline(v=c(0.5,0.8,0.9),lty=3)
dev.off()


###########################################
# Dynamics:

temp<-getURL("https://raw.githubusercontent.com/PrisonRodeo/GSERM-2018-git/master/Data/AIDS-Africa-97-01.csv")
AIDS<-read.csv(text=temp, header=TRUE)
rm(temp)

summary(AIDS)

# Annual mean / ranges plot:

HIVAnn <- ddply(AIDS, .(year), summarise,
             MeanHIV = mean(exp(lnAIDS),na.rm=TRUE),
             MinHIV = min(exp(lnAIDS)),
             MaxHIV = max(exp(lnAIDS)))

pdf("MeanHIVs-R.pdf",6,5)
par(mar=c(4,4,2,2))
with(HIVAnn, plot(year, MeanHIV, t="l",lwd=2,ylim=c(0,50),
               xlab="Year",ylab="HIV Rate"))
with(HIVAnn, points(year, MeanHIV, pch=19))
with(HIVAnn, segments(year,MinHIV,year,MaxHIV,lwd=2))
legend("topleft",bty="n",legend=c("Points Are All-Country Means;",
                                  "lines indicate min-max ranges"))
dev.off()

# Panel unit root tests:

lnAIDS<-cbind(AIDS$ccode,AIDS$year,AIDS$lnAIDS)
purtest(lnAIDS,exo="trend",test=c("levinlin"))
purtest(lnAIDS,exo="trend",test=c("hadri"))
purtest(lnAIDS,exo="trend",test=c("ips"))

# Models....

AIDS.lm<-lm(lnAIDS~lnAIDSlag+warlag+popden+refsin,data=AIDS)

AIDS.fe<-plm(lnAIDS~lnAIDSlag+warlag+popden+refsin,data=AIDS,effect="individual",model="within")

