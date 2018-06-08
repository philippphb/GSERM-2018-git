#########################################################
# GSERM - St. Gallen (2018)
#
# Survival model extensions...
#
########################################################
# Load packages (install as needed), set options:


library(RCurl)
library(RColorBrewer)
library(colorspace)
library(foreign)
library(gtools)
library(boot)
library(plyr)
library(texreg)
library(statmod)
library(survival)
library(pscl)
library(smcure)
library(nltm)
library(coxme)

options(scipen = 6) # bias against scientific notation
options(digits = 2) # show fewer decimal places

###########################
# Cure Models...
#
# Plot:

t<-seq(0.1,50,by=0.1)
pi<-0.5 # P(cured) = 0.5
lambda<-0.1
S.exp<-exp(-(lambda*t))
S.mix<-pi+((1-pi)*(exp(-lambda*t)))
S.nomix<-exp(log(pi)*(1-(exp(-lambda*t))))

# Plot:

pdf("CureTypeS.pdf",6,5)
par(mar=c(4,4,2,2))
plot(t,S.exp,t="l",lwd=3,xlab="Time",ylab="Survival",
     main=expression(paste("Exponential Hazards with ",
          lambda," = 0.1 and ",pi," = 0.5")))
lines(t,S.mix,lwd=3,lty=2,col="red")
lines(t,S.nomix,lwd=3,lty=3,col="blue")
abline(h=0.5,lty=5,col="grey")
legend("topright",bty="n",lwd=3,lty=c(1,2,3),
       col=c("black","red","blue"),
       c("No Cured Fraction","Mixture Cure","Non-Mixture Cure"))
dev.off()

# Simulate some cured-fraction data:

set.seed(7222009)
X<-rnorm(500)
Z<-rbinom(500,1,0.5)
T<-rweibull(500,shape=1.2,scale=1/(exp(0.5+1*X)))
C<-rbinom(500,1,(0.4-0.3*Z)) # Z increases cure probability
S<-Surv(T,C)

pdf("CureSimKM.pdf",10,5)
par(mar=c(4,4,2,2))
plot(survfit(S~1),mark.time=FALSE,lwd=c(3,1,1),
     xlab="Time",ylab="Survival")
dev.off()

coxph(S~X)
coxph(S~X+Z)

data.cure<-cbind(X,Z,T,C)
data.cure<-data.frame(data.cure)
cure.fit<-smcure(S~X,cureform=~Z,data=data.cure,model="ph")

cure.hat<-predictsmcure(cure.fit,c(rep(mean(X),times=2)),
                        c(0,1),model="ph")
cure.pic<-plotpredictsmcure(cure.hat,type="S",model="ph",
            lwd=c(3,3))

# Plot:

pdf("CureHats.pdf",10,5)
par(mar=c(4,4,2,2))
plotpredictsmcure(cure.hat,type="S",model="ph",
                  lwd=c(3,3),main="Predicted Survival")
legend("topright",inset=0.04,bty="n",lty=c(1,2),cex=1.2,
      c("X at mean, Z = 0","X at mean, Z = 1"))
dev.off()


# Fitting real cure models...sorta.

CFURL<-"https://raw.githubusercontent.com/PrisonRodeo/GSERM-2018-git/master/Data/ceasefiresTC.csv"
temp<-getURL(CFURL)
CF<-read.csv(textConnection(temp))
rm(CFURL,temp)

CF<-CF[complete.cases(CF),]

CF.S<-Surv(CF$peace,CF$uncensor)

pdf("CF-KM.pdf",10,5)
par(mar=c(4,4,2,2))
plot(survfit(CF.S~1),mark.time=FALSE,lwd=c(3,1,1),
     xlab="Time (in months)",ylab="Survival")
dev.off()

CF.cox<-coxph(CF.S~tie+imposed+lndeaths+contig+onedem+twodem,
              data=CF,method="efron")
CF.cox

CF.cure1.fit<-smcure(CF.S~tie+lndeaths+imposed,
                    cureform=~contig,data=CF,model="ph",
                    link="logit",emmax=500)
# LNB:

LNBURL<-"https://raw.githubusercontent.com/PrisonRodeo/GSERM-2018-git/master/Data/LNB.csv"
temp<-getURL(LNBURL)
LNB<-read.csv(textConnection(temp))
rm(LNBURL,temp)

# Deal w/missing data + sort:

LNB<-LNB[complete.cases(LNB),]
LNB<-LNB[order(LNB$dyad,LNB$year),]

LNB.S<-Surv(LNB$count1-1,LNB$count1,LNB$buofmzmid)
LNB.altS<-Surv(LNB$count1,LNB$buofmzmid)

pdf("LNB-KM.pdf",10,5)
par(mar=c(4,4,2,2))
plot(survfit(LNB.S~1,id=dyad,data=LNB),mark.time=FALSE,lwd=c(3,1,1),
     xlab="Time (in years)",ylab="Survival")
legend("bottomright",inset=0.04,bty="n",cex=1.2,
       c("Kaplan-Meier Survival", "Estimate, Long et al. (2007)"))
dev.off()

LNB.cox<-coxph(LNB.S~relcap+major+jdem+border+wartime+s_wt_glo+
               medarb+noagg+arbcom+organ+milinst+cluster(dyad),
               data=LNB,method="breslow")
LNB.cox

# # Runs for a long time:
LNB.cure<-smcure(LNB.altS~relcap+major+jdem+border+wartime+s_wt_glo+
                   medarb+noagg+arbcom+organ+milinst,
                   cureform=~border,model="ph",data=LNB)

# # Also does not work:
# LNB.cure1<-nltm(LNB.S~relcap+major+jdem+border,
#                 nlt.model="PHPHC",data=LNB)

# Stata code for strsmix:
# 
# stset count1, id(episode) f(buofmzmid==1)
# gen h0=0
# strsmix major jdem border wartime, bhazard(h0) distribution(weibull) \\\
#    link(logistic) k1(relcap major jdem border wartime s_wt_glo medarb \\\
#    noagg arbcom organ milinst)

#################################################
# Frailty (random effects) models.
#
# Simulate some frailty-type data...

set.seed(7222009)
G<-1:40        # "groups"
F<-rnorm(40)   # frailties
data<-data.frame(cbind(G,F))
data<-data[rep(1:nrow(data),each=20),]
data$X<-rbinom(nrow(data),1,0.5)
data$T<-rexp(nrow(data),rate=exp(0+1*data$X+(2*data$F)))
data$C<-rbinom(nrow(data),1,0.5)
data<-data[order(data$F),]

S<-Surv(data$T,data$C)

Fcolors<-diverge_hcl(length(F))[rank(F)]

pdf("FrailtyKM.pdf",6,5)
par(mar=c(4,4,2,2))
plot(survfit(S~strata(data$G)),col=Fcolors,mark=20,
     xlab="ln(Time)",ylab="Survival",log="x",xlim=c(0.0001,100))
legend("bottomleft",bty="n",cex=0.9,inset=0,
       c("Reds indicate strata","with larger frailties;",
         "blues with smaller ones"))
dev.off()

cox.noF<-coxph(S~X,data=data)
summary(cox.noF)

weib.noF<-survreg(S~X,data=data,dist="weib")
summary(weib.noF)

cox.F<-coxph(S~X+frailty.gaussian(F),data=data)
summary(cox.F)

weib.F<-survreg(S~X+frailty.gaussian(F),data=data,dist="weib")
summary(weib.F)

# Predicted survival plot:
# 
# plot(survfit(cox.noF),log="x",mark.time=FALSE,lwd=c(3,1,1),
#      xlab="ln(Time)",ylab="Fitted Survival")
# lines(survfit(cox.F),col="red",log="x",mark.time=FALSE,lwd=c(3,1,1))

# Examples using leaders data...

leadURL<-"https://raw.githubusercontent.com/PrisonRodeo/GSERM-2018-git/master/Data/leaders.csv"
temp<-getURL(leadURL)
lead<-read.csv(textConnection(temp))
rm(temp)

lead<-lead[lead$year<2004,]

lead.S<-Surv(lead$tenstart,lead$tenure,lead$tenureend)

Rs<-as.matrix(lead[,13:17])
lead$region<-factor((Rs %*% 1:ncol(Rs))+1,
                    labels=c("NorthAm",colnames(Rs)))
rm(Rs)

lead.F<-coxph(lead.S~female*region+frailty.gamma(ccode),data=lead)
summary(lead.F)

pdf("leadHats.pdf",6,5)
par(mar=c(4,4,2,2))
plot(survfit(lead.F,se.fit=TRUE),conf.int=TRUE,mark.time=FALSE,
     log="x",lwd=c(2,1,1),col="red",xlab="Time (in days)",
     ylab="Survival")
lines(survfit(lead.S~1),conf.int=FALSE,col="black",
      mark.time=FALSE,lwd=2)
legend("bottomleft",bty="n",inset=0.04,lty=1,lwd=3,col=c("red","black"),
       c("Predicted Survival","Baseline (Univariate) Survival"))
dev.off()

# Mixed-effects

lead.coxME<-coxme(lead.S~female + (1 | ccode/female),data=lead)
lead.coxME

# Stratified vs. Frailty, etc.

pdf("lead-KM.pdf",6,5)
par(mar=c(4,4,2,2))
plot(survfit(lead.S~1,id=leadid,data=lead),mark.time=FALSE,lwd=c(3,1,1),
     xlab="Time (in days)",ylab="Survival",log="x")
dev.off()

# Plot strata by country

pdf("leadKMcountries.pdf",6,5)
par(mar=c(4,4,2,2))
plot(survfit(lead.S~strata(ccode),id=leadid,data=lead),
     col=brewer.pal(9,"Set1"),log="x",mark.time=FALSE,
     xlab="Time (in days)", ylab="Survival")
dev.off()

# Plot strata by region:

pdf("leadKMregions.pdf",6,5)
par(mar=c(4,4,2,2))
plot(survfit(lead.S~strata(region),id=leadid,data=lead),
     col=brewer.pal(6,"Set1"),lwd=2,log="x",mark.time=FALSE,
     xlab="Time (in days)", ylab="Survival")
legend("bottomleft",inset=0.02,bty="n",col=brewer.pal(6,"Set1"),
       c("N. America","Latin America","Europe","Africa","Asia",
         "Middle East"),lty=1,lwd=2)
dev.off()

lead.Fstrat<-coxph(lead.S~female*strata(region)+
                     frailty.gamma(ccode),data=lead)
summary(lead.Fstrat)

lead.stratCl<-coxph(lead.S~female*strata(region)+
                      cluster(ccode),data=lead)
summary(lead.stratCl)

lead.FstratCl<-coxph(lead.S~female*strata(region)+frailty.gamma(ccode)+
                       cluster(ccode),data=lead)
# boom


