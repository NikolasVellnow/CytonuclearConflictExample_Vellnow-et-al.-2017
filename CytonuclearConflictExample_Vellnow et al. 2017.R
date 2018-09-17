#### Start ############################################################
rm(list=ls())
library(lme4)
library(arm)
library(blmeco)
library(doBy)
library(pbkrtest)
library(ggplot2)
library(ggthemes)
library(FactoMineR)
library(pscl)
library(lmerTest)
library(plyr)
library(dplyr)
library(parallel)


# data <- read.delim("/Users/niko/Dropbox/Basel//Cytonuclear Conflict/3 Iso -excl.txt")
data <- read.delim("C:/Users/nikolas/Dropbox/Basel/Cytonuclear Conflict/3 Iso -excl.txt")

data$line.cross <- as.factor(data$X.lines)
data$line.cross <- revalue(data$line.cross, c("22x8"="8x22",
                                              "29x28"="28x29",
                                              "47x6"="6x47",
                                              "57x49"="49x57",
                                              "61x44"="44x61",
                                              "69x39"="39x69",
                                              "84x75"="75x84"))
data$cytotype <- as.factor(data$line)
data$cytotypedummy <- as.factor(data$linedummy)
data$line.cross.cytotype <- as.factor(data$Xlines...line)
data$pairID <- as.factor(paste(data$line.cross,data$cross.. , data$rep, sep = "/"))
data$parentID <- as.factor(paste(data$cytotype,data$cross.. , data$rep, sep = "/"))
data$line.rep <- as.factor(data$cross..)
data$cytotype.line.rep <- as.factor(paste(data$cytotype,data$cross.., sep = "/"))
data$cross.rep <- as.factor(data$rep)
data$ID <- as.factor(data$ID)
data$mtGroup <- as.factor(data$mtGroup)
data$a <- as.numeric(data$a)
data$t <- as.numeric(data$t)
data$sem <- as.numeric(data$sem)
data$o <- as.numeric(data$o)
data$sa <- as.numeric(data$t/(data$o+data$t))

dataEye <- subset(data,mean.e!="NA")

table(data$mtGroup, data$line.cross)
# 4 of the crosses have not more than 8 replicates per cytotype...
# namely: 1x13, 29x28, 35x81, 40x76

# How many "pairs" with only one parent? load "dplyr"
pairID_parentID <- data %>% group_by(pairID, parentID) %>% tally()

# How many replicates per pair?
pair_replicates <- data %>% group_by(pairID) %>% tally()
mean(pair_replicates$n)
# 3.128
range(pair_replicates$n)
# [1] 1 4

# How many offspring per mother?
d <- data %>% group_by(parentID) %>% tally()
mean(d$n)
# 1.622407
range(d$n)
# [1] 1 2

#### Exploratory post hoc test for ovary size ########################

# effect.sizes <- read.delim("/Users/niko/Dropbox/Basel//Cytonuclear Conflict/Effect sizes.txt")
effect.sizes <- read.delim("C:/Users/nikolas/Dropbox/Basel/Cytonuclear Conflict/Effect sizes.txt")
effect.sizes$w <- I(1/(effect.sizes$se.d.Hedge)^2)

effect.sizes <- subset(effect.sizes,Transform!="original")
ble <- subset(effect.sizes,trait=="o")

# Make nicer with colors for all traits
ggplot(effect.sizes, aes(x = line.cross, y = d.Hedge, color=trait)) +
  geom_point(aes(factor=trait),position=position_dodge(width=0.75), size=2.0)+
  theme(panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"))+
  geom_errorbar(aes(factor=trait, ymin = low.CI.d.Hedge, ymax = high.CI.d.Hedge), width=0.0, size=0.8,position=position_dodge(width=0.75))+
  geom_hline(yintercept=0)+
  xlab("Line cross")+
  ylab("Hedge's d (+/- 95% CI)")

# Make nicer with colors for ovary
ggplot(ble, aes(x = line.cross, y = d.Hedge, color=cross.site)) +
  geom_point(aes(factor=cross.site),position=position_dodge(width=1.0), size=4.0)+
  theme(panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"))+
  geom_errorbar(aes(factor=cross.site, ymin = low.CI.d.Hedge, ymax = high.CI.d.Hedge), width=0.0, size=1.0,position=position_dodge(width=1.0))+
  geom_hline(yintercept=0)+
  xlab("Line cross")+
  ylab("Hedge's d (+/- 95% CI)")

# Make nicer with symbols for ovary
g <-ggplot(ble, aes(x = line.cross, y = d.Hedge, shape=cross.site)) +
  geom_point(aes(factor=cross.site),position=position_dodge(width=1.0), size=5.0)+
  theme(panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"))+
  geom_errorbar(aes(factor=cross.site, ymin = low.CI.d.Hedge, ymax = high.CI.d.Hedge), width=0.0, size=1.0,position=position_dodge(width=1.0))+
  geom_hline(yintercept=0)+
  xlab("Line cross")+
  ylab("Hedge's d (+/- 95% CI)")

g+theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

plot(ble$d.Hedge~ble$crossed.pop)
plot(effect.sizes$d.Hedge~effect.sizes$crossed.pop)

### First try for post hoc model
posthoc1 <- lm(d.Hedge ~ crossed.pop, data=ble)
summary(posthoc1)
anova(posthoc1, test=F)

plot(posthoc1)

# weigh it by the inverse of the variance (error variance of Hedges D)
posthoc2 <- lm(d.Hedge ~ crossed.pop, data=ble, weights=w)
summary(posthoc2)
anova(posthoc2, test=F)

plot(posthoc2)

# test for distance between sample sites
plot(ble$d.Hedge~ble$distance.meters)
posthoc3 <- lm(d.Hedge ~ distance.meters, data=ble, weights=w)
summary(posthoc3)
anova(posthoc3, test=F)

# test within vs. among site
plot(ble$d.Hedge~ble$cross.site) # very different variance...
posthoc4 <- lm(d.Hedge ~ cross.site, data=ble, weights=w)
summary(posthoc4)
anova(posthoc4, test=F)

plot(posthoc4)


# weighted t-test
library(weights)
wtd.t.test(ble$d.Hedge, weight=ble$w, bootse=TRUE, bootn=10000)
#$test
# [1] "One Sample Weighted T-Test"

# $coefficients
# t.value           df      p.value 
# 4.044698959 14.000000000  0.001205698 

# $additional
# Difference        Mean Alternative    Std. Err 
# 0.5671370   0.5671370   0.0000000   0.1402174


ggplot(ble, aes(x = cross.site, y = d.Hedge)) +
  geom_point(size = 4.0) +
  theme(panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"))+
  geom_hline(yintercept=0)+
  xlab("Line cross")+
  ylab("Hedge's d (+/- 95% CI)")


# make lmm with weights for all traits (and line cross as random factor)

nc <- detectCores()
clus <- makeCluster(rep("localhost",nc))

posthoc5 <- lmer(d.Hedge ~ cross.site*trait + (1|line.cross), weights=w, REML=F, data=effect.sizes)
summary(posthoc5)
anova(posthoc5)

posthoc6 <- lmer(d.Hedge ~ cross.site + trait + (1|line.cross), weights=w, REML=F, data=effect.sizes)
summary(posthoc6)

# interaction effect
summary(PBmodcomp(posthoc5, posthoc6,  nsim = 500, ref = NULL, seed=NULL, cl = clus, details = 0))
# Parametric bootstrap test; time: 69.46 sec; samples: 10000 extremes: 5556;
# large : d.Hedge ~ cross.site * trait + (1 | line.cross)
# small : d.Hedge ~ cross.site + trait + (1 | line.cross)
# stat     df    ddf p.value
# PBtest   5.7698                0.5556
# Gamma    5.7698                0.5530
# Bartlett 4.0203 5.0000         0.5465
# F        1.1540 5.0000 2.3238  0.5090
# LRT      5.7698 5.0000         0.3293

posthoc7 <- lmer(d.Hedge ~ trait + (1|line.cross), weights=w, REML=F, data=effect.sizes)
summary(posthoc7)
# cross.site effect
summary(PBmodcomp(posthoc6, posthoc7,  nsim = 500, ref = NULL, seed=NULL, cl = clus, details = 0))
# Parametric bootstrap test; time: 63.40 sec; samples: 10000 extremes: 552;
# large : d.Hedge ~ cross.site + trait + (1 | line.cross)
# small : d.Hedge ~ trait + (1 | line.cross)
# stat     df    ddf p.value  
# PBtest   4.1247               0.05529 .
# Gamma    4.1247               0.05283 .
# Bartlett 3.6448 1.0000        0.05624 .
# F        4.1247 1.0000 17.188 0.05802 .
# LRT      4.1247 1.0000        0.04226 *

posthoc8 <- lmer(d.Hedge ~ cross.site + (1|line.cross), weights=w, REML=F, data=effect.sizes)
summary(posthoc8)
# trait effect
summary(PBmodcomp(posthoc6, posthoc8,  nsim = 500, ref = NULL, seed=NULL, cl = clus, details = 0))
# Parametric bootstrap test; time: 60.50 sec; samples: 10000 extremes: 95;
# large : d.Hedge ~ cross.site + trait + (1 | line.cross)
# small : d.Hedge ~ cross.site + (1 | line.cross)
# stat      df    ddf  p.value   
# PBtest   19.9415                0.009599 **
# Gamma    19.9415                0.009761 **
# Bartlett 15.1565  5.0000        0.009714 **
# F         3.9883  5.0000 2.3585 0.181984   
# LRT      19.9415  5.0000        0.001282 **


posthoc9 <- lmer(d.Hedge ~ 1 + (1|line.cross), weights=w, REML=F, data=effect.sizes)
summary(posthoc9)
posthoc9a <- lmer(d.Hedge ~ -1 + (1|line.cross), weights=w, REML=F, data=effect.sizes)
summary(posthoc9a)
anova(posthoc9a, posthoc9)


stopCluster(clus)



#### Testis size ##############################################################

g1 <-ggplot(data,aes(factor(line.cross),t))+
  geom_boxplot(aes(fill = factor(cytotypedummy)))+
  geom_point(aes(factor=cytotypedummy),position=position_dodge(width=0.75))

g1 + scale_fill_discrete(guide=FALSE)

model1 <- lmer(scale(log(t))~ 1+ (1|line.cross/cytotype/line.rep/rep) + (1|plate) , data=data)
summary(model1)

model2 <- lmer(scale(log(t))~ 1+ (1|line.cross/cytotype/line.rep) + (1|plate), data=data)
summary(model2)

summary(PBmodcomp(model1, model2,  nsim = 500, ref = NULL, seed=NULL, cl = NULL, details = 0))

model3 <- lmer(scale(log(t))~ 1+ (1|line.cross/cytotype) + (1|plate), data=data)
summary(model3)

summary(PBmodcomp(model2, model3,  nsim = 500, ref = NULL, seed=NULL, cl = NULL, details = 0))

model4 <- lmer(scale(log(t))~ 1+ (1|line.cross) + (1|plate), data=data)
summary(model4)

summary(PBmodcomp(model3, model4,  nsim = 500, ref = NULL, seed=NULL, cl = NULL, details = 0))

model5 <- lmer(scale(log(t))~ 1 + (1|plate), data=data)
summary(model5)

summary(PBmodcomp(model4, model5,  nsim = 500, ref = NULL, seed=NULL, cl = NULL, details = 0)) 

# Does it differ when I compare it to the full model

model6 <- lmer(scale(log(t))~ 1+ (1|cytotype/cross.rep/rep) + (1|plate) , data=data)
summary(model6)

summary(PBmodcomp(model1, model6,  nsim = 500, ref = NULL, seed=NULL, cl = NULL, details = 0))

# Results do not seem to be very different...

##### Ovary size ######################################################################

ggplot(data,aes(factor(line.cross),o))+
  geom_boxplot(aes(fill = factor(cytotypedummy)))+
  geom_point(aes(factor=cytotypedummy),position=position_dodge(width=0.75))

model1 <- lmer(scale(log(o))~ 1+ (1|line.cross/cytotype/line.rep/rep) + (1|plate), data=data)
summary(model1)

# checking model fit
par(mfrow=(c(2,2)))
qqnorm(resid(model1))
qqline(resid(model1))

qqnorm(ranef(model1)$line.cross[,1])
qqline(ranef(model1)$line.cross[,1])

plot(fitted(model1), resid(model1))
abline(h=0)

data$fitted <- fitted(model1)
plot(data$fitted, jitter(scale(log(data$o))))
abline(0,1)

par(mfrow=(c(1,1)))

# autocorrelation
par(mfrow=(c(1,2)))
acf(resid(model1))
acf(resid(model1), type="p")
par(mfrow=(c(1,1)))

# seems fine

# check all random effects
par(mfrow=(c(2,2)))
qqnorm(ranef(model1)$line.cross[,1])
qqline(ranef(model1)$line.cross[,1])

qqnorm(ranef(model1)$cytotype[,1])
qqline(ranef(model1)$cytotype[,1])

qqnorm(ranef(model1)$cross..[,1])
qqline(ranef(model1)$cross..[,1])

qqnorm(ranef(model1)$rep[,1])
qqline(ranef(model1)$rep[,1])

par(mfrow=(c(1,1)))

# all seems fine


model2 <- lmer(scale(log(o))~ 1+ (1|line.cross/cytotype/line.rep) + (1|plate), data=data)
summary(model2)
summary(PBmodcomp(model1, model2,  nsim = 500, ref = NULL, seed=NULL, cl = NULL, details = 0))

model3 <- lmer(scale(log(o))~ 1+ (1|line.cross/cytotype) + (1|plate), data=data)
summary(model3)
summary(PBmodcomp(model2, model3,  nsim = 500, ref = NULL, seed=NULL, cl = clus, details = 0))

model4 <- lmer(scale(log(o))~ 1+ (1|line.cross) + (1|plate), data=data)
summary(model4)
summary(PBmodcomp(model3, model4,  nsim = 500, ref = NULL, seed=NULL, cl = NULL, details = 3))

# Also significant when comparing full to full - cytotype model?

model1b <- lmer(scale(log(o))~ 1+ (1|line.cross/cross.rep/rep) + (1|plate), data=data)
summary(model1b)
summary(PBmodcomp(model1, model1b,  nsim = 500, ref = NULL, seed=NULL, cl = NULL, details = 0))

# What Cytotype is responsible for that? Post hoc exploration

model4b <- lmer(scale(log(o))~ cytotype + (1|line.cross) , data=data)
summary(model4b) # Don't know if this is correct
model4bb <- lmer(scale(log(o))~ 1 + (1|line.cross) , data=data)
summary(model4bb)
summary(PBmodcomp(model4b, model4bb,  nsim = 500, ref = NULL, seed=NULL, cl = NULL, details = 0))

# What mitochondrial sequence group is responsible for that? Post hoc exploration

modelmt <- lmer(scale(log(o))~ mtGroup + (1|line.cross) , data=subset(data, mtGroup!="NA"))
summary(modelmt)
modelmt2 <- lmer(scale(log(o))~ 1 + (1|line.cross) , data=subset(data, mtGroup!="NA"))
summary(modelmt2)
summary(PBmodcomp(modelmt, modelmt2,  nsim = 500, ref = NULL, seed=NULL, cl = NULL, details = 0))


# make subset with promising crosses
dat <- subset(data, line.cross== "46x67" | line.cross=="14x26" | line.cross=="22x8" | line.cross=="61x44" | line.cross=="69x39")
# model on that subset
model4c <- lmer(scale(log(o))~ 1 + (1|line.cross/cytotype) , data=dat)
summary(model4c)

model5c <- lmer(scale(log(o))~ 1 + (1|line.cross) , data=dat)
summary(model5c)
summary(PBmodcomp(model4c, model5c,  nsim = 500, ref = NULL, seed=NULL, cl = NULL, details = 0))
# This is of course likely to bring false positives and when the significance theshold is divided by 15 (# of crosses)
# 0.05/15 = 0.0033 it is not significant anymore (although that might be too conservative for an explorative study)

# make subsets with promising crosses and check it cross by cross
dat46 <- subset(data, line.cross== "46x67")
dat14 <- subset(data, line.cross=="14x26")
dat22 <- subset(data, line.cross=="22x8")
dat61 <- subset(data, line.cross=="61x44")
dat69 <- subset(data, line.cross=="69x39")

# Cross 46x67
model46a <- lmer(scale(log(o))~ cytotype + (1|cross..) + (1|rep) , data=dat46)
summary(model46a)
model46b <- lmer(scale(log(o))~ 1 + (1|cross..) + (1|rep) , data=dat46)
summary(model46b)
summary(PBmodcomp(model46a, model46b,  nsim = 500, ref = NULL, seed=NULL, cl = NULL, details = 0))

# Cross 14x26
model14a <- lmer(scale(log(o))~ cytotype + (1|cross..)+ (1|rep) , data=dat14)
summary(model14a)
model14b <- lmer(scale(log(o))~ 1 + (1|cross..)+ (1|rep) , data=dat14)
summary(model14b)
summary(PBmodcomp(model14a, model14b,  nsim = 500, ref = NULL, seed=NULL, cl = NULL, details = 0))

# Cross 22x8
model22a <- lmer(scale(log(o))~ cytotype + (1|cross..)+ (1|rep) , data=dat22)
summary(model22a)
model22b <- lmer(scale(log(o))~ 1 + (1|cross..)+ (1|rep) , data=dat22)
summary(model22b)
summary(PBmodcomp(model22a, model22b,  nsim = 500, ref = NULL, seed=NULL, cl = NULL, details = 0))

# Cross 61x44
model61a <- lmer(scale(log(o))~ cytotype + (1|cross..) + (1|rep) , data=dat61)
summary(model61a)
model61b <- lmer(scale(log(o))~ 1 + (1|cross..) + (1|rep) , data=dat61)
summary(model61b)
summary(PBmodcomp(model61a, model61b,  nsim = 500, ref = NULL, seed=NULL, cl = NULL, details = 0))

# Cross 69x39
model69a <- lmer(scale(log(o))~ cytotype + (1|cross..) + (1|rep) , data=dat69)
summary(model69a)
model69b <- lmer(scale(log(o))~ 1 + (1|cross..) + (1|rep) , data=dat69)
summary(model69b)
summary(PBmodcomp(model69a, model69b,  nsim = 500, ref = NULL, seed=NULL, cl = NULL, details = 0))

# 22x8 should probably not be looked at in further experiments

model5 <- lmer(scale(log(o))~ 1 + (1|plate) , data=data)
summary(model5)
summary(PBmodcomp(model4, model5,  nsim = 500, ref = NULL, seed=NULL, cl = NULL, details = 0))

####### Ovary with bodyarea as covariate #######################################################################################

model1 <- lmer(scale(log(o))~ log(a)+ (1|line.cross/cytotype/line.rep/rep) + (1|plate), data=data)
summary(model1)

# checking model fit
par(mfrow=(c(2,2)))
qqnorm(resid(model1))
qqline(resid(model1))

qqnorm(ranef(model1)$line.cross[,1])
qqline(ranef(model1)$line.cross[,1])

plot(fitted(model1), resid(model1))
abline(h=0)

data$fitted <- fitted(model1)
plot(data$fitted, jitter(scale(log(data$o))))
abline(0,1)

par(mfrow=(c(1,1)))

# autocorrelation
par(mfrow=(c(1,2)))
acf(resid(model1))
acf(resid(model1), type="p")
par(mfrow=(c(1,1)))

# seems fine

# check all random effects
par(mfrow=(c(2,2)))
qqnorm(ranef(model1)$line.cross[,1])
qqline(ranef(model1)$line.cross[,1])

qqnorm(ranef(model1)$cytotype[,1])
qqline(ranef(model1)$cytotype[,1])

qqnorm(ranef(model1)$cross..[,1])
qqline(ranef(model1)$cross..[,1])

qqnorm(ranef(model1)$rep[,1])
qqline(ranef(model1)$rep[,1])

par(mfrow=(c(1,1)))

# all seems fine

model2 <- lmer(scale(log(o))~ log(a)+ (1|line.cross/cytotype/line.rep) + (1|plate), data=data)
summary(model2)
summary(PBmodcomp(model1, model2,  nsim = 500, ref = NULL, seed=NULL, cl = NULL, details = 0))

model3 <- lmer(scale(log(o))~ log(a)+ (1|line.cross/cytotype) + (1|plate), data=data)
summary(model3)
summary(PBmodcomp(model2, model3,  nsim = 500, ref = NULL, seed=NULL, cl = NULL, details = 0))

model4 <- lmer(scale(log(o))~ log(a)+ (1|line.cross) + (1|plate), data=data)
summary(model4)
summary(PBmodcomp(model3, model4,  nsim = 500, ref = NULL, seed=NULL, cl = NULL, details = 3))

model5 <- lmer(scale(log(o))~ log(a) + (1|plate) , data=data)
summary(model5)
summary(PBmodcomp(model4, model5,  nsim = 500, ref = NULL, seed=NULL, cl = NULL, details = 0))

# Also significant when comparing full to full - cytotype model?

model1b <- lmer(scale(log(o))~ log(a)+ (1|line.cross/cross.rep/rep) + (1|plate), data=data)
summary(model1b)
summary(PBmodcomp(model1, model1b,  nsim = 500, ref = NULL, seed=NULL, cl = NULL, details = 0))

