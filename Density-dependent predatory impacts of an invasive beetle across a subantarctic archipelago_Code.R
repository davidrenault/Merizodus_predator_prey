########### Context Independent Alien Beetle ###############
library(lme4)
library(ggplot2)
library(MuMIn)
library(car)
library(glmmTMB)
library(pscl)
library(frair)
library(car)
library(glmmTMB)
library(devtools)
library(dplyr)
library(cowplot)
devtools::install_github("tomwenseleers/export")
library(export)


##### DatasetCall #####
REP <- read.csv2("E:/Storage/Dataset.csv", sep=";", dec=".", header = TRUE, stringsAsFactors = FALSE)

hist(REP$FemaleToMaleRatio)
hist(REP$MeanBodyMass)

REP$noneaten = as.numeric(as.character(REP$noneaten))
REP$alive = as.numeric(as.character(REP$alive))
REP$dead = as.numeric(as.character(REP$dead))
REP$attacked = as.numeric(as.character(REP$attacked))
REP$nonattacked = as.numeric(as.character(REP$nonattacked))
REP$Pop = as.factor(as.character(REP$Pop))
REP$NPrey = as.numeric(as.character(REP$NPrey))
REP$NPredateur = as.numeric(as.character(REP$NPredateur))
REP$Id = as.factor(as.character(REP$Id))


#Data visualisation
hist(REP$eaten)
hist((REP$eaten/REP$NPredator))
REP$PropEaten<-(REP$eaten/REP$NPrey)
hist(REP$PropEaten)
REP$PropAlive<-(REP$alive/REP$NPrey)
hist(REP$PropAlive)
REP$Propattacked<-(REP$attacked/REP$NPrey)
hist(REP$Propattacked)


hist(REP$ResidenceTime)
hist(log(REP$ResidenceTime))
hist(REP$ExpTime)
hist(log(REP$ExpTime))
hist(REP$NPrey)
hist(REP$NPredateur)
hist(REP$FemaleToMaleRatio)
hist(REP$MeanBodyMass)
REP$FemaleToMaleRatioPerc<-(REP$FemaleToMaleRatio*100)
hist(REP$FemaleToMaleRatioPerc)
hist(REP$MeanBodyMass)
#for the random factor to account for population id
REP$ExpTimeFactor = as.factor(as.character(REP$ExpTime))


#################################################################################################
######################################### Models ################################################
#################################################################################################


####### EATEN SIMPLE
#Random factors
EatenSimpleScaled <- glmer(cbind(eaten, noneaten)~ ((ResidenceTime + NPrey + I(NPrey^2) + NPredators +I(NPredators^2)+ ExpTime + MeanBodyMass + FemaleToMaleRatioPerc+NPredators*ResidenceTime+ResidenceTime*FemaleToMaleRatioPerc+ResidenceTime*MeanBodyMass))+(1|Id), data=REP, family = binomial(link="logit"), na.action = na.fail)
EatenSimpleGLM <- glm(cbind(eaten, noneaten)~ ((ResidenceTime + NPrey + I(NPrey^2) + NPredators+I(NPredators^2) + ExpTime + MeanBodyMass + FemaleToMaleRatioPerc+NPredators*ResidenceTime+ResidenceTime*FemaleToMaleRatioPerc+ResidenceTime*MeanBodyMass)), data=REP, family = binomial(link="logit"), na.action = na.fail)
anova(EatenSimpleScaled, EatenSimpleGLM)
#pop id as random factor should be kept


r.squaredGLMM(EatenSimpleScaled)
#best model selection
DredgeEaten<-dredge(EatenSimpleScaled, trace=2)
DredgeEaten
BestEatenSimple<-(get.models(DredgeEaten, 1)[[1]])
summary(BestEatenSimple)
Anova(BestEatenSimple)
r.squaredGLMM(BestEatenSimple)
hist(residuals(BestEatenSimple))

library(DHARMa)
sim2 <- simulateResiduals(BestEatenSimple)
testZeroInflation(sim2)
#no zero inflation
hist(REP$PropEaten)
testDispersion(sim2)
#no zero-inflation


########Dead Alive
aliveSimpleScaled <- glmer(cbind(dead, alive)~ ((ResidenceTime + NPrey + I(NPrey^2) + I(NPredators^2)+ NPredators + ExpTime + MeanBodyMass + FemaleToMaleRatioPerc+NPredators*ResidenceTime+ResidenceTime*FemaleToMaleRatioPerc+ResidenceTime*MeanBodyMass))+(1|Id), data=REP, family = binomial(link="logit"), na.action = na.fail)
aliveSimpleGLM <- glm(cbind(dead, alive)~ ((ResidenceTime + NPrey + I(NPrey^2) + NPredators+ I(NPredators^2) + ExpTime + MeanBodyMass + FemaleToMaleRatioPerc+NPredators*ResidenceTime+ResidenceTime*FemaleToMaleRatioPerc+ResidenceTime*MeanBodyMass)), data=REP, family = binomial(link="logit"), na.action = na.fail)

anova(aliveSimpleScaled, aliveSimpleGLM)
#random factor 


r.squaredGLMM(aliveSimpleScaled)
#selection of the best model
Dredgealive<-dredge(aliveSimpleScaled, trace=2)
Dredgealive
BestaliveSimple<-(get.models(Dredgealive, 1)[[1]])
summary(BestaliveSimple)
Anova(BestaliveSimple)
r.squaredGLMM(BestaliveSimple)
hist(residuals(BestaliveSimple))
#pas trop mal
library(DHARMa)
sim2 <- simulateResiduals(BestaliveSimple)
testZeroInflation(sim2) 
#no zero inflation



######### attacked / non attacked
attackedSimpleScaled <- glmer(cbind(attacked, nonattacked)~ ((ResidenceTime + NPrey + I(NPrey^2) + I(NPredators^2)+ NPredators + ExpTime + MeanBodyMass + FemaleToMaleRatioPerc+NPredators*ResidenceTime+ResidenceTime*FemaleToMaleRatioPerc+ResidenceTime*MeanBodyMass))+(1|Id), data=REP, family = binomial(link="logit"), na.action = na.fail)
attackedSimpleGLM <- glm(cbind(attacked, nonattacked)~ ((ResidenceTime + NPrey + I(NPrey^2)  + NPredators+I(NPredators^2) + ExpTime + MeanBodyMass + FemaleToMaleRatioPerc+NPredators*ResidenceTime+ResidenceTime*FemaleToMaleRatioPerc+ResidenceTime*MeanBodyMass)), data=REP, family = binomial(link="logit"), na.action = na.fail)

anova(attackedSimpleScaled, attackedSimpleGLM)
#random factor

r.squaredGLMM(attackedSimpleScaled)
summary(attackedSimpleScaled)
hist(residuals(attackedSimpleScaled))

#best model selection
Dredgeattacked<-dredge(attackedSimpleScaled, trace=2)
Dredgeattacked
BestattackedSimple<-(get.models(Dredgeattacked, 1)[[1]])
summary(BestattackedSimple)
Anova(BestattackedSimple)
r.squaredGLMM(BestattackedSimple)
hist(residuals(BestattackedSimple))
library(DHARMa)
sim2 <- simulateResiduals(BestattackedSimple)
testZeroInflation(sim2)
#no zero inflation

