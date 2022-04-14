### CI as Outcome in Amish Only Analysis for Supplement
### Please contact Michael Osterman by email (mdo17@case.edu) with any questions

## Loading packages
library(ROCR)
library(pROC)
library(DescTools)
library(tidyverse)

## Read in phenotype and score data, clean for analysis

pheno = read.table("PHENO_CI.txt", header=T)

Jansen = read.table("Jansen Results/Amish_CI_Jansen.best", header=T)
pheno$Jansen_PRS = Jansen$PRS

pheno$age = as.numeric(pheno$age)
pheno$dx = as.numeric(pheno$dx)
pheno$apoe = as.factor(pheno$apoe)

pheno$Jansen_PRS = pheno$Jansen_PRS - mean(pheno$Jansen_PRS)
pheno$Jansen_PRS = pheno$Jansen_PRS * 1/sd(pheno$Jansen_PRS)

## Amish Age 75+

Amish_age75 = subset(pheno, age>=75)

## Buiding Models with and without Covariates

Amish_Cov_mod = glm(dx ~ sex + age, data=Amish_age75)
summary(Amish_Cov_mod)
Cstat(Amish_Cov_mod, Amish_age75$dx)

Amish_APOE_mod = glm(dx ~ apoe, data=Amish_age75)
summary(Amish_APOE_mod)
Cstat(Amish_APOE_mod, Amish_age75$dx)

Amish_Cov_APOE_mod = glm(dx ~ sex + age + apoe, data=Amish_age75)
summary(Amish_Cov_APOE_mod)
Cstat(Amish_Cov_APOE_mod, Amish_age75$dx)

Amish_Jansen_PRS_mod1 = glm(dx ~ Jansen_PRS, data=Amish_age75)
summary(Amish_Jansen_PRS_mod1)
Cstat(Amish_Jansen_PRS_mod1, Amish_age75$dx)

Amish_Jansen_PRS_mod2 = glm(dx ~ sex + age + Jansen_PRS, data=Amish_age75)
summary(Amish_Jansen_PRS_mod2)
Cstat(Amish_Jansen_PRS_mod2, Amish_age75$dx)

Amish_Jansen_PRS_mod3 = glm(dx ~ sex + age + Jansen_PRS + apoe, data=Amish_age75)
summary(Amish_Jansen_PRS_mod3)
Cstat(Amish_Jansen_PRS_mod3, Amish_age75$dx)

Amish_Jansen_PRS_mod4 = glm(dx ~ Jansen_PRS + apoe, data=Amish_age75)
summary(Amish_Jansen_PRS_mod4)
Cstat(Amish_Jansen_PRS_mod4, Amish_age75$dx)

## Comparing Risk Scores in Age 75+

Amish_cases_75 = subset(Amish_age75, dx==1)
Amish_controls_75 = subset(Amish_age75, dx==0)
boxplot(Amish_cases_75$Kunkle_PRS, Amish_controls_75$Kunkle_PRS, names = c("Amish Cases", "Amish Controls"), main = "Comparison of Kunkle PRS Estimates by Cognitive Impairment Status, Age 75+")
boxplot(Amish_cases_75$Jansen_PRS, Amish_controls_75$Jansen_PRS, names = c("Amish Cases", "Amish Controls"), main = "Comparison of Jansen PRS Estimates by Cognitive Impairment Status, Age 75+")
plot(Amish_age75$Kunkle_PRS, Amish_age75$Jansen_PRS)
cor(Amish_age75$Kunkle_PRS, Amish_age75$Jansen_PRS)

t.test(Amish_cases_75$Kunkle_PRS, Amish_controls_75$Kunkle_PRS)
t.test(Amish_cases_75$Jansen_PRS, Amish_controls_75$Jansen_PRS)

## Creating figure for supplement

pheno$group = as.factor(pheno$dx)
pheno_complete = subset(pheno, !is.na(dx))
pheno_complete$group = factor(pheno_complete$group, levels = c("1", "0"))

tiff("FigureS3.tiff", units = "in", width = 7, height = 4, res = 300)
ggplot(pheno_complete, aes(x=group, y=Jansen_PRS, fill=group)) +
  geom_violin() +
  geom_boxplot(width=0.2) +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Spectral") +
  ggtitle("Distribution of Polygenic Risk Scores in Amish by CI Status, Age 75+") +
  theme(plot.title = element_text(hjust=0.5)) +
  xlab("Group") +
  ylab("Polygenic Risk Score") +
  scale_x_discrete(labels=c("CI Amish \n n=345", "CU Amish \n n=656", "Non-Amish Cases \n n=544", "Non-Amish Controls \n n=416"))
dev.off()

