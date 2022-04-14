### GRS Analysis, GRS Figure
### Please contact Michael Osterman by email (mdo17@case.edu) with any questions


## Loading packages
library(ROCR)
library(pROC)
library(DescTools)
library(tidyverse)

## Read in phenotype and score data, clean for analysis
pheno = read.table("PHENO_Final.txt", header=T)

scores = read.table("GRS_Age_75.best", header=T)

pheno$GRS = scores$GRS
summary(pheno$GRS)
sd(pheno$GRS)
pheno$GRS = pheno$GRS - mean(pheno$GRS)
pheno$GRS = pheno$GRS * 1/sd(pheno$GRS)

pheno$age = as.numeric(pheno$age)
pheno$dx = as.numeric(pheno$dx)
pheno$apoe = as.factor(pheno$apoe)

## Subsetting pheno file

EU_pheno = subset(pheno, Amish == 0)
Amish_pheno = subset(pheno, Amish == 1)

### Analysis with APOE allele counts

## Amish Age 75+ w/ APOE allele counts

Amish_age75 = subset(Amish_pheno, age>=75)

Amish75_Cov_mod = glm(dx ~ sex + age, data=Amish_age75)
summary(Amish75_Cov_mod)
Cstat(Amish75_Cov_mod, Amish_age75$dx)

Amish75_APOE_mod = glm(dx ~ apoe2 + apoe4, data=Amish_age75)
summary(Amish75_APOE_mod)
exp(coef(Amish75_APOE_mod))
exp(confint(Amish75_APOE_mod))
Cstat(Amish75_APOE_mod, Amish_age75$dx)

Amish75_Cov_APOE_mod = glm(dx ~ sex + age + apoe2 + apoe4, data=Amish_age75)
summary(Amish75_Cov_APOE_mod)
Cstat(Amish75_Cov_APOE_mod, Amish_age75$dx)

Amish75_GRS_mod1 = glm(dx ~ GRS, data=Amish_age75)
summary(Amish75_GRS_mod1)
Cstat(Amish75_GRS_mod1, Amish_age75$dx)

Amish75_GRS_mod2 = glm(dx ~ sex + age + GRS, data=Amish_age75)
summary(Amish75_GRS_mod2)
Cstat(Amish75_GRS_mod2, Amish_age75$dx)

Amish75_GRS_mod3 = glm(dx ~ sex + age + GRS + apoe2 + apoe4, data=Amish_age75)
summary(Amish75_GRS_mod3)
Cstat(Amish75_GRS_mod3, Amish_age75$dx)
exp(coef(Amish75_GRS_mod3))
exp(confint(Amish75_GRS_mod3))

Amish75_GRS_mod4 = glm(dx ~ GRS + apoe2 + apoe4, data=Amish_age75)
summary(Amish75_GRS_mod4)
Cstat(Amish75_GRS_mod4, Amish_age75$dx)

## EU Age 75+ w/ APOE allele counts

EU_age75 = subset(EU_pheno, age>=75)

EU75_Cov_mod = glm(dx ~ sex + age, data=EU_age75)
summary(EU75_Cov_mod)
Cstat(EU75_Cov_mod, EU_age75$dx)

EU75_APOE_mod = glm(dx ~ apoe2 + apoe4, data=EU_age75)
summary(EU75_APOE_mod)
Cstat(EU75_APOE_mod, EU_age75$dx)
exp(coef(EU75_APOE_mod))
exp(confint(EU75_APOE_mod))

EU75_Cov_APOE_mod = glm(dx ~ sex + age + apoe2 + apoe4, data=EU_age75)
summary(EU75_Cov_APOE_mod)
Cstat(EU75_Cov_APOE_mod, EU_age75$dx)

EU75_GRS_mod1 = glm(dx ~ GRS, data=EU_age75)
summary(EU75_GRS_mod1)
Cstat(EU75_GRS_mod1, EU_age75$dx)

EU75_GRS_mod2 = glm(dx ~ sex + age + GRS, data=EU_age75)
summary(EU75_GRS_mod2)
Cstat(EU75_GRS_mod2, EU_age75$dx)

EU75_GRS_mod3 = glm(dx ~ sex + age + GRS + apoe2 + apoe4, data=EU_age75)
summary(EU75_GRS_mod3)
Cstat(EU75_GRS_mod3, EU_age75$dx)
exp(coef(EU75_GRS_mod3))
exp(confint(EU75_GRS_mod3))

EU75_GRS_mod4 = glm(dx ~ GRS + apoe2 + apoe4, data=EU_age75)
summary(EU75_GRS_mod4)
Cstat(EU75_GRS_mod4, EU_age75$dx)

## Amish Age 75+ Analysis w/ APOE as factor

Amish_age75 = subset(Amish_pheno, age>=75)

Amish75_Cov_mod = glm(dx ~ sex + age, data=Amish_age75)
summary(Amish75_Cov_mod)
Cstat(Amish75_Cov_mod, Amish_age75$dx)

Amish75_APOE_mod = glm(dx ~ apoe, data=Amish_age75)
summary(Amish75_APOE_mod)
Cstat(Amish75_APOE_mod, Amish_age75$dx)

Amish75_Cov_APOE_mod = glm(dx ~ sex + age + apoe, data=Amish_age75)
summary(Amish75_Cov_APOE_mod)
Cstat(Amish75_Cov_APOE_mod, Amish_age75$dx)

Amish75_GRS_mod1 = glm(dx ~ GRS, data=Amish_age75)
summary(Amish75_GRS_mod1)
Cstat(Amish75_GRS_mod1, Amish_age75$dx)

Amish75_GRS_mod2 = glm(dx ~ sex + age + GRS, data=Amish_age75)
summary(Amish75_GRS_mod2)
Cstat(Amish75_GRS_mod2, Amish_age75$dx)

Amish75_GRS_mod3 = glm(dx ~ sex + age + GRS + relevel(apoe, ref = "3|3"), data=Amish_age75)
summary(Amish75_GRS_mod3)
Cstat(Amish75_GRS_mod3, Amish_age75$dx)
exp(coef(Amish75_GRS_mod3))
exp(confint(Amish75_GRS_mod3))

Amish75_GRS_mod4 = glm(dx ~ GRS + apoe, data=Amish_age75)
summary(Amish75_GRS_mod4)
Cstat(Amish75_GRS_mod4, Amish_age75$dx)

## EU Age 75+ Analysis

EU_age75 = subset(EU_pheno, age>=75)

EU75_Cov_mod = glm(dx ~ sex + age, data=EU_age75)
summary(EU75_Cov_mod)
Cstat(EU75_Cov_mod, EU_age75$dx)

EU75_APOE_mod = glm(dx ~ apoe, data=EU_age75)
summary(EU75_APOE_mod)
Cstat(EU75_APOE_mod, EU_age75$dx)

EU75_Cov_APOE_mod = glm(dx ~ sex + age + apoe, data=EU_age75)
summary(EU75_Cov_APOE_mod)
Cstat(EU75_Cov_APOE_mod, EU_age75$dx)

EU75_GRS_mod1 = glm(dx ~ GRS, data=EU_age75)
summary(EU75_GRS_mod1)
Cstat(EU75_GRS_mod1, EU_age75$dx)

EU75_GRS_mod2 = glm(dx ~ sex + age + GRS, data=EU_age75)
summary(EU75_GRS_mod2)
Cstat(EU75_GRS_mod2, EU_age75$dx)

EU75_GRS_mod3 = glm(dx ~ sex + age + GRS + relevel(apoe, ref = "3|3"), data=EU_age75)
summary(EU75_GRS_mod3)
Cstat(EU75_GRS_mod3, EU_age75$dx)
exp(coef(EU75_GRS_mod3))
exp(confint(EU75_GRS_mod3))

EU75_GRS_mod4 = glm(dx ~ GRS + apoe, data=EU_age75)
summary(EU75_GRS_mod4)
Cstat(EU75_GRS_mod4, EU_age75$dx)

### Comparing Risk Scores
Amish_cases_75 = subset(Amish_age75, dx==1)
Amish_controls_75 = subset(Amish_age75, dx==0)
EU_cases_75 = subset(EU_age75, dx==1)
EU_controls_75 = subset(EU_age75, dx==0)

t.test(Amish_cases_75$GRS, Amish_controls_75$GRS)
t.test(EU_cases_75$GRS, EU_controls_75$GRS)
t.test(Amish_cases_75$GRS, EU_cases_75$GRS)
t.test(EU_cases_75$GRS, Amish_controls_75$GRS)
t.test(Amish_cases_75$GRS, EU_controls_75$GRS)
t.test(Amish_controls_75$GRS, EU_controls_75$GRS)

### GRS Plot for Manuscript

Amish_cases_75$group = "Affected Amish"
Amish_controls_75$group = "Unaffected Amish"
EU_cases_75$group = "Non-Amish Cases"
EU_controls_75$group = "Non-Amish Controls"

pheno_75_combined = rbind(Amish_cases_75, EU_cases_75, Amish_controls_75, EU_controls_75)
pheno_75_combined$group = as.factor(pheno_75_combined$group)
pheno_75_combined$group = factor(pheno_75_combined$group, levels = c("Affected Amish", "Unaffected Amish", "Non-Amish Cases", "Non-Amish Controls"))

tiff("Figure1.tiff", units = "in", width = 7, height = 4, res = 300)
ggplot(pheno_75_combined, aes(x=group, y=GRS, fill=group)) +
  geom_violin() +
  geom_boxplot(width=.2) +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Spectral") +
  ggtitle("Distribution of Genetic Risk Scores in Amish and Non-Amish Groups, Age 75+") +
  theme(plot.title = element_text(hjust=0.5), axis.text = element_text(size=12)) +
  xlab("Group") +
  ylab("Genetic Risk Score") +
  scale_x_discrete(labels=c("Affected Amish \n n=137", "Unaffected Amish \n n=954", "Non-Amish Cases \n n=544", "Non-Amish Controls \n n=416"))
dev.off()
