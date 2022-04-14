### PRS Analysis, PRS Figure, Age Distribution Figures, and APOE*Group Interaction Term Analysis
### Please contact Michael Osterman by email (mdo17@case.edu) with any questions


## Loading packages
library(ROCR)
library(pROC)
library(DescTools)
library(tidyverse)

## Read in phenotype and score data, clean for analysis
pheno = read.table("PHENO_Final.txt", header=T)

scores = read.table("PRS_Age_75.best", header=T)

pheno$PRS = scores$PRS
summary(pheno$PRS)
sd(pheno$PRS)
pheno$PRS = pheno$PRS - mean(pheno$PRS)
pheno$PRS = pheno$PRS * 1/sd(pheno$PRS)

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

Amish75_PRS_mod1 = glm(dx ~ PRS, data=Amish_age75)
summary(Amish75_PRS_mod1)
Cstat(Amish75_PRS_mod1, Amish_age75$dx)

Amish75_PRS_mod2 = glm(dx ~ sex + age + PRS, data=Amish_age75)
summary(Amish75_PRS_mod2)
Cstat(Amish75_PRS_mod2, Amish_age75$dx)

Amish75_PRS_mod3 = glm(dx ~ sex + age + PRS + apoe2 + apoe4, data=Amish_age75)
summary(Amish75_PRS_mod3)
Cstat(Amish75_PRS_mod3, Amish_age75$dx)
exp(coef(Amish75_PRS_mod3))
exp(confint(Amish75_PRS_mod3))

Amish75_PRS_mod4 = glm(dx ~ PRS + apoe2 + apoe4, data=Amish_age75)
summary(Amish75_PRS_mod4)
Cstat(Amish75_PRS_mod4, Amish_age75$dx)

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

EU75_PRS_mod1 = glm(dx ~ PRS, data=EU_age75)
summary(EU75_PRS_mod1)
Cstat(EU75_PRS_mod1, EU_age75$dx)

EU75_PRS_mod2 = glm(dx ~ sex + age + PRS, data=EU_age75)
summary(EU75_PRS_mod2)
Cstat(EU75_PRS_mod2, EU_age75$dx)

EU75_PRS_mod3 = glm(dx ~ sex + age + PRS + apoe2 + apoe4, data=EU_age75)
summary(EU75_PRS_mod3)
Cstat(EU75_PRS_mod3, EU_age75$dx)
exp(coef(EU75_PRS_mod3))
exp(confint(EU75_PRS_mod3))

EU75_PRS_mod4 = glm(dx ~ PRS + apoe2 + apoe4, data=EU_age75)
summary(EU75_PRS_mod4)
Cstat(EU75_PRS_mod4, EU_age75$dx)

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

Amish75_PRS_mod1 = glm(dx ~ PRS, data=Amish_age75)
summary(Amish75_PRS_mod1)
Cstat(Amish75_PRS_mod1, Amish_age75$dx)

Amish75_PRS_mod2 = glm(dx ~ sex + age + PRS, data=Amish_age75)
summary(Amish75_PRS_mod2)
Cstat(Amish75_PRS_mod2, Amish_age75$dx)

Amish75_PRS_mod3 = glm(dx ~ sex + age + PRS + relevel(apoe, ref = "3|3"), data=Amish_age75)
summary(Amish75_PRS_mod3)
Cstat(Amish75_PRS_mod3, Amish_age75$dx)
exp(coef(Amish75_PRS_mod3))
exp(confint(Amish75_PRS_mod3))

Amish75_PRS_mod4 = glm(dx ~ PRS + apoe, data=Amish_age75)
summary(Amish75_PRS_mod4)
Cstat(Amish75_PRS_mod4, Amish_age75$dx)

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

EU75_PRS_mod1 = glm(dx ~ PRS, data=EU_age75)
summary(EU75_PRS_mod1)
Cstat(EU75_PRS_mod1, EU_age75$dx)

EU75_PRS_mod2 = glm(dx ~ sex + age + PRS, data=EU_age75)
summary(EU75_PRS_mod2)
Cstat(EU75_PRS_mod2, EU_age75$dx)

EU75_PRS_mod3 = glm(dx ~ sex + age + PRS + relevel(apoe, ref = "3|3"), data=EU_age75)
summary(EU75_PRS_mod3)
Cstat(EU75_PRS_mod3, EU_age75$dx)
exp(coef(EU75_PRS_mod3))
exp(confint(EU75_PRS_mod3))

EU75_PRS_mod4 = glm(dx ~ PRS + apoe, data=EU_age75)
summary(EU75_PRS_mod4)
Cstat(EU75_PRS_mod4, EU_age75$dx)

### Comparing Risk Scores
Amish_cases_75 = subset(Amish_age75, dx==1)
Amish_controls_75 = subset(Amish_age75, dx==0)
EU_cases_75 = subset(EU_age75, dx==1)
EU_controls_75 = subset(EU_age75, dx==0)

t.test(Amish_cases_75$PRS, Amish_controls_75$PRS)
t.test(EU_cases_75$PRS, EU_controls_75$PRS)
t.test(Amish_cases_75$PRS, EU_cases_75$PRS)
t.test(EU_cases_75$PRS, Amish_controls_75$PRS)
t.test(Amish_cases_75$PRS, EU_controls_75$PRS)
t.test(Amish_controls_75$PRS, EU_controls_75$PRS)

### PRS Plot for Manuscript

Amish_cases_75$group = "Affected Amish"
Amish_controls_75$group = "Unaffected Amish"
EU_cases_75$group = "Non-Amish Cases"
EU_controls_75$group = "Non-Amish Controls"

pheno_75_combined = rbind(Amish_cases_75, EU_cases_75, Amish_controls_75, EU_controls_75)
pheno_75_combined$group = as.factor(pheno_75_combined$group)
pheno_75_combined$group = factor(pheno_75_combined$group, levels = c("Affected Amish", "Unaffected Amish", "Non-Amish Cases", "Non-Amish Controls"))

tiff("Figure2.tiff", units = "in", width = 7, height = 4, res = 300)
ggplot(pheno_75_combined, aes(x=group, y=PRS, fill=group)) +
  geom_violin() +
  geom_boxplot(width=.2) +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Spectral") +
  ggtitle("Distribution of Polygenic Risk Scores in Amish and Non-Amish Groups, Age 75+") +
  theme(plot.title = element_text(hjust=0.5), axis.text = element_text(size=12)) +
  xlab("Group") +
  ylab("Polygenic Risk Score") +
  scale_x_discrete(labels=c("Affected Amish \n n=137", "Unaffected Amish \n n=954", "Non-Amish Cases \n n=544", "Non-Amish Controls \n n=416"))
dev.off()


### Age Distribution Plots for Supplement

Amish_cases$group = "Affected Amish"
Amish_controls$group = "Unaffected Amish"
EU_cases$group = "Non-Amish Cases"
EU_controls$group = "Non-Amish Controls"
all_age = rbind(Amish_cases, Amish_controls, EU_cases, EU_controls)
table(all_age$group)
tiff("FigureS1.tiff", units = "in", width = 7, height = 4, res = 300)
ggplot(all_age, aes(x=group, y=age, fill=group)) +
  geom_violin() +
  geom_boxplot(width=0.2) +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Spectral") +
  ggtitle("Distribution of Age by Group and Alzheimer Disease Status") +
  theme(plot.title = element_text(hjust=0.5)) +
  xlab("Group") +
  ylab("Age (Years)") +
  scale_x_discrete(labels=c("Affected Amish \n n=152", "Unaffected Amish \n n=1690", "Non-Amish Cases \n n=1177", "Non-Amish Controls \n n=1126"))
dev.off()

tiff("FigureS2.tiff", units = "in", width = 7, height = 4, res = 300)
table(age_75$group)
ggplot(age_75, aes(x=group, y=age, fill=group)) +
  geom_violin() +
  geom_boxplot(width=0.2) +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Spectral") +
  ggtitle("Distribution of Age by Group and Alzheimer Disease Status, Age 75+") +
  theme(plot.title = element_text(hjust=0.5)) +
  xlab("Group") +
  ylab("Age (Years)") +
  scale_x_discrete(labels=c("Affected Amish \n n=137", "Unaffected Amish \n n=954", "Non-Amish Cases \n n=544", "Non-Amish Controls \n n=416"))
dev.off()

### APOE*group interaction term analysis

Combined_APOE_count_no_PRS_model = glm(dx ~ sex + age + apoe2 + apoe4 + apoe2*Amish + apoe4*Amish, data=pheno_75_combined)
summary(Combined_APOE_count_no_PRS_model)
Cstat(Combined_APOE_count_no_PRS_model, pheno$dx)
exp(coef(Combined_APOE_count_no_PRS_model))
exp(confint(Combined_APOE_count_no_PRS_model))
