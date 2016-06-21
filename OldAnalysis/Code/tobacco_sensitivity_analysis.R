### Last edited: March 14, 2015
### Kellie Ottoboni
### Sensitivity analyses for smoking & excess mortality

set.seed(38)
library(dplyr)
library(Hmisc)

setwd("../Data")
salt <- read.table("omnibus_data.csv", sep = "\t", header = TRUE)
source("../Code/modeling.R")
library(ModelMatch)

### Clean up data, impute missing predictors within each year
salt <- mutate(salt, popF_prop = popF/(popF + popM), popM_prop = popM/(popF + popM))
salt$popF_prop[is.na(salt$popF)] <- 0.5; salt$popM_prop[is.na(salt$popM)] <- 0.5
salt$etohboth <- ifelse(is.na(salt$etohboth), salt$popM_prop*salt$etohM + salt$popF_prop*salt$etohF, salt$etohboth)

salt <- arrange(salt, country, year)
salt_filt <- dplyr::select(salt, country, e0_M, e0_F, year, Na_M, Na_F, etohboth, etohM, etohF, pc_gdp, popM_prop, popF_prop, smoking_M, smoking_F)

salt2010           <- filter(salt_filt, year == 2010)
etohM_mod <- lm(etohM~ e0_M + e0_F + Na_M + Na_F + pc_gdp + smoking_M + smoking_F, data = salt2010)
etohF_mod <- lm(etohF~ e0_M + e0_F + Na_M + Na_F + pc_gdp + smoking_M + smoking_F, data = salt2010)
twn <- which(is.na(salt2010$etohM))
salt2010[twn, "etohM"] <- predict(etohM_mod, salt2010[twn,]); salt2010[twn, "etohF"] <- predict(etohF_mod, salt2010[twn,])
salt2010[twn, "etohboth"] <- (salt2010$popM_prop*salt2010$etohM + salt2010$popF_prop*salt2010$etohF)[twn]

salt1990           <- filter(salt_filt, year == 1990)
etohboth_mod <- lm(etohboth ~ e0_M + e0_F + Na_M + Na_F + pc_gdp + smoking_M + smoking_F, data = salt1990)
twn <- which(is.na(salt1990$etohboth))
salt1990[twn, "etohboth"] <- predict(etohboth_mod, salt1990[twn,])
salt1990$etohM     <- salt1990$etohboth * (salt2010$etohM)/(salt2010$popM_prop*salt2010$etohM+salt2010$popF_prop*salt2010$etohF)
salt1990$etohF     <- salt1990$etohboth * (salt2010$etohF)/(salt2010$popM_prop*salt2010$etohM+salt2010$popF_prop*salt2010$etohF)


### Take difference between 2010 and 1990, split into male and female.
countries <- salt2010$country
salt_diff <- salt2010[,-1]-salt1990[,-1]
salt_diff <- cbind(countries, salt_diff)
salt_diff <- filter(salt_diff, !is.na(e0_M))



### original data
male <- dplyr::select(salt_diff, e0_M, etohM, pc_gdp)
male_smoking <- salt_diff$smoking_M
colnames(male) <- gsub("_M", "", colnames(male))
female <- dplyr::select(salt_diff, e0_F, etohF, pc_gdp)
female_smoking <- salt_diff$smoking_F
colnames(female) <- gsub("_F", "", colnames(female))

# predict on etoh and pc_gdp
mod_M1 <- choose_model(male)
pred_M1 <- predict(mod_M1, male)
mod_F1 <- choose_model(female)
pred_F1 <- predict(mod_F1, female)
# predict on etoh only
mod_M2 <- choose_model(male[,-3])
pred_M2 <- predict(mod_M2, male[,-3])
mod_F2 <- choose_model(female[,-3])
pred_F2 <- predict(mod_F2, female[,-3])
# predict on pc_gdp only
mod_M3 <- choose_model(male[,-2])
pred_M3 <- predict(mod_M3, male[,-2])
mod_F3 <- choose_model(female[,-2])
pred_F3 <- predict(mod_F3, female[,-2])

### Permutation test


res_M1 <- permu_pearson(prediction = pred_M1, response = male$e0, treatment = male_smoking, iters = 10000)
res_F1 <- permu_pearson(prediction = pred_F1, response = female$e0, treatment = female_smoking, iters = 10000)

res_M2 <- permu_pearson(prediction = pred_M2, response = male$e0, treatment = male_smoking, iters = 10000)
res_F2 <- permu_pearson(prediction = pred_F2, response = female$e0, treatment = female_smoking, iters = 10000)


res_M3 <- permu_pearson(prediction = pred_M3, response = male$e0, treatment = male_smoking, iters = 10000)
res_F3 <- permu_pearson(prediction = pred_F3, response = female$e0, treatment = female_smoking, iters = 10000)




res_table1 <- data.frame(rbind(c(res_M1$estimate, res_M1$pvalue), c(res_F1$estimate, res_F1$pvalue),c(res_M2$estimate, res_M2$pvalue), c(res_F2$estimate, res_F2$pvalue),c(res_M3$estimate, res_M3$pvalue), c(res_F3$estimate, res_F3$pvalue)))
rownames(res_table1) <- c("Male, both predictors", "Female, both predictors", "Male, ETOH", "Female, ETOH", "Male, pc_gdp", "Female, pc_gdp"); colnames(res_table1) <- c("Estimate", "Upper p", "Lower p", "Two-sided p")




### Subsampled data
subsamp <- sample(1:nrow(male), size = floor(1.5*nrow(male)), replace = TRUE)
male <- rbind(male, male[subsamp,] + matrix(ncol = 3, rnorm(3*length(subsamp), mean = 0, sd = .05)))
female <- rbind(female, female[subsamp,] + matrix(ncol = 3, rnorm(3*length(subsamp), mean = 0, sd = .05)))
female_smoking <- c(female_smoking, female_smoking[subsamp])
male_smoking<-c(male_smoking, male_smoking[subsamp])


# predict on etoh and pc_gdp
mod_M1 <- choose_model(male)
pred_M1 <- predict(mod_M1, male)
mod_F1 <- choose_model(female)
pred_F1 <- predict(mod_F1, female)
# predict on etoh only
mod_M2 <- choose_model(male[,-3])
pred_M2 <- predict(mod_M2, male[,-3])
mod_F2 <- choose_model(female[,-3])
pred_F2 <- predict(mod_F2, female[,-3])
# predict on pc_gdp only
mod_M3 <- choose_model(male[,-2])
pred_M3 <- predict(mod_M3, male[,-2])
mod_F3 <- choose_model(female[,-2])
pred_F3 <- predict(mod_F3, female[,-2])

### Permutation test


res_M1 <- permu_pearson(prediction = pred_M1, response = male$e0, treatment = male_smoking, iters = 10000)
res_F1 <- permu_pearson(prediction = pred_F1, response = female$e0, treatment = female_smoking, iters = 10000)

res_M2 <- permu_pearson(prediction = pred_M2, response = male$e0, treatment = male_smoking, iters = 10000)
res_F2 <- permu_pearson(prediction = pred_F2, response = female$e0, treatment = female_smoking, iters = 10000)


res_M3 <- permu_pearson(prediction = pred_M3, response = male$e0, treatment = male_smoking, iters = 10000)
res_F3 <- permu_pearson(prediction = pred_F3, response = female$e0, treatment = female_smoking, iters = 10000)




res_table2 <- data.frame(rbind(c(res_M1$estimate, res_M1$pvalue), c(res_F1$estimate, res_F1$pvalue),c(res_M2$estimate, res_M2$pvalue), c(res_F2$estimate, res_F2$pvalue),c(res_M3$estimate, res_M3$pvalue), c(res_F3$estimate, res_F3$pvalue)))
rownames(res_table2) <- c("Male, both predictors", "Female, both predictors", "Male, ETOH", "Female, ETOH", "Male, pc_gdp", "Female, pc_gdp"); colnames(res_table2) <- c("Estimate", "Upper p", "Lower p", "Two-sided p")


print(res_table1)
print(res_table2)

