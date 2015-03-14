### Last edited: March 14, 2015
### Kellie Ottoboni
### Test for association between change in tobacco use and change in life expectancy, from 1990 to 2010

set.seed(38)
library(dplyr)
library(Hmisc)

setwd("../Data")
salt <- read.table("omnibus_data.csv", sep = "\t", header = TRUE)


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

male <- dplyr::select(salt_diff, e0_M, etohM, pc_gdp)
male_smoking <- salt_diff$smoking_M
colnames(male) <- gsub("_M", "", colnames(male))
female <- dplyr::select(salt_diff, e0_F, etohF, pc_gdp)
female_smoking <- salt_diff$smoking_F
colnames(female) <- gsub("_F", "", colnames(female))
# 
# subsamp <- sample(1:nrow(male), size = floor(3*nrow(male)), replace = TRUE)
# male <- rbind(male, male[subsamp,] + matrix(ncol = 3, rnorm(3*length(subsamp), mean = 0, sd = .05)))
# female <- rbind(female, female[subsamp,] + matrix(ncol = 3, rnorm(3*length(subsamp), mean = 0, sd = .05)))
# female_smoking <- c(female_smoking, female_smoking[subsamp])
# male_smoking<-c(male_smoking, male_smoking[subsamp])


### Modeling step - test a bunch of different models and pick the one with the best in-sample fit

source("../Code/modeling.R")
library(devtools); install_github("kellieotto/ModelMatch/ModelMatch"); library(ModelMatch)
mod_M <- choose_model(male)
pred_M <- predict(mod_M, male)
mod_F <- choose_model(female)
pred_F <- predict(mod_F, female)


### Permutation test


res_M <- permu_pearson(prediction = pred_M, response = male$e0, treatment = male_smoking, iters = 10000)
res_F <- permu_pearson(prediction = pred_F, response = female$e0, treatment = female_smoking, iters = 10000)



res_table <- data.frame(rbind(c(res_M$estimate, res_M$pvalue), c(res_F$estimate, res_F$pvalue)))
rownames(res_table) <- c("Male", "Female"); colnames(res_table) <- c("Estimate", "Upper p", "Lower p", "Two-sided p")
print(res_table)



### Plots
library(ggplot2)
setwd("../Analysis")

pdf("smoking_lifeexp.pdf")
dat <- data.frame("sex"     = c(rep("Male",nrow(male)), rep("Female",nrow(female))),
                  "smoking"  = as.numeric(c(male_smoking, female_smoking)),
                  "le"      = c(male$e0, female$e0),
                  "country" = rep(salt_diff$countries, 2),
                  "ex_mort" = c(pred_M-male$e0, pred_F-female$e0),
                  stringsAsFactors = FALSE)
qplot(smoking, le, data = dat, facets = sex~., colour = sex) + theme(legend.position = "none") + labs(x = expression(paste(Delta, " Smoking (cigarettes/smoker/day)")), y = expression(paste(Delta, " Life Expectancy (years)")), title = "Change from 1990 to 2010") + geom_text(aes(label=ifelse(le<0,paste(country),"")) , hjust=-0.2, size = 3.5)
dev.off()



pdf("smoking_exmort.pdf")
qplot(smoking, ex_mort, data = dat, facets = sex~., colour = sex) + theme(legend.position = "none") + labs(x = expression(paste(Delta, " Smoking (cigarettes/smoker/day)")), y = "Excess Mortality (years)", title = "Change from 1990 to 2010") + geom_text(aes(label=ifelse(abs(ex_mort) > 1, paste(country), "")), hjust=-0.2, size = 3.5)
dev.off()

pdf("lifeex_exmort_smokinganalysis.pdf")
qplot(le, ex_mort, data = dat, facets = sex~., colour = sex) + theme(legend.position = "none") + labs(x = expression(paste(Delta, " Life Expectancy (years)")), y = "Excess Mortality (years)", title = "Change in Life Expectancy")  + geom_text(aes(label=ifelse(abs(ex_mort) > 1, paste(country), "")), hjust=1.1, size = 3.5)
dev.off()


