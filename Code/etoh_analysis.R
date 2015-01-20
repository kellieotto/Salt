### Last edited: January 19, 2015
### Kellie Ottoboni
### Test for association between change in Na+ consumption and change in life expectancy, from 1990 to 2010

library(dplyr)
library(Hmisc)
detach(package:e1071)

setwd("../Data")
salt <- read.table("omnibus_data.csv", sep = "\t", header = TRUE)


### Clean up data, impute missing predictors within each year
salt$alcohol_M <- ifelse(is.na(salt$etohM), salt$etohboth, salt$etohM)
salt$alcohol_F <- ifelse(is.na(salt$etohF), salt$etohboth, salt$etohF)

salt <- arrange(salt, country, year)
salt_filt <- select(salt, country, e0_M:rgdpe,emp:alcohol_F)
salt_filt <- select(salt_filt, -(etohboth:etohF), -pwt_year)

salt2010           <- filter(salt_filt, year == 2010)
salt2010$avh       <- as.numeric(impute(salt2010$avh, median))
salt2010$hc        <- as.numeric(impute(salt2010$hc, median))
salt2010$ctfp      <- as.numeric(impute(salt2010$ctfp, median))
salt2010$rtfpna    <- as.numeric(impute(salt2010$rtfpna, median))


salt1990           <- filter(salt_filt, year == 1990)
salt1990$avh       <- as.numeric(impute(salt1990$avh, median))
salt1990$hc        <- as.numeric(impute(salt1990$hc, median))
salt1990$ctfp      <- as.numeric(impute(salt1990$ctfp, median))
salt1990$rtfpna    <- as.numeric(impute(salt1990$rtfpna, median))


### Take difference between 2010 and 1990, split into male and female.
countries <- salt2010$country
salt_diff <- salt2010[,-1]-salt1990[,-1]
salt_diff <- cbind(countries, salt_diff)
salt_diff <- filter(salt_diff, !is.na(e0_M) & !is.na(alcohol_M))

male <- select(salt_diff, -countries, -e0_F, -year, -Na_F, -alcohol_M, -alcohol_F)
male_etoh <- salt_diff$alcohol_M
colnames(male) <- gsub("_M", "", colnames(male))
female <- select(salt_diff, -countries, -e0_M, -year, -Na_M, -alcohol_F, -alcohol_M)
female_etoh <- salt_diff$alcohol_F
colnames(female) <- gsub("_F", "", colnames(female))



### Modeling step - test a bunch of different models and pick the one with the best in-sample fit

source("../Code/perm_tests.R")
mod_M <- choose_model(male)
pred_M <- predict(mod_M, male)
mod_F <- choose_model(female)
pred_F <- predict(mod_F, female)


### Permutation test

res_M <- permu_pearson(prediction = pred_M, outcome = male$e0, treatment = male_etoh, iters = 10000)
res_F <- permu_pearson(prediction = pred_F, outcome = female$e0, treatment = female_etoh, iters = 10000)


### Confidence intervals - by inverting the permutation test

ci_M <- permu_CI_pearson(prediction = pred_M, outcome = male$e0, treatment = male_etoh, iters = 10000, side = "both")
ci_F <- permu_CI_pearson(prediction = pred_F, outcome = female$e0, treatment = female_etoh, iters = 10000, side = "both")

cat("MALE:\n", "Estimate:", res_M$estimate, "\nP-value", res_M$pvalue, "\nConfidence Interval", ci_M)
cat("FEMALE:\n", "Estimate:", res_F$estimate, "\nP-value", res_F$pvalue, "\nConfidence Interval", ci_F)



### Plots
library(ggplot2)
setwd("../Analysis")

pdf("etoh_lifeexp.pdf")
dat <- data.frame("sex"     = c(rep("Male",nrow(male)), rep("Female",nrow(female))),
                  "etoh"  = as.numeric(c(male_etoh, female_etoh)),
                  "le"      = c(male$e0, female$e0),
                  "country" = rep(salt_diff$countries, 2),
                  "ex_mort" = c(pred_M-male$e0, pred_F-female$e0),
                  stringsAsFactors = FALSE)
qplot(etoh, le, data = dat, facets = sex~., colour = sex) + theme(legend.position = "none") + labs(x = "Alcohol Consumed (L/year)", y = "Life Expectancy (years)", title = "Change from 1990 to 2010") + geom_text(aes(label=ifelse(le<0,paste(country),"")) , hjust=-0.2, size = 3.5)
dev.off()


pdf("permdistribution_etoh.pdf")
dat2 <- data.frame("sex"     = c(rep("Male",length(res_M$distr)), rep("Female",length(res_F$distr))),
                   "distr"  = as.numeric(c(res_M$distr, res_F$distr)),
                   stringsAsFactors = FALSE)
vline.dat <- data.frame("sex" = c(rep("Male",3), rep("Female",3)),
                        "vline"  = as.numeric(c(res_M$estimate, ci_M, res_F$estimate, ci_F)),
                        "lty"   = rep(c("dashed", "dotted", "dotted"), 2),
                        stringsAsFactors = FALSE)
qplot(distr, data = dat2, geom = "histogram", facets = sex~., colour = sex, fill = sex) + geom_vline(aes(xintercept=vline, linetype=lty), data = vline.dat, show_guide = TRUE, color="black")+ theme(legend.position = "none") + labs(x = "Pearson Correlation", y = "", title = "Permutation Distribution")
dev.off()

pdf("etoh_exmort.pdf")
qplot(etoh, ex_mort, data = dat, facets = sex~., colour = sex) + theme(legend.position = "none") + labs(x = "Ethanol (L/year)", y = "Excess Mortality (years)", title = "Change in Alcohol Consumption") + geom_text(aes(label=ifelse(ex_mort < -2 | ex_mort > 2, paste(country), "")), hjust=-0.2, size = 3.5)
dev.off()

pdf("lifeex_exmort_etohmodels.pdf")
qplot(le, ex_mort, data = dat, facets = sex~., colour = sex) + theme(legend.position = "none") + labs(x = "Life Expectancy (years)", y = "Excess Mortality (years)", title = "Change in Life Expectancy")  + geom_text(aes(label=ifelse(ex_mort < -2 | ex_mort > 2, paste(country), "")), hjust=1.1, size = 3.5)
dev.off()



