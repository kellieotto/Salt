### Last edited: February 16, 2015
### Kellie Ottoboni
### Test for association between change in Na+ consumption and change in life expectancy, from 1990 to 2010

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
salt_filt <- dplyr::select(salt, country, e0_M, e0_F, year, Na_M, Na_F, etohboth, etohM, etohF, pc_gdp, popM_prop, popF_prop)

salt2010           <- filter(salt_filt, year == 2010)
etohM_mod <- lm(etohM~ e0_M + e0_F + Na_M + Na_F + pc_gdp, data = salt2010)
etohF_mod <- lm(etohF~ e0_M + e0_F + Na_M + Na_F + pc_gdp, data = salt2010)
twn <- which(is.na(salt2010$etohM))
salt2010[twn, "etohM"] <- predict(etohM_mod, salt2010[twn,]); salt2010[twn, "etohF"] <- predict(etohF_mod, salt2010[twn,])
salt2010[twn, "etohboth"] <- (salt2010$popM_prop*salt2010$etohM + salt2010$popF_prop*salt2010$etohF)[twn]

salt1990           <- filter(salt_filt, year == 1990)
etohboth_mod <- lm(etohboth ~ e0_M + e0_F + Na_M + Na_F + pc_gdp, data = salt1990)
twn <- which(is.na(salt1990$etohboth))
salt1990[twn, "etohboth"] <- predict(etohboth_mod, salt1990[twn,])
salt1990$etohM     <- salt1990$etohboth * (salt2010$etohM)/(salt2010$popM_prop*salt2010$etohM+salt2010$popF_prop*salt2010$etohF)
salt1990$etohF     <- salt1990$etohboth * (salt2010$etohF)/(salt2010$popM_prop*salt2010$etohM+salt2010$popF_prop*salt2010$etohF)


### Take difference between 2010 and 1990, split into male and female.
countries <- salt2010$country
salt_diff <- salt2010[,-1]-salt1990[,-1]
salt_diff <- cbind(countries, salt_diff)
salt_diff <- filter(salt_diff, !is.na(e0_M))

male <- dplyr::select(salt_diff, -countries, -e0_F, -year, -Na_F, -Na_M, -etohF, -etohboth)
male_Na <- salt_diff$Na_M
colnames(male) <- gsub("_M", "", colnames(male))
female <- dplyr::select(salt_diff, -countries, -e0_M, -year, -Na_M, -Na_F, -etohM, -etohboth)
female_Na <- salt_diff$Na_F
colnames(female) <- gsub("_F", "", colnames(female))




### Modeling step - test a bunch of different models and pick the one with the best in-sample fit

source("../Code/modeling.R")
library(devtools); install_github("kellieotto/ModelMatch/ModelMatch"); library(ModelMatch)
mod_M <- choose_model(male)
pred_M <- predict(mod_M, male)
mod_F <- choose_model(female)
pred_F <- predict(mod_F, female)


### Permutation test


res_M <- permu_pearson(prediction = pred_M, response = male$e0, treatment = male_Na, iters = 10000)
res_F <- permu_pearson(prediction = pred_F, response = female$e0, treatment = female_Na, iters = 10000)


### Confidence intervals - by inverting the permutation test



ci_M <- permu_CI_pearson(prediction = pred_M, response = male$e0, treatment = male_Na, iters = 10000, side = "both", verbosity = TRUE)
ci_F <- permu_CI_pearson(prediction = pred_F, response = female$e0, treatment = female_Na, iters = 10000, side = "both", verbosity = TRUE)

res_table <- data.frame(rbind(c(res_M$estimate, res_M$pvalue, ci_M), c(res_F$estimate, res_F$pvalue, ci_F)))
rownames(res_table) <- c("Male", "Female"); colnames(res_table) <- c("Estimate", "Upper p", "Lower p", "Two-sided p", "Lower CI", "Upper CI")
print(res_table)



### Plots
library(ggplot2)
setwd("../Analysis")

pdf("sodium_lifeexp.pdf")
dat <- data.frame("sex"     = c(rep("Male",nrow(male)), rep("Female",nrow(female))),
                  "sodium"  = as.numeric(c(male_Na, female_Na)),
                  "le"      = c(male$e0, female$e0),
                  "country" = rep(salt_diff$countries, 2),
                  "ex_mort" = c(pred_M-male$e0, pred_F-female$e0),
                  stringsAsFactors = FALSE)
qplot(sodium, le, data = dat, facets = sex~., colour = sex) + theme(legend.position = "none") + labs(x = expression(paste(Delta, " Sodium (g/day)")), y = expression(paste(Delta, " Life Expectancy (years)")), title = "Change from 1990 to 2010") + geom_text(aes(label=ifelse(le<0,paste(country),"")) , hjust=-0.2, size = 3.5)
dev.off()


pdf("permdistribution.pdf")
dat2 <- data.frame("sex"     = c(rep("Male",length(res_M$distr)), rep("Female",length(res_F$distr))),
                  "distr"  = as.numeric(c(res_M$distr, res_F$distr)),
                  stringsAsFactors = FALSE)
vline.dat <- data.frame("sex" = c(rep("Male",3), rep("Female",3)),
                        "vline"  = as.numeric(c(res_M$estimate, ci_M, res_F$estimate, ci_F)),
                        "lty"   = rep(c("dashed", "dotted", "dotted"), 2),
stringsAsFactors = FALSE)
qplot(distr, data = dat2, geom = "histogram", facets = sex~., colour = sex, fill = sex) + geom_vline(aes(xintercept=vline, linetype=lty), data = vline.dat, show_guide = TRUE, color="black")+ theme(legend.position = "none") + labs(x = "Pearson Correlation", y = "", title = "Permutation Distribution")
dev.off()

pdf("sodium_exmort.pdf")
qplot(sodium, ex_mort, data = dat, facets = sex~., colour = sex) + theme(legend.position = "none") + labs(x = expression(paste(Delta, " Sodium (g/day)")), y = "Excess Mortality (years)", title = "Change in Sodium") + geom_text(aes(label=ifelse(ex_mort < -2, paste(country), "")), hjust=-0.2, size = 3.5)
dev.off()

pdf("lifeex_exmort.pdf")
qplot(le, ex_mort, data = dat, facets = sex~., colour = sex) + theme(legend.position = "none") + labs(x = expression(paste(Delta, " Life Expectancy (years)")), y = "Excess Mortality (years)", title = "Change in Life Expectancy")  + geom_text(aes(label=ifelse(ex_mort < -2, paste(country), "")), hjust=1.1, size = 3.5)
dev.off()

