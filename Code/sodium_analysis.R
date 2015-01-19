### Last edited: January 18, 2015
### Kellie Ottoboni
### Test for association between change in Na+ consumption and change in life expectancy, from 1990 to 2010

library(dplyr)
library(Hmisc)

setwd("~Kellie/Dropbox/Causal Inference/Salt/")
salt <- read.table("omnibus_data.csv", sep = "\t", header = TRUE)


### Clean up data, impute missing predictors within each year
salt$alcohol_M <- ifelse(is.na(salt$etohM), salt$etohboth, salt$etohM)
salt$alcohol_F <- ifelse(is.na(salt$etohF), salt$etohboth, salt$etohF)

salt <- arrange(salt, country, year)
salt_filt <- select(salt, country, e0_M:rgdpe,emp:alcohol_F)
salt_filt <- select(salt_filt, -(etohboth:etohF), -pwt_year)

salt2010 <- filter(salt_filt, year == 2010)
salt2010$avh       <- as.numeric(impute(salt2010$avh, median))
salt2010$hc        <- as.numeric(impute(salt2010$hc, median))
salt2010$ctfp      <- as.numeric(impute(salt2010$ctfp, median))
salt2010$rtfpna    <- as.numeric(impute(salt2010$rtfpna, median))
salt2010$alcohol_M <- as.numeric(impute(salt2010$alcohol_M, median))
salt2010$alcohol_F <- as.numeric(impute(salt2010$alcohol_F, median))

salt1990 <- filter(salt_filt, year == 1990)
salt1990$avh       <- as.numeric(impute(salt1990$avh, median))
salt1990$hc        <- as.numeric(impute(salt1990$hc, median))
salt1990$ctfp      <- as.numeric(impute(salt1990$ctfp, median))
salt1990$rtfpna    <- as.numeric(impute(salt1990$rtfpna, median))
salt1990$alcohol_M <- as.numeric(impute(salt1990$alcohol_M, median))
salt1990$alcohol_F <- as.numeric(impute(salt1990$alcohol_F, median))

### Take difference between 2010 and 1990, split into male and female.
countries <- salt2010$country
salt_diff <- salt2010[,-1]-salt1990[,-1]
salt_diff <- cbind(countries, salt_diff)
salt_diff <- filter(salt_diff, !is.na(e0_M))

male <- select(salt_diff, -countries, -e0_F, -year, -Na_F, -Na_M, -alcohol_F)
male_Na <- salt_diff$Na_M
colnames(male) <- gsub("_M", "", colnames(male))
female <- select(salt_diff, -countries, -e0_M, -year, -Na_M, -Na_F, -alcohol_M)
female_Na <- salt_diff$Na_F
colnames(female) <- gsub("_F", "", colnames(female))




### Modeling step - a function to test a bunch of different models and pick the one with the best in-sample fit
library(rpart)
library(ipred)
library(randomForest)
library(e1071)
library(gam)
library(gbm)
library(nnet)
choose_model <- function(dat){
  mod_lm <- lm(e0~., data = dat); mod_lm$fitted
  mod_poly <- lm(dat$e0~poly(as.matrix(dat[,-1], degree = 2, raw = TRUE))); predict(mod_poly)
  mod_tree <- rpart(e0~., data=dat); predict(mod_tree)
  mod_bag <- bagging(e0~., data=dat); predict(mod_bag)
  mod_rf <- randomForest(e0~., data=dat); predict(mod_rf)
  mod_svm <- svm(e0~., data = dat)
  mod_gam <- gam(e0~., data = dat)
  mod_boost <- gbm(e0~., data = dat, bag.fraction = 1, train.fraction = 1)
  mod_nnet <- nnet(e0/6 ~ ., data = dat, size = 2)
  models <- list(mod_lm, mod_poly, mod_tree, mod_bag, mod_rf, mod_svm, mod_gam, mod_nnet, mod_boost)
  mse <- sapply(models[1:length(models)-1], function(x) mean((predict(x) - dat$e0)^2))
  mse <- c(mse, mean((predict(mod_nnet)*6 - dat$e0)^2))
  mse <- c(mse, mean((predict(mod_boost, dat, n.trees = 100) - dat$e0)^2))
  best <- which.min(mse)[1]
  return(models[[best]])
}

mod_M <- choose_model(male)
pred_M <- predict(mod_M, male)
mod_F <- choose_model(female)
pred_F <- predict(mod_F, female)


### Permutation test
permu_pearson <- function(prediction, outcome, treatment, iters){
  resid <- prediction-outcome
  pearson_r <- cor(resid, treatment)
  distr <- replicate(iters, {
    tr <- sample(treatment)
    cor(resid, tr)
  })
  pval <- c("p_upper" = sum(distr >= pearson_r)/iters,
            "p_lower" = sum(distr <= pearson_r)/iters,  
            "twosided" = sum(abs(distr) >= abs(pearson_r))/iters)	
  return(list("estimate" = pearson_r, "distr" = distr, "pvalue" = pval))
}

res_M <- permu_pearson(prediction = pred_M, outcome = male$e0, treatment = male_Na, iters = 10000)
res_F <- permu_pearson(prediction = pred_F, outcome = female$e0, treatment = female_Na, iters = 10000)


### Confidence intervals - by bootstrapping
detach("package:dplyr")
bootstrap_CI_pearson <- function(prediction, outcome, treatment, iters, alpha = 0.05, side = "both"){
  ### side must be a string: either "both", "upper", or "lower"
  resid <- prediction-outcome
  distr <- replicate(iters, {
    resid_boot <- sample(resid, length(resid), replace=TRUE)
    cor(resid_boot, treatment)
  })
  if(side == "both"){
    alpha <- alpha/2
    quant <- c(floor(length(distr)*alpha), length(distr)-floor(length(distr)*alpha))
    print(quant)}else{
      if(side == "upper"){
        quant <- length(distr) - floor(length(distr)*alpha)
        print(quant)}else{
          quant <- floor(length(distr)*alpha)
          print(quant)}
        }
  return(sort(distr)[quant])      
}
  

permu_pearson_shift <- function(prediction, outcome, treatment, iters, shift){
  resid <- prediction-outcome
  pearson_r <- cor(resid, treatment)
  distr <- replicate(iters, {
    tr <- sample(treatment)
    cor(resid, tr) + shift
  })
  pval <- c("p_upper" = sum(distr >= pearson_r)/iters,
            "p_lower" = sum(distr <= pearson_r)/iters,  
            "twosided" = sum(abs(distr) >= abs(pearson_r))/iters)  
  return(list("estimate" = pearson_r, "pvalue"=pval))
}



permu_CI_pearson <- function(prediction, outcome, treatment, iters, alpha = 0.05, side = "both"){
  ### side must be a string: either "both", "upper", or "lower"
  resid <- prediction-outcome
  if(side == "both"){
    d1 <- permu_CI_pearson(prediction, outcome, treatment, iters, alpha = alpha/2, side = "lower")
    d2 <- permu_CI_pearson(prediction, outcome, treatment, iters, alpha = alpha/2, side = "upper")
    return(c(d1,d2))
  }
  which_p <- which(c("lower", "upper", "both") %in% side)
  pval <- rep(1,3)
  d <- ifelse(side == "upper", 0.005, -0.005); shift <- d
  while(pval[which_p] >= alpha){
    res <- permu_pearson_shift(prediction, outcome, treatment, iters, shift = shift)
    pval <- res$pvalue; shift <- shift+d
    cat(shift-d, "\n", pval, "\n")
  }   
  return(res$estimate +  shift)
}

ci_M <- permu_CI_pearson(prediction = pred_M, outcome = male$e0, treatment = male_Na, iters = 10000, side = "both")
ci_F <- permu_CI_pearson(prediction = pred_F, outcome = female$e0, treatment = female_Na, iters = 10000, side = "both")

cat("MALE:\n", "Estimate:", res_M$estimate, "\nP-value", res_M$pvalue, "\nConfidence Interval", ci_M)
cat("FEMALE:\n", "Estimate:", res_F$estimate, "\nP-value", res_F$pvalue, "\nConfidence Interval", ci_F)




### Plots
library(ggplot2)

pdf("sodium_lifeexp.pdf")
dat <- data.frame("sex"     = c(rep("Male",nrow(male)), rep("Female",nrow(female))),
                  "sodium"  = as.numeric(c(male_Na, female_Na)),
                  "le"      = c(male$e0, female$e0),
                  "country" = rep(salt_diff$countries, 2),
                  "ex_mort" = c(pred_M-male$e0, pred_F-female$e0),
                  stringsAsFactors = FALSE)
qplot(sodium, le, data = dat, facets = sex~., colour = sex) + theme(legend.position = "none") + labs(x = "Sodium (g/day)", y = "Life Expectancy (years)", title = "Change from 1990 to 2010") + geom_text(aes(label=ifelse(le<0,paste(country),"")) , hjust=-0.2, size = 3.5)
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
qplot(sodium, ex_mort, data = dat, facets = sex~., colour = sex) + theme(legend.position = "none") + labs(x = "Sodium (g/day)", y = "Excess Mortality (years)", title = "Change in Sodium") + geom_text(aes(label=ifelse(ex_mort < -2, paste(country), "")), hjust=-0.2, size = 3.5)
dev.off()

pdf("lifeex_exmort.pdf")
qplot(le, ex_mort, data = dat, facets = sex~., colour = sex) + theme(legend.position = "none") + labs(x = "Life Expectancy (years)", y = "Excess Mortality (years)", title = "Change in Life Expectancy")  + geom_text(aes(label=ifelse(ex_mort < -2, paste(country), "")), hjust=1.1, size = 3.5)
dev.off()

