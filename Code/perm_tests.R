library(rpart)
library(ipred)
library(randomForest)
library(e1071)
library(gam)
library(gbm)
library(nnet)
choose_model <- function(dat){
  ### dat = dataframe of inputs, including the outcome and all predictors. e0 is the outcome
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




bootstrap_CI_pearson <- function(prediction, outcome, treatment, iters, alpha = 0.05, side = "both"){
  ### side must be a string: either "both", "upper", or "lower"
  resid <- prediction-outcome
  distr <- replicate(iters, {
    resid_boot <- base::sample(resid, length(resid), replace=TRUE)
    cor(resid_boot, treatment)
  })
  if(side == "both"){
    alpha <- alpha/2
    quant <- c(floor(length(distr)*alpha), length(distr)-floor(length(distr)*alpha))
    }else{
      if(side == "upper"){
        quant <- length(distr) - floor(length(distr)*alpha)
        }else{quant <- floor(length(distr)*alpha)}
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