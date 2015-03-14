# library(rpart)
# library(ipred)
# library(randomForest)
# library(e1071)
# library(gam)
# library(gbm)
# library(nnet)
library(caret)


# make it an option to return the info for all models
choose_model <- function(dat, bs = 10){
  ### dat = dataframe of inputs, including the outcome and all predictors. e0 is the outcome
  ### bs  = number of bootstrap replicates (default 10)
  suppressWarnings({
  ctrl <- trainControl(method = "boot", number = bs)
  mod_lm <- train(e0~., data = dat, method = "lm", trControl = ctrl, tuneLength = 10)
  mod_poly <- lm(dat$e0~poly(as.matrix(dat[,-1], degree = 2, raw = TRUE)))
  mod_tree <- train(e0~., data=dat, method = "rpart", trControl = ctrl, tuneLength = 10)
  mod_bag <- train(e0~., data = dat, method = "treebag", tuneLength = 10, trControl = ctrl)
  mod_rf <- train(e0~., data = dat, method = "rf", tuneLength =10, trControl = ctrl)
  mod_boost <- train(e0~., data = dat, method = "gbm", tuneLength = 20, bag.fraction = 1, trControl = ctrl, verbose = FALSE)
  mod_nnet <- train(e0/6~., data = dat, method = "nnet", tuneLength = 10, trControl = ctrl, trace = FALSE)
  })
  
  models <- list(mod_lm, mod_poly, mod_tree, mod_bag, mod_rf, mod_boost, mod_nnet)
  names(models) <- c("Linear regression", "Polynomial regression", "CART", "Bagged CART", "Random forest", "Stochastic gradient boosting", "Neural net")
  mse_calc <- function(x) mean((predict(x) - dat$e0)^2)
  mse <- sapply(models[1:length(models)-1], mse_calc)
  mse <- c(mse, mean((predict(mod_nnet)*6 - dat$e0)^2))
  names(mse) <- names(models)
  best <- which.min(mse)[1]
  print(mse)
  return(models[[best]])
}


