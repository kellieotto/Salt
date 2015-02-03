# library(rpart)
# library(ipred)
# library(randomForest)
# library(e1071)
# library(gam)
# library(gbm)
# library(nnet)
library(caret)
choose_model <- function(dat, cv = 5){
  ### dat = dataframe of inputs, including the outcome and all predictors. e0 is the outcome
  ### cv  = number of folds for cross-validation (default 5)
  ctrl <- trainControl(method = "cv", number = cv)
  mod_lm <- train(e0~., data = dat, method = "leapSeq", trControl = ctrl, tuneLength = 10); predict(mod_lm, dat)
  mod_poly <- lm(dat$e0~poly(as.matrix(dat[,-1], degree = 2, raw = TRUE))); predict(mod_poly)
  mod_tree <- train(e0~., data=dat, method = "rpart", trControl = ctrl, tuneLength = 10)
  mod_bag <- train(e0~., data = dat, method = "treebag", tuneLength = 10, trControl = ctrl)
  mod_rf <- train(e0~., data = dat, method = "rf", tuneLength = 16, trControl = ctrl)
  mod_svm <- train(e0~., data = dat, method = "svmRadial", tuneGrid = expand.grid("C" = 2^seq(-2, 8),"sigma"=2*seq(.1,1, by = 1/10)), trControl = ctrl)
  mod_boost <- train(e0~., data = dat, method = "gbm", tuneLength = 20, bag.fraction = 1, trControl = ctrl, verbose = FALSE)
  mod_nnet <- train(e0~., data = dat, method = "nnet", tuneLength = 10, trControl = ctrl, verbose = FALSE)
  
  
  
  ### Can't seem to get GAMs to work?
  #mod_gam <- gam(formula(paste("e0~", paste(sapply(colnames(male), function(x) paste("s(",x,")")), collapse = "+"))), data = dat)
  #mod_gam <- train(e0~., data = dat, method = "gam", trControl = ctrl)
  
  
  
  models <- list(mod_lm, mod_poly, mod_tree, mod_bag, mod_rf, mod_svm, mod_gam, mod_nnet, mod_boost)
  mse <- sapply(models[1:length(models)], function(x) mean((predict(x) - dat$e0)^2))
  best <- which.min(mse)[1]
  return(models[[best]])
}


