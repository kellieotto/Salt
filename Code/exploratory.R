### January 10, 2015
### Kellie Ottoboni
### Exploratory data analysis for Na+

setwd("../Data")
salt <- read.table("omnibus_data.csv", sep = "\t", header = TRUE)
dim(salt) # 74 29
str(salt)

### Clean up data
salt$year <- as.factor(salt$year)
salt$alcoholM <- ifelse(is.na(salt$etohM), salt$etohboth, salt$etohM)
salt$alcoholF <- ifelse(is.na(salt$etohF), salt$etohboth, salt$etohF)

### Split males and females into high/low sodium groups
hist(salt$Na_M)
clust2_M <- kmeans(salt$Na_M, centers = 2, nstart = 25)
ord_M <- order(salt$Na_M)
plot(salt$Na_M[ord_M], col = clust2_M$cluster[ord_M])
threshold_M <- salt$Na_M[ord_M][which(diff(clust2_M$cluster[ord_M]) != 0)]
clust2_F <- kmeans(salt$Na_F, centers = 2, nstart = 25)
ord_F <- order(salt$Na_F)
plot(salt$Na_F[ord_F], col = clust2_F$cluster[ord_F])
threshold_F <- salt$Na_F[ord_F][which(diff(clust2_F$cluster[ord_F]) != 0)]
salt$Tr_M <- salt$Na_M >= threshold_M
salt$Tr_F <- salt$Na_F >= threshold_F
sum(salt$Tr_M != salt$Tr_F) # 4 country-year observations with discordant male/female sodium

### Model step
library(ipred)
mod_M <- bagging(e0_M ~ year + alcoholM + pc_gdp + pop + rgdpe + emp + avh + hc + cgdpe + cgdpo + ck + ctfp + rgdpna + rkna + rtfpna + labsh + sqrt_pc_gdp + log_pc_gdp, data = salt, na.action=na.omit)
mod_F<- bagging(e0_F ~ year + alcoholF + pc_gdp + pop + rgdpe + emp + avh + hc + cgdpe + cgdpo + ck + ctfp + rgdpna + rkna + rtfpna + labsh + sqrt_pc_gdp + log_pc_gdp, data = salt, na.action=na.omit)
library(randomForest)
mod_M2 <- randomForest(e0_M ~ year + alcoholM + pc_gdp + pop + rgdpe + emp + avh + hc + cgdpe + cgdpo + ck + ctfp + rgdpna + rkna + rtfpna + labsh + sqrt_pc_gdp + log_pc_gdp, data = salt, na.action = na.omit)
mod_F2 <- randomForest(e0_F ~ year + alcoholF + pc_gdp + pop + rgdpe + emp + avh + hc + cgdpe + cgdpo + ck + ctfp + rgdpna + rkna + rtfpna + labsh + sqrt_pc_gdp + log_pc_gdp, data = salt, na.action = na.omit)

mean((salt$e0_M - predict(mod_M, salt))^2, na.rm=T) #[1] 3.291907
mean((salt$e0_M - predict(mod_M2, salt))^2, na.rm=T) #[1] 0.7062297 - RF has much better MSE but can't handle missing values. We should come back to this

pred_M <- predict(mod_M, salt); pred_F <- predict(mod_F, salt)



### ModelMatch step - males
library(devtools); install_github("kellieotto/ModelMatch/ModelMatch"); library(ModelMatch)
bad <- which(is.na(salt$e0_M)); dat <- salt[-bad,]
strata <- Matches(treatment=dat$Tr_M, prediction=pred_M[-bad])
res_M <- permu_test_mean(strata, prediction = pred_M[-bad], treatment = dat$Na_M, response = dat$e0_M)
ci_M <- permu_CI_mean(groups = strata, prediction = pred_M[-bad], response = dat$e0_M, treatment = dat$Na_M, side = "both", verbose = TRUE)

### ModelMatch step - females

bad <- which(is.na(salt$e0_F)); dat <- salt[-bad,]
strata <- Matches(treatment=dat$Tr_F, prediction=pred_F[-bad])
res_F <- permu_test_mean(strata, prediction = pred_F[-bad], treatment = dat$Na_F, response = dat$e0_F)
ci_F <- permu_CI_mean(groups = strata, prediction = pred_F[-bad], response = dat$e0_F, treatment = dat$Na_F, side = "both", verbose = TRUE, alpha = 0.05)


### Plot
par(mfrow = c(1,2))
hist(res_M$perm_distribution/length(strata), main = "Life Expectancy, Males", xlab = "Na+ Effect"); abline(v=res_M$diff_means, col = "red")
abline(v = ci_M[1], lty = 2); abline(v = ci_M[2], lty = 2)
hist(res_F$perm_distribution/length(strata), main = "Life Expectancy, Females", xlab = "Na+ Effect"); abline(v=res_F$diff_means, col = "red")
abline(v = ci_F[1], lty = 2); abline(v = ci_F[2], lty = 2)
par(xpd=TRUE)
legend("topright",lty = c(1,2), col = c("red","black"), legend = c("Estimate", "95% CI"))
