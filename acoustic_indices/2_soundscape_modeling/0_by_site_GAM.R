# Pre-modeling analysis of acoustic indices
# links to help:
# https://noamross.github.io/mgcv-esa-workshop/ and https://noamross.github.io/gams-in-r-course/
# 1) GAM basics https://gsp.humboldt.edu/olm/R/05_03_GAM.html
# 2) GAM coef simulations https://www.maths.ed.ac.uk/~swood34/mgcv/check-select.pdf
# 3) GAM book pdf in Mendeley"
#       - pg182 mgcv description of model fits
#       - pg328 describes some prediciton
# 4) GL + GAM general description: https://christophm.github.io/interpretable-ml-book/extend-lm.html#gam
# 5) GLM vs GAM opinion: https://stats.stackexchange.com/questions/380426/when-to-use-a-gam-vs-glm
# 6) some GAM implementation: https://m-clark.github.io/generalized-additive-models/application.html
# 7) posterior estimate with GAM: https://stats.stackexchange.com/questions/190348/can-i-use-bootstrapping-to-estimate-the-uncertainty-in-a-maximum-value-of-a-gam
# 8) https://kevintshoemaker.github.io/NRES-746/Generalized%20Additive%20Models%20(GAMs).pdf



# TODO:

library(data.table)
library(dplyr)
library(caret)
library(ggplot2)
library(mgcv)
library(tidyverse)
library(outliers)
library(gratia)

wd = '//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/paper1-AcousticIndices/results/'
# abgqi_df = fread(paste0(wd, 'ABGQI_inference/averages/site_avg_ABGQI.csv')) %>%
#   select(-contains('var')) %>%
#   select(site, wavs, contains('mean'))
# corrected ABGQI data (May 31 2022)
abgqi_df = fread(paste0(wd, 'paired_ARUs/corrected_data/corrected_site_abgqi_31May2022.csv'))
indices_df = fread(paste0(wd, 'acoustic_indices_aggregation/averages/site_avg_acoustic_indices.csv'))
df = merge(indices_df, abgqi_df, by = 'site')

# responses
responses = c("ACI", "ADI", "AEI", "NDSI", 'NDSI_A', 'NDSI_B', "BI", 
              "H", "Ht", "Hs", "M", 'R', "sfm", "rugo", 
              "zcr_max", "zcr_mean", "zcr_min")

predictors = c('Anthropophony', 'Biophony', 'Geophony', 
               'Quiet', 'Interference', 'Unidentified')

response_exam = function(df, response){
  gathered_df = df %>%
    gather(variable, value, -y)
  
  # y ~ x
  ggplot(gathered_df, aes(x = value, y = y)) +
    geom_point(alpha = 0.2) +
    facet_wrap(~variable, scales = 'free') +
    stat_smooth(method = 'gam', se = FALSE) +
    labs(title = paste0(response, 'with GAM curve', x = "Predictor value"))
}
toLogit<-function(x){
  cx<-ifelse(x==0,0.00001,
             ifelse(x==1,0.99999,x))
  lgx<-log(cx)-log(1-cx)
  return(lgx)
}


# Drop sites with any 0 values
mod_df = abgqi_df
mod_df[mod_df == 0] = NA
mod_df = mod_df[complete.cases(mod_df),] # n = -27
mod_df$ARU = substr(mod_df$site, 4,5)

# # Based on Wilcoxon tests and corrections we will remove LGs
# mod_df = mod_df %>%
#   filter(ARU == 'am')
# 
# # results in dropping:
# 1247 - 994 # 253 sites were LG
# sum(abgqi_df$wavs) - sum(mod_df$wavs) # 155,847 recordings
# 
# # AM has
# nrow(mod_df) # 994 sites
# sum(mod_df$wavs) #585,214 recordings



###########################################
# NDSI
i = responses[4]
temp_y = data.frame(site = indices_df$site, NDSI = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
   select(-c(site, ARU))

# visualize pairs of corr
temp_df %>%
  GGally::ggpairs(aes(alpha = 0.05), progress = FALSE)

# add transformed NDSI
# temp_df$logit_NDSI = toLogit((temp_df$NDSI + 1)/2)
temp_df$posNDSI = temp_df$NDSI + 1
temp_df$beta_NDSI = (temp_df$NDSI + 1)/2
temp_df$sq_beta_NDSI = (temp_df$beta_NDSI)^2

 
###### TESTING using: https://www.youtube.com/watch?v=Ukfvd8akfco
mod1 = gam(NDSI ~ 
              Anthropophony+ 
              Biophony + 
              Geophony + 
              Quiet + 
              Interference,
            data = temp_df)
summary(mod1)

mod2 = gam(beta_NDSI ~ 
            s(Anthropophony) + 
            s(Biophony) + 
            s(Geophony) + 
            s(Quiet) + 
            s(Interference),
          data = temp_df,
          family = betar(),
          method = 'REML')
summary(mod2) # model summary

AIC(mod1, mod2)


mod3 = gam(beta_NDSI ~ 
             s(Anthropophony) + 
             s(Biophony) + 
             s(Quiet) + 
             s(Interference),
           data = temp_df,
           family = betar(),
           method = 'REML')
summary(mod3)

AIC(mod2, mod3)

draw(mod3) # nice visualization tool from gratia for partial effects

k.check(mod3) 

#appraise(mod2, method = "simulate") 
# QQ should follow 1:1 wihtin the reference band (do data follow assumptions of model well) some do not
# residuals should be approx gaussian 
par(mfrow = c(2, 2))
vis.gam(mod3, theta = 65)

#appraise(mod8, method = "simulate") # only works with gaussian data
gam.check(mod3, rep = 500) # initial Anthro k is too small there is more unexplained information not captured
                           # large k doesn't allow for more modeling freedom
                           # for rep = 500 does the QQ roughly fall within confidence interval? Not at the upper tail
                           # histogram is roughly gaussian, Resp vs fitted is roughly 1:1, resids vs linear pred is potentially problematic

variance_comp(mod8)


# PREDICTIONS
library(itsadug)
pred = get_predictions(mod3, 
                       cond = list(Anthropophony = seq(min(temp_df$Anthropophony, na.rm = TRUE),
                                                            max(temp_df$Anthropophony, na.rm = TRUE),
                                                            length.out = nrow(temp_df)), se = TRUE))

pred = pred %>%
  mutate(NDSI = 1)

temp_df %>%
  ggplot(aes(x = Anthropophony, y = NDSI)) +
    geom_point(alpha = 0.3) +
  geom_ribbon(data = pred, alpha = 0.4, aes(ymin = fit - CI, ymax = fit + CI), show.legend = F, fill = 'forestgreen') +
  geom_line(data = pred, aes(y = fit), show.legend = FALSE, color = 'forestgreen')




###########################################
# ACI
i = responses[1]
temp_y = data.frame(site = indices_df$site, ACI = indices_df[[i]])
temp_df = merge(x = temp_y, y = mod_df) %>% 
  select(-c(site,wavs))

# visualize pairs of corr
temp_df %>%
  GGally::ggpairs(aes(alpha = 0.05))

# add transformed 
temp_df$log_ACI = log(temp_df$ACI)
temp_df$sqrt_ACI = sqrt(temp_df$ACI)
outlier(temp_df$ACI)
temp_df = temp_df %>%
  filter(ACI < outlier(ACI))

faraway::halfnorm(temp_df$ACI)

# fit a model
gaussian_fit = gam(ACI ~ 
                 s(Anthropophony_mean)+ s(Biophony_mean) + s(Geophony_mean) + s(Quiet_mean) + s(Interference_mean),
               data = temp_df,
               family = gaussian(),
               method = 'REML')

# gaussian_fit2 = gam(ACI ~ s(Anthropophony_mean)+ s(Biophony_mean) + s(Geophony_mean) + s(Quiet_mean) + s(Interference_mean) +
#                       te(Anthropophony_mean, Biophony_mean) + te(Anthropophony_mean, Quiet_mean) +
#                       te(Biophony_mean, Quiet_mean) + te(Geophony_mean, Interference_mean),
#                     data = temp_df,
#                     family = gaussian(),
#                     method = 'REML')

lmfit = lm(ACI ~ Anthropophony_mean + Biophony_mean + Geophony_mean + Quiet_mean + Interference_mean,
                            data = temp_df)

# DIAGNOSTICS
mod = gaussian_fit

#fit # dof = dof of each covariate + 1. A dof ~ 1 is approximately linear
plot(mod, residuals = TRUE, pch = 19, cex = 0.5, pages = 1, shade = T, scale = F)
summary(mod)
AIC(mod)
par(mfrow = c(2, 2))
gam.check(mod)
vis.gam(mod, theta = 65)
anova(mod)
























# generate predictions using covariate values
# https://www.maths.ed.ac.uk/~swood34/mgcv/check-select.pdf
y_hat = predict(mod, type = 'lpmatrix')

# simulate coefficients drawing from posterior
library(MASS)
eq = seq(min(boot_df[,1]), max(boot_df[,1]), length = 1000)
br = mvrnorm(1000, coef(fit), vcov(fit))
max.eq = rep(NA, 1000)
for (i in 1:1000){ 
  fv = y_hat %*% br[i,]
  max.eq[i] = eq[fv==max(fv)]
}

ci = quantile(max.eq[complete.cases(max.eq)], c(.025,.975))
ci

# PREDICTION OPTION
y_hat = predict.gam(fit, type = 'response')
preds = data.frame(y = temp_df[,1], y_hat)
ggplot(preds, aes(x = y_hat, y = y)) +
  geom_point() +
  geom_abline(slope=1, intercept=0)

pred = data.frame(
  Anthropophony_mean = temp_df$Anthropophony_mean,
  Biophony_mean = temp_df$Biophony_mean,
  Geophony_mean = temp_df$Geophony_mean,  
  Quiet_mean = temp_df$Quiet_mean,
  Interference_mean = temp_df$Interference_mean,
  Unidentified_mean = temp_df$Unidentified_mean,
  NDSI = temp_df$NDSI,
  predicted_values = predict(fit, newdata = temp_df)) %>%
  gather(variable, value, -predicted_values, -NDSI)
  
  
ggplot(pred, aes(x = value)) +
  geom_point(aes(y = NDSI), size = 1, alpha = 0.5) +
  geom_smooth(aes(y = predicted_values), colour = "red", alpha = 0.3) +
  facet_wrap(~variable)

# 
# # bootstrap resampling w/ replacement
# set.seed(143)
# n_boot = 3
# replicates = lapply(1:n_boot, function(x) temp_df[sample(nrow(temp_df), nrow(temp_df), replace = TRUE),])

# 
# summary(fit)
# gam.check(fit)
# 
# # diagnostics
# plot.gam(fit, residuals = TRUE, # residuals
#          shade = TRUE, # 95% CI for the mean shape of effects
#          shade.col = "lightblue", seWithMean = TRUE, pages = 1, pch =1, cex = 1, all.terms = TRUE)
# plot(fit1, pages = 1, all.terms = TRUE) # partial effects plots
# 
# 
# # small p-values indicate residuals are not randomly distributed (i.e., not enough basis fxs)
# gam.check(fit1)
# # if p is significant, try increasing k: gam(y ~ s(x1, k = 12))
# 
# # values over ~0.8 for worst may mean issues with collinearity
# concurvity(fit, full = TRUE)
# concurvity(fit, full = FALSE) # use to id where the relationships are UNID and Int is highest
# 
# AIC(fit) # 1074.22
# AIC(fit1) # 1056.504
# 
# # set up pred df
# temp_df$y_pred = predict(fit, newdata = temp_df)
# 
# temp_df %>%
#   gather(covariate, value, -NDSI, -y_pred) %>%
#   ggplot(aes(x = value, y = y_pred)) +
#   geom_point() +
#   geom_point(aes(y = NDSI), colour = 'red', alpha = 0.2) +
#   facet_wrap(~covariate, scales = 'free_x')
# 
# temp_df %>%
#   mutate(y_diff = NDSI - y_pred) %>%
#   gather(covariate, value, -NDSI, -y_pred, -y_diff) %>%
#   ggplot(aes(x = value, y = y_diff)) +
#   geom_point() +
#   facet_wrap(~covariate, scales = 'free_x')
# 
# ###########################################
# # ACI
# i = responses[1]
# temp_y = data.frame(site = indices_df$site, ACI = indices_df[[i]])
# temp_df = merge(x = temp_y, y = mod_df) %>% 
#   select(-c(site,wavs))
# 
# ggplot(temp_df, aes(x = 1/(ACI)^2)) + geom_histogram()
# ggplot(temp_df, aes(x = log(ACI))) + geom_histogram()
# 
# # visualize pairs of corr
# df %>%
#   GGally::ggpairs(aes(alpha = 0.05))
# 
# temp_df %>% 
#   mutate(y = log(y)) %>%
#   response_exam(i)
# 
# # FIT GAM
# fit1 = gam(ACI ~ s(Anthropophony_mean) + s(Anthropophony_mean) + s(Biophony_mean) + s(Geophony_mean) + 
#              s(Quiet_mean) + s(Interference_mean) + s(Unidentified_mean), 
#            data = temp_df)
# summary(fit1)
# 
# plot.gam(fit1, residuals = TRUE, # residuals
#          shade = TRUE, # 95% CI for the mean shape of effects
#          shade.col = "lightblue", seWithMean = TRUE, pages = 1, pch =1, cex = 1, all.terms = TRUE)
# 
# plot(fit1, pages = 1, all.terms = TRUE) # partial effects plots
# 
# # small p-values indicate residuals are not randomly distributed (i.e., not enough basis fxs)
# gam.check(fit1)
# # if p is significant, try increasing k: gam(y ~ s(x1, k = 12))
# 
# # FIT 2
# # geophony is weakly significant, increase k
# fit2 = gam(log(ACI) ~ s(Anthropophony_mean) + s(Anthropophony_mean) + s(Biophony_mean) + s(Geophony_mean, k = 12) + 
#              s(Quiet_mean) + s(Interference_mean) + s(Unidentified_mean), 
#            data = temp_df)
# summary(fit2)
# gam.check(fit2)
# 
# # values over ~0.8 for worst may mean issues with collinearity
# concurvity(fit2, full = TRUE)
# concurvity(fit2, full = FALSE) # use to id where the relationships are UNID and Int is highest
# 
# y_pred = predict(fit2, temp_df[,2:7])
# 
# 
# # FIT 3 TRYING WITH TRANSFORMED ACI
# fit3 = gam(log(ACI) ~ s(log(Anthropophony_mean)) + s(log(Anthropophony_mean)) + s(log(Biophony_mean)) + 
#              s(log(Geophony_mean)) + s(log(Quiet_mean)) + s(log(Interference_mean)) + s(log(Unidentified_mean)), 
#            data = temp_df)
# summary(fit3)
# gam.check(fit3)
# plot(fit3, pages = 1, all.terms = TRUE) # partial effects plots
# autoplot(fit3)
