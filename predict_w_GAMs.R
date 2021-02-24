# Redo the analysis using GAMMs w/ the package mgcv (instead of STL)
# Advatnages:
#   - No need to interpolate NAs
#   - Clear definition of uncertainty
#   - Easy to model residual variation against covariates

# Current strategy: use nested smooths to estimate the cyclic annual pattern + trend,
# then predict from only cyclic

# Other variables to incorporate:
#   Weather? (Precip, temp, solar radiation) > forecastable
#   Lagged covariates?

library(tidyverse)
library(mgcv)
library(neonstore)

# Load data
phenoDat <- read_csv("phenology-targets.csv.gz") %>% 
  mutate(type = "Historic") %>% 
  filter(!is.na(gcc_90))

# Neon data
site_names <- c("HARV", "BART", "SCBI", "STEI", "UKFS", "GRSM", "DELA", "CLBJ")
par_good <- read_csv("data/PAR.csv")


phenoDat <- left_join(phenoDat, par_good, by = c("time" = "date",
                                                 "siteID"))


newdat <- data.frame(time = max(phenoDat$time) + 1:35) %>% 
  mutate(yday = lubridate::yday(time), year = lubridate::year(time),
         date = time,
         time = as.numeric(time))

sites <- unique(phenoDat$siteID)
predict_list <- list()

# Loop over sites
for (i in 1:length(sites)) {
  cat("Predicting for site", sites[[i]], "\n")
  
  thisDat <- phenoDat %>% 
    filter(siteID == sites[[i]]) %>% 
    # as_tsibble(index = time, key = siteID) %>% 
    mutate(yday = lubridate::yday(time), year = lubridate::year(time),
           time = as.numeric(time)) %>% 
    filter(!is.na(gcc_90), !is.na(par))
  newdat$par <- mean(thisDat$par)


  mod <- gamm(gcc_90 ~ s(yday, bs = "cc", k = 150) + s(time, bs = "cr", k = 10),
              data = thisDat, method = "REML",
              correlation = corAR1(form = ~ 1 | year))
  
  mod_wpar <- gamm(gcc_90 ~ s(yday, bs = "cc", k = 150) + s(time, bs = "cr", k = 10) + par,
              data = thisDat, method = "REML",
              correlation = corAR1(form = ~ 1 | year))
  
  predicted <- predict.gam(mod$gam, newdat,
                           se.fit = T, type = "iterms", terms = "s(yday)")
  predicted_wpar <- 
    predict.gam(mod_wpar$gam, newdat,
                se.fit = T, type = "iterms", terms = c("s(yday)", "par"))

  predict_list[[i]] <- data.frame(
    gcc_90 = unlist(predicted$fit + attr(predicted, "constant")),
    gcc_90_SE = predicted$se.fit,
    gcc_90_wpar = predicted_wpar$fit[,1] + predicted_wpar$fit[,2] + 
                  attr(predicted_wpar, "constant"),
    gcc_90_SE_wpar = sqrt(predicted_wpar$se.fit[,1]^2 + predicted_wpar$se.fit[,2]^2),
    time = newdat$date,
    siteID = sites[[i]]
  )
}

predict_df <- do.call(rbind, predict_list)
colnames(predict_df) <- c("gcc_90", "gcc_90_SE", 
                          "gcc_90_wpar", "gcc_90_SE_wpar",
                          "time", "siteID")

# Visualize predictions
predict_plots <- list()
predict_plots_wpar <- list()
for (i in 1:length(sites)) {
  predict_plots[[i]] <- predict_df %>% 
    mutate(type = "Predict") %>% 
    bind_rows(phenoDat[, c("gcc_90", "type", "time", "siteID")]) %>% 
    filter(siteID == sites[i]) %>% 
    ggplot(aes(time, gcc_90, col = type)) +
    geom_point() +
    geom_errorbar(aes(time, ymin = gcc_90 - 1.96*gcc_90_SE, 
                      ymax = gcc_90 + 1.96*gcc_90_SE, col = type, alpha = 0.1)) +
    scale_alpha_continuous(guide = FALSE) +
    ggtitle(sites[[i]])
  
  predict_plots_wpar[[i]] <- predict_df %>% 
    mutate(type = "Predict", gcc_90 = gcc_90_wpar, gcc_90_SE = gcc_90_SE_wpar) %>% 
    bind_rows(phenoDat[, c("gcc_90", "type", "time", "siteID")]) %>% 
    filter(siteID == sites[i]) %>% 
    ggplot(aes(time, gcc_90, col = type)) +
    geom_point() +
    geom_errorbar(aes(time, ymin = gcc_90 - 1.96*gcc_90_SE, 
                      ymax = gcc_90 + 1.96*gcc_90_SE, col = type, alpha = 0.1)) +
    scale_alpha_continuous(guide = FALSE) +
    ggtitle(sites[[i]])
}




predict_df %>% 
  select(time, siteID, gcc_90, gcc_90_SE) %>% 
  rename(mean = gcc_90, sd = gcc_90_SE) %>% 
  pivot_longer(cols = c("mean", "sd"), 
               names_to = "statistic",
               values_to = "gcc_90") %>% 
  write_csv("prediction-2021-02-22-greenbears.csv")
  
  

