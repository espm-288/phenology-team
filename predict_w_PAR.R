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

par_means <- par_good %>% 
  mutate(yday = lubridate::yday(date)) %>% 
  group_by(yday, siteID) %>% 
  summarize(par = mean(par, na.rm = T))

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
  
  this_newdat <- newdat %>% mutate(siteID = sites[[i]])

  PAR_gam <- gamm(par ~ s(yday, bs = "cc", k = 150) + s(time, bs = "cr", k = 10),
                  data = thisDat, method = "REML",
                  correlation = corAR1(form = ~ 1 | year))
  predicted_par <- predict.gam(PAR_gam$gam, newdata = this_newdat, 
                                 type = "iterms", terms = "s(yday)")
  this_newdat$par <- as.numeric(attr(predicted_par, "constant") + predicted_par)
  
  
  mod_wpar <- gamm(gcc_90 ~ s(yday, bs = "cc", k = 150) + s(time, bs = "cr", k = 10) + par,
              data = thisDat, method = "REML",
              correlation = corAR1(form = ~ 1 | year))
  
  predicted_wpar <- 
    predict.gam(mod_wpar$gam, this_newdat,
                se.fit = T, type = "iterms", terms = c("s(yday)", "par"))

  predict_list[[i]] <- data.frame(
    gcc_90 = predicted_wpar$fit[,1] + predicted_wpar$fit[,2] + 
                  attr(predicted_wpar, "constant"),
    gcc_90_SE = sqrt(predicted_wpar$se.fit[,1]^2 + predicted_wpar$se.fit[,2]^2),
    time = newdat$date,
    siteID = sites[[i]]
  )
}

predict_df <- do.call(rbind, predict_list)
colnames(predict_df) <- c("gcc_90", "gcc_90_SE", "time", "siteID")

# Visualize predictions
predict_plots <- list()
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
}




predict_df %>% 
  select(time, siteID, gcc_90, gcc_90_SE) %>% 
  rename(mean = gcc_90, sd = gcc_90_SE) %>% 
  pivot_longer(cols = c("mean", "sd"), 
               names_to = "statistic",
               values_to = "gcc_90") %>% 
  mutate(forecast = 1, data_assimilation = 0) %>% 
  write_csv(paste0("phenology-2021-03-29-greenbears_PAR.csv"))

  

