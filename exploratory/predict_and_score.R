# How many days should we back-predict?
n_predict <- 35

library(tidyverse)
library(mgcv)
library(neonstore)
source("score_it.R")

# Load data
phenoDat <- read_csv("phenology-targets.csv.gz") %>% 
  mutate(type = "Historic") %>% 
  filter(!is.na(gcc_90))

# Neon data
site_names <- c("HARV", "BART", "SCBI", "STEI", "UKFS", "GRSM", "DELA", "CLBJ")
par_good <- read_csv("data/PAR.csv")


phenoDat <- left_join(phenoDat, par_good, by = c("time" = "date",
                                                 "siteID"))

phenoDat_future <- phenoDat %>% 
  filter(time %in% (max(phenoDat$time) - n_predict + 1):max(phenoDat$time))

phenoDat <- phenoDat %>% 
  filter(time <= max(phenoDat$time) - n_predict)
newdat <- data.frame(time = max(phenoDat$time) + 1:n_predict) %>% 
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
  
  predict_list[[i]] <- data.frame(
    gcc_90 = unlist(predicted$fit + attr(predicted, "constant")),
    gcc_90_SE = predicted$se.fit,
    time = newdat$date,
    siteID = sites[[i]]
  )
}

predict_df <- do.call(rbind, predict_list)
colnames(predict_df) <- c("gcc_90", "gcc_90_SE", 
                          "time", "siteID")



predict_df %>% 
  select(time, siteID, gcc_90, gcc_90_SE) %>% 
  rename(mean = gcc_90, sd = gcc_90_SE) %>% 
  pivot_longer(cols = c("mean", "sd"), 
               names_to = "statistic",
               values_to = "gcc_90") %>% 
  write_csv("phenology-2021-2-22-greenbearsPRACTICE.csv")


score_it(targets_file = "phenology-targets.csv.gz", 
         forecast_files = "phenology-2021-03-10-greenbears_gams.csv", 
         target_variables = "gcc_90")

score <- read_csv("scores/scores-phenology-2021-03-10-greenbears_gams.csv.gz")

uncertainty_df <- predict_df %>% 
  filter(statistic == "sd") %>% 
  select(-statistic) %>% 
  rename(predict_sd = gcc_90)

comparison_df <- predict_df %>% 
  filter(statistic == "mean") %>% 
  rename(prediction = gcc_90) %>% 
  left_join(phenoDat_future, by = c("siteID", "time")) %>% 
  left_join(uncertainty_df, by = c("siteID", "time"))

comparison_df %>% 
  ggplot() +
  geom_ribbon(mapping = aes(time, ymin = prediction - 1.96*predict_sd, 
                            ymax = prediction + 1.96*predict_sd),
              alpha = 0.2) +
  geom_line(mapping = aes(time, prediction)) +
  geom_point(mapping = aes(time, gcc_90)) +
  facet_wrap(~siteID)

site_scores <- score %>% 
  group_by(siteID) %>% 
  summarize(score = mean(score, na.rm = T))

