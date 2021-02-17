# Redo the analysis using GAMMs w/ the package mgcv (instead of STL)
# Advatnages:
#   - No need to interpolate NAs
#   - Clear definition of uncertainty
#   - Easy to model residual variation against covariates

# Current strategy: use nested smooths to estimate the cyclic annual pattern + trend,
# then predict from only cyclic

library(tidyverse)
library(mgcv)

# Load data
phenoDat <- read_csv("phenology-targets.csv.gz") %>% 
  mutate(type = "Historic") %>% 
  filter(!is.na(gcc_90))

newdat <- data.frame(time = max(phenoDat$time) + 1:100) %>% 
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
    filter(!is.na(gcc_90))


  mod <- gamm(gcc_90 ~ s(yday, bs = "cc", k = 150) + s(time, bs = "cr", k = 10),
              data = thisDat, method = "REML",
              correlation = corAR1(form = ~ 1 | year))
  
  predicted <- predict.gam(mod$gam, newdat,
                           se.fit = T, type = "iterms", terms = "s(yday)")

  predict_list[[i]] <- data.frame(
    gcc_90 = predicted$fit + attr(predicted, "constant"),
    gcc_90_SE = predicted$se.fit,
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

