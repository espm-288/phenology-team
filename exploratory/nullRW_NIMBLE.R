library("nimble")
library("ecoforecastR")
library("ncdf4")
library("tidyverse")
library("tidybayes")

source("randomWalkNullModelFunction.R")
###Random WalkNull Model Calculations
####Note: Currently this is not set up to run iteratively because I am not sure how the challenge is planning on doing this.
####Note (continued): Hopefully someone who knows about how this will be done can use this code to do that

generate_plots <- TRUE
team_name <- "EFInull"

download.file("https://data.ecoforecast.org/targets/phenology/phenology-targets.csv.gz",
              "phenology-targets.csv.gz")

phenoDat <- read.csv("phenology-targets.csv.gz",header=TRUE)
sites <- unique(as.character(phenoDat$siteID))

forecast_length <- 35
predictions <- array(NA, dim = c(forecast_length, length(sites), 2000))
parameters <- array(NA, dim = c(length(sites),1000))


#'Generic random walk state-space model is JAGS format.  We use this model for
#'both the oxygen and temperature null forecasts
RWcode <- nimbleCode({
  # Priors
  x[1] ~ dnorm(x_ic,tau_add)
  tau_obs[1] <- 1 / pow(sd_obs[1], 2)
  y[1] ~ dnorm(x[1],tau_obs[1])
  sd_add  ~ dunif(0.000001, 100)
  tau_add <- 1/ pow(sd_add, 2)
  # Process Model
  for(t in 2:N){
    x[t] ~ dnorm(x[t-1], tau_add)
    tau_obs[t] <- 1 / pow(sd_obs[t], 2)
    y[t] ~ dnorm(x[t], tau_obs[t])
  }
})

forecast_saved <- NULL

for(s in 1:length(sites)){
  
  message(paste0("forecasting site: ",sites[s]))
  
  forecast_length <- 35
  
  
  
  sitePhenoDat <- phenoDat[phenoDat$siteID==sites[s],]
  sitePhenoDat$time <- lubridate::as_date(sitePhenoDat$time)
  
  start_forecast <- max(sitePhenoDat$time) + lubridate::days(1)
  
  sitePhenoDat <- sitePhenoDat
  full_time <- tibble::tibble(time = seq(min(sitePhenoDat$time), max(sitePhenoDat$time) + lubridate::days(forecast_length), by = "1 day"))
  forecast_start_index <- which(full_time$time == max(sitePhenoDat$time) + lubridate::days(1))
  d <- tibble::tibble(time = sitePhenoDat$time,
                      p=as.numeric(sitePhenoDat$gcc_90),
                      p.sd=as.numeric(sitePhenoDat$gcc_sd))
  d <- dplyr::full_join(d, full_time)
  
  ggplot(d, aes(x = time, y = p)) +
    geom_point()
  
  
  #gap fill the missing precisions by assigning them the average sd for the site
  d$p.sd[!is.finite(d$p.sd)] <- NA
  d$p.sd[is.na(d$p.sd)] <- mean(d$p.sd,na.rm=TRUE)
  d$p.sd[d$p.sd == 0.0] <- min(d$p.sd[d$p.sd != 0.0])
  
  # Linear interpolation of missing ps
  init_x <- approx(x = d$time[!is.na(d$p)], y = d$p[!is.na(d$p)], xout = d$time, rule = 2)$y
  
  
  # Create model
  RWmodel <-nimbleModel(
    RWcode,
    constants = list(
      N = length(d$p)
    ),
    data = list(
      y = d$p,
      sd_obs = d$p.sd,
      x_ic = 0.3
    ), inits = list(
      sd_add = sd(diff(d$p[!is.na(d$p)])),
      x = init_x
    )
  )
  RWmcmcConf <- configureMCMC(RWmodel)
  RWmcmcConf$addMonitors(c("x", "y"))
  RWmcmc <- buildMCMC(RWmcmcConf)
  
  complist <- compileNimble(RWmcmc, RWmodel)
  
  m <- runMCMC(complist$RWmcmc, niter = 20000, nburnin = 10000, thin = 5)
  
  #Use TidyBayes package to clean up the output
  model_output <- m %>%
    spread_draws(y[day]) %>%
    filter(.chain == 1) %>%
    rename(ensemble = .iteration) %>%
    mutate(time = full_time$time[day]) %>%
    ungroup() %>%
    select(time, y, ensemble)
  
  
  
  if(generate_plots){
    #Pull in the observed data for plotting
    obs <- tibble(time = d$time,
                  obs = d$p)
    
    
    #Post past and future
    model_output %>%
      group_by(time) %>%
      summarise(mean = mean(y),
                upper = quantile(y, 0.975),
                lower = quantile(y, 0.025),.groups = "drop") %>%
      ggplot(aes(x = time, y = mean)) +
      geom_line() +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = "lightblue", fill = "lightblue") +
      geom_point(data = obs, aes(x = time, y = obs), color = "red") +
      labs(x = "Date", y = "oxygen")
    
    ggsave(paste0("phenology_",sites[s],"_figure.pdf"), device = "pdf")
  }
  
  #Filter only the forecasted dates and add columns for required variable
  forecast_saved_tmp <- model_output %>%
    filter(time >= start_forecast) %>%
    rename(gcc_90 = y) %>%
    mutate(data_assimilation = 0,
           forecast = 1,
           obs_flag = 2,
           siteID = sites[s]) %>%
    mutate(forecast_iteration_id = start_forecast) %>%
    mutate(forecast_project_id = team_name)
  
  
  predictions[ ,s , ] <- forecast_saved_tmp %>%
    pivot_wider(names_from = ensemble, values_from = gcc_90) %>%
    select(-c("data_assimilation","forecast", "obs_flag", "siteID", "forecast_iteration_id", "forecast_project_id","time")) %>%
    as.matrix()
  
  
  # Combined with the previous sites
  forecast_saved <- rbind(forecast_saved, forecast_saved_tmp)
  
}

forecast_time <- unique(forecast_saved$time)

##Put in EFI standard form (based on EFI standards logistic-metadata-example vignette)
##Forecast Identifiers (please change to what is needed)

forecast_project_id <- team_name
forecast_model_id <- "v0.1"
forecast_iteration_id <- Sys.time()

#plot(predictions[1,], type = 'l', ylim = range(c(predictions)))
#for(i in 2:nrow(predictions)){
#  points(predictions[i,], type = "l")
#}