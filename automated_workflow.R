# automated_workflow.R
# author: greenbears team
# This script automates our whole pipeline for producing forecasts of NEON
# phenology data and then publishing them to the NEON challenge page.

library(tidyverse)
library(mgcv)
library(neonstore)
library(dplyr)
library(aws.s3)

#### Step 1. Get today's data ####

source("downloadPhenoCam.R")
source("calculatePhenoCamUncertainty.R")

# Selected Sites for Challenge
siteIDs <- c("NEON.D01.HARV.DP1.00033","NEON.D01.BART.DP1.00033","NEON.D02.SCBI.DP1.00033",
             "NEON.D05.STEI.DP1.00033","NEON.D06.UKFS.DP1.00033","NEON.D07.GRSM.DP1.00033",
             "NEON.D08.DELA.DP1.00033","NEON.D11.CLBJ.DP1.00033")

site_names <- c("HARV", "BART", "SCBI", "STEI", "UKFS", "GRSM", "DELA", "CLBJ")

allData <- data.frame(matrix(nrow = 0, ncol = 5))

message(paste0("Downloading and generating phenology targets ", Sys.time()))

for(i in 1:length(siteIDs)){
  siteName <- siteIDs[i]
  message(siteName)
  if(siteName != "NEON.D11.CLBJ.DP1.00033"){
    URL_gcc90 <- paste('https://phenocam.sr.unh.edu/data/archive/',siteName,"/ROI/",siteName,"_DB_1000_1day.csv",sep="") ##URL for daily summary statistics
    URL_temp <- paste('https://phenocam.sr.unh.edu/data/archive/',siteName,"/ROI/",siteName,"DP4.00001.001.csv",sep="") ##URL for daily summary statistics
    URL_individual <- paste('https://phenocam.sr.unh.edu/data/archive/',siteName,"/ROI/",siteName,"_DB_1000_roistats.csv",sep="") ##URL for individual image metrics
  }else{
    URL_gcc90 <- paste('https://phenocam.sr.unh.edu/data/archive/',siteName,"/ROI/",siteName,"_DB_2000_1day.csv",sep="") ##URL for daily summary statistics
    URL_individual <- paste('https://phenocam.sr.unh.edu/data/archive/',siteName,"/ROI/",siteName,"_DB_2000_roistats.csv",sep="") ##URL for individual image metrics
  }
  phenoData <- download.phenocam(URL = URL_gcc90)
  dates <- unique(phenoData$date)
  phenoData_individual <- download.phenocam(URL=URL_individual,skipNum = 17)
  gcc_sd <- calculate.phenocam.uncertainty(dat=phenoData_individual,dates=dates) ##Calculates standard deviations on daily gcc90 values
  
  subPhenoData <- phenoData %>% 
    mutate(siteID = stringr::str_sub(siteName, 10, 13), 
           time = date) %>% 
    select(time, siteID, gcc_90)
  subPhenoData <- cbind(subPhenoData,gcc_sd)
  
  allData <- rbind(allData,subPhenoData)
  
}
readr::write_csv(allData, "phenology-targets.csv.gz")



#### Step 2. Fit a model and predict for all sites ####
# The code in this section is copied from predict_w_GAMs.R

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
  
  predicted <- predict.gam(mod$gam, newdat,
                           se.fit = T, type = "iterms", terms = "s(yday)")

  predict_list[[i]] <- data.frame(
    gcc_90 = unlist(predicted$fit + attr(predicted, "constant")),
    gcc_90_SE = predicted$se.fit,
    time = newdat$date,
    siteID = sites[[i]]
  )
}

#### Step 3. Format output ####

predict_df <- do.call(rbind, predict_list)
colnames(predict_df) <- c("gcc_90", "gcc_90_SE", "time", "siteID")

final_predict_df <- predict_df %>% 
  select(time, siteID, gcc_90, gcc_90_SE) %>% 
  rename(mean = gcc_90, sd = gcc_90_SE) %>% 
  pivot_longer(cols = c("mean", "sd"), 
               names_to = "statistic",
               values_to = "gcc_90") %>% 
  mutate(forecast = 1, data_assimilation = 0)

#### Step 4. Publish ####
name <- paste0("submissions/phenology-", Sys.Date(), "-greenbears_gams.csv")

write_csv(final_predict_df, name)

if (aws.s3::put_object(file = name, bucket = "submissions", region="data", base_url = "ecoforecast.org")) {
  cat("SUBMITTED!!!!!!!")
}




