library(tidyverse)
library(fpp3)
library(forecast)
library(lubridate)



# Load data
phenoDat <- read_csv("phenology-targets.csv.gz")
sites <- unique(phenoDat$siteID)
plots <- list()
resultlist <- list()

# For-loop to forecast by siteID
for(i in 1:length(sites)){
  thisDat <- filter(phenoDat, siteID == sites[i],
    !(month(time) == 2 & mday(time) == 29))
  
  min_time <- min(thisDat$time[!is.na(thisDat$gcc_90)])
  max_time <- max(thisDat$time[!is.na(thisDat$gcc_90)])
  thatDat <- thisDat %>% filter(time >= min_time,
                                  time <= max_time)
  #linear interpolation for na values
  j <- 0
  while (j < nrow(thatDat)) {
    j <- j + 1
    
    if (is.na(thatDat$gcc_90[j])) {
      start_index <- j - 1
      while (is.na(thatDat$gcc_90[j]) &&
             j < nrow(thatDat)) {
        j <- j + 1
      }
      end_index <- j
      
      thatDat$gcc_90[(start_index + 1):(end_index - 1)] <-
        approx(x = thatDat$time[c(start_index, end_index)],
               y = thatDat$gcc_90[c(start_index, end_index)], 
               xout = thatDat$time[(start_index + 1):(end_index - 1)])$y
    }
  }
  
  pheno_ts <- ts(data = thatDat$gcc_90, 
                 start = c(year(thatDat$time[1]), yday(thatDat$time[1])),
                 end = c(year(thatDat$time[nrow(thatDat)]), yday(thatDat$time[nrow(thatDat)])), 
                 frequency = 365)
  
  
  nday <- 35
  result <- capture.output(fc <- stlf(pheno_ts, h = nday, s.window = 2) %>% 
    summary(print = FALSE) %>% 
    rename(gcc_90 = `Point Forecast`, lower = `Lo 95`, upper = `Hi 95`))
  
    
  fc$time <- thatDat$time[nrow(thatDat)] + 1:nday
  fc$type <- "Forecast"
  fc$siteID <- sites[i]
  resultlist[[i]] <- fc
  
  thatDat$type <- "Historic"
  thatDat$lower <- NA
  thatDat$upper <- NA
  
  plots[[i]] <- ggplot(mapping = aes(time, gcc_90, col = type)) +
    geom_line(data = thatDat) +
    geom_line(data = fc) +
    geom_errorbar(aes(ymin = lower, ymax = upper), data = fc, alpha = 0.1) +
    ggtitle(sites[i])
  
  

}

plots

forcastdata <- do.call(rbind, resultlist) %>% 
  mutate(sd = (gcc_90 - lower)/1.96)

forcastdata %>% 
  select(time, siteID, gcc_90, sd) %>% 
  rename(mean = gcc_90) %>% 
  pivot_longer(cols = c("mean", "sd"), 
               names_to = "statistic",
               values_to = "gcc_90") %>% 
  write_csv(paste0("phenology-", ("2021-03-08"), "-greenbears_stl.csv"))
