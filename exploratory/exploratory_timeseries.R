library(tidyverse)
library(fpp3)
library(forecast)

# Load data
phenoDat <- read_csv("phenology-targets.csv.gz") %>% 
  filter(siteID %in% c("CLBJ", "HARV")) %>% 
  as_tsibble(index = time, key = siteID)


# Linear interpolation code
min_time <- min(phenoDat$time[!is.na(phenoDat$gcc_90)])
max_time <- max(phenoDat$time[!is.na(phenoDat$gcc_90)])
phenoDat <- phenoDat %>% filter(time >= min_time,
                                time <= max_time)

all_dates <- seq(min(phenoDat$time), max(phenoDat$time), by = 1)
na_inds <- which(is.na(phenoDat$gcc_90))

i <- 0
while (i < nrow(phenoDat)) {
  i <- i + 1
  
  if (is.na(phenoDat$gcc_90[i])) {
    start_index <- i - 1
    while (is.na(phenoDat$gcc_90[i]) &&
           i < nrow(phenoDat)) {
      i <- i + 1
    }
    end_index <- i
    
    phenoDat$gcc_90[(start_index + 1):(end_index - 1)] <-
      approx(x = phenoDat$time[c(start_index, end_index)],
             y = phenoDat$gcc_90[c(start_index, end_index)], 
             xout = phenoDat$time[(start_index + 1):(end_index - 1)])$y
  }
}


# Play with timeseries tools
autoplot(phenoDat, gcc_90)


# Default STL
dcmp <- phenoDat %>%
  model(STL(gcc_90))


autoplot(components(dcmp))


# Only year trend
dcmp <- phenoDat %>%
  model(STL(gcc_90 ~ season("year")))

autoplot(components(dcmp))

library(lubridate)
# Lowercase version

phenoDat_filtered <- phenoDat %>% 
  filter(!(month(time) == 2 & mday(time) == 29))
  

pheno_ts <- ts(data = phenoDat_filtered$gcc_90, 
     start = c(year(phenoDat$time[1]), yday(phenoDat$time[1])),
     end = c(year(phenoDat$time[nrow(phenoDat)]), yday(phenoDat$time[nrow(phenoDat)])), 
     frequency = 365, 
     names = phenoDat_filtered$siteID)

fit <- mstl(pheno_ts, s.window = 365, t.window = 13)
summary(fit)
autoplot(fit)


nday <- 366
fc <- stlf(pheno_ts, h = nday, s.window = 2) %>% 
  summary(print = F) %>% 
  rename(gcc_90 = `Point Forecast`, lower = `Lo 95`, upper = `Hi 95`)
fc$time <- phenoDat$time[nrow(phenoDat)] + 1:nday
fc$type <- "Forecast"

phenoDat$type <- "Historic"
phenoDat$lower <- NA
phenoDat$upper <- NA

ggplot(mapping = aes(time, gcc_90, col = type)) +
  geom_line(data = phenoDat) +
  geom_line(data = fc) +
  geom_errorbar(aes(ymin = lower, ymax = upper), data = fc, alpha = 0.1)


# plot(fc$time, fc$upper - fc$lower)




