library(tidyverse)
library(fpp3)


phenoDat <- read_csv("phenology-targets.csv.gz") %>% 
  filter(siteID == "CLBJ") %>% 
  as_tsibble(index = time, key = siteID)

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


autoplot(phenoDat, gcc_90)


# Default STL
dcmp <- phenoDat %>%
  model(STL(gcc_90 ~ season("year")))


autoplot(components(dcmp))


# Only year trend
dcmp <- phenoDat %>%
  model(STL(gcc_90 ~ season("year")))


autoplot(components(dcmp))

