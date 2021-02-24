# How many days should we back-predict?
n_predict <- 35

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
