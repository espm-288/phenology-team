library(tidyverse)
library(mgcv)
library(neonstore)

# Load data
phenoDat <- read_csv("phenology-targets.csv.gz") %>% 
  mutate(type = "Historic") %>% 
  filter(!is.na(gcc_90))

site_names <- c("HARV", "BART", "SCBI", "STEI", "UKFS", "GRSM", "DELA", "CLBJ")


# PAR photosynthetically active radiation
neon_download("DP1.00024.001", site = site_names, file_regex = ".*30min.*\\.csv",
              start_date = min(phenoDat$time), end_date = max(phenoDat$time),
              .token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJiZW4uZ29sZHN0ZWluQGJlcmtlbGV5LmVkdSIsInNjb3BlIjoicmF0ZTpwdWJsaWMiLCJpc3MiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnLyIsImV4cCI6MTc0ODkwNDQ2MCwiaWF0IjoxNTkxMjI0NDYwLCJlbWFpbCI6ImJlbi5nb2xkc3RlaW5AYmVya2VsZXkuZWR1In0.dDKtuk-ZZriLnNOkvKXG-IowZii7uhWNRr13xcw5FwXI1k0-4tQSW3oxKjPbfJF6sG9fRokJbJFqhVZRTgj_KA")


# Relative humidity
neon_download("DP1.00098.001", site = site_names, file_regex = ".*30min.*\\.csv",
              start_date = min(phenoDat$time), end_date = max(phenoDat$time),
              .token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJiZW4uZ29sZHN0ZWluQGJlcmtlbGV5LmVkdSIsInNjb3BlIjoicmF0ZTpwdWJsaWMiLCJpc3MiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnLyIsImV4cCI6MTc0ODkwNDQ2MCwiaWF0IjoxNTkxMjI0NDYwLCJlbWFpbCI6ImJlbi5nb2xkc3RlaW5AYmVya2VsZXkuZWR1In0.dDKtuk-ZZriLnNOkvKXG-IowZii7uhWNRr13xcw5FwXI1k0-4tQSW3oxKjPbfJF6sG9fRokJbJFqhVZRTgj_KA")

# Temperature
neon_download("DP1.00005.001", site = site_names, file_regex = ".*30min.*\\.csv",
              start_date = min(phenoDat$time), end_date = max(phenoDat$time),
              .token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJiZW4uZ29sZHN0ZWluQGJlcmtlbGV5LmVkdSIsInNjb3BlIjoicmF0ZTpwdWJsaWMiLCJpc3MiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnLyIsImV4cCI6MTc0ODkwNDQ2MCwiaWF0IjoxNTkxMjI0NDYwLCJlbWFpbCI6ImJlbi5nb2xkc3RlaW5AYmVya2VsZXkuZWR1In0.dDKtuk-ZZriLnNOkvKXG-IowZii7uhWNRr13xcw5FwXI1k0-4tQSW3oxKjPbfJF6sG9fRokJbJFqhVZRTgj_KA")

# Plant phenophase
neon_download("DP1.10055.001", site = site_names, 
              start_date = min(phenoDat$time), end_date = max(phenoDat$time),
              .token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJiZW4uZ29sZHN0ZWluQGJlcmtlbGV5LmVkdSIsInNjb3BlIjoicmF0ZTpwdWJsaWMiLCJpc3MiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnLyIsImV4cCI6MTc0ODkwNDQ2MCwiaWF0IjoxNTkxMjI0NDYwLCJlbWFpbCI6ImJlbi5nb2xkc3RlaW5AYmVya2VsZXkuZWR1In0.dDKtuk-ZZriLnNOkvKXG-IowZii7uhWNRr13xcw5FwXI1k0-4tQSW3oxKjPbfJF6sG9fRokJbJFqhVZRTgj_KA")




neonstore::neon_store(product = "DP1.00024.001")
neonstore::neon_store(product = "DP1.00098.001")
neonstore::neon_store(product = "DP1.00005.001")
neonstore::neon_store(product = "DP1.10055.001")


