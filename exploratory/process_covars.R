library(tidyverse)


rh_table <- neonstore::neon_table("RH_30min-expanded")
rh_good <- rh_table %>% 
  select(startDateTime, RHMean, tempRHMean, siteID) %>% 
  mutate(date = lubridate::date(startDateTime)) %>% 
  group_by(date, siteID) %>% 
  summarize(rh = mean(RHMean), temp_mean = mean(tempRHMean))
write_csv(rh_good, "data/relative_humidity.csv")

rh_good %>% 
  filter(siteID == "BART", lubridate::year(date) == 2019) %>% 
  ggplot(aes(date, temp_mean)) + geom_point()




neon_download("DP1.00024.001", site = site_names, file_regex = ".*30min.*\\.csv",
              start_date = min(phenoDat$time), end_date = max(phenoDat$time),
              .token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJiZW4uZ29sZHN0ZWluQGJlcmtlbGV5LmVkdSIsInNjb3BlIjoicmF0ZTpwdWJsaWMiLCJpc3MiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnLyIsImV4cCI6MTc0ODkwNDQ2MCwiaWF0IjoxNTkxMjI0NDYwLCJlbWFpbCI6ImJlbi5nb2xkc3RlaW5AYmVya2VsZXkuZWR1In0.dDKtuk-ZZriLnNOkvKXG-IowZii7uhWNRr13xcw5FwXI1k0-4tQSW3oxKjPbfJF6sG9fRokJbJFqhVZRTgj_KA")
neonstore::neon_store(product = "DP1.00024.001")
par_table <- neonstore::neon_table("PARPAR_30min-expanded")
par_good <- par_table %>% 
  select(startDateTime, PARMean, siteID) %>% 
  mutate(date = lubridate::date(startDateTime)) %>% 
  group_by(date, siteID) %>% 
  summarize(par = mean(PARMean))
write_csv(par_good, "data/PAR.csv")

par_good %>% 
  filter(siteID == "BART") %>% 
  ggplot(aes(date, par)) + geom_point()




phe_table <- neonstore::neon_table("phe_perindividual-basic")
par_good <- par_table %>% 
  select(startDateTime, PARMean, siteID) %>% 
  mutate(date = lubridate::date(startDateTime)) %>% 
  group_by(date, siteID) %>% 
  summarize(par = mean(PARMean))
write_csv(par_good, "data/PAR.csv")

par_good %>% 
  filter(siteID == "BART") %>% 
  ggplot(aes(date, par)) + geom_point()
