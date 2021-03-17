



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




neonstore::neon_download()
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
