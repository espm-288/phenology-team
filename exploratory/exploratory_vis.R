


joined <- left_join(allData, temptbl, by = c("time" = "date", "siteID")) %>% 
  left_join(humidtbl, by = c("time" = "date", "siteID")) %>% 
  left_join(radtbl, by = c("time" = "date", "siteID")) %>% 
  filter(!is.na(wssTempTripleMean))
# 
# joined <- joined %>% 
#   mutate()





joined %>% 
  mutate(month = lubridate::month(time), year = lubridate::year(time)) %>% 
  group_by(month, year, siteID) %>% 
  summarize(temp = mean(wssTempTripleMaximum), gcc = mean(gcc_90)) %>% 
  mutate(time = lubridate::date(paste(year, month, 1, sep = "-"))) %>% 
  ggplot(aes(time, temp, col = siteID)) +
  geom_point() + geom_line()


joined %>% 
  ggplot(aes(time, scale(wssTempTripleMean), col = siteID)) +
  geom_line() +
  geom_point(aes(time, scale(gcc_90), col = siteID)) +
  facet_wrap(~siteID)



