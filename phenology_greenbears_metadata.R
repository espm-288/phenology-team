## April 12, working on metadata script before adding to automated script. 

devtools::install_github("eco4cast/neon4cast")

library(neon4cast)
library(neonstore)


write_meta_template("submissions/phenology-", Sys.Date(), "-greenbears_par.csv")

generate_metadata(forecast_file, 
                  metadata_yaml, 
                  forecast_issue_time, 
                  forecast_iteration_id)