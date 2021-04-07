library(tidyverse)
library(scoringRules)



read_forecast <- function(file_in, 
                          grouping_variables = c("siteID", "time"),
                          target_variables = c("oxygen", 
                                               "temperature", 
                                               "richness",
                                               "abundance", 
                                               "nee",
                                               "le", 
                                               "vswc",
                                               "gcc_90"),
                          reps_col = "ensemble",
                          ...){
  print(file_in)
  no_forecast <- FALSE
  if(any(vapply(c("[.]csv", "[.]csv\\.gz"), grepl, logical(1), file_in))){  
    # if file is csv zip file
    out <- read_csv(file_in, guess_max = 1e6, ...) 
    
    
  } else if(grepl("[.]nc", file_in)){ #if file is nc
    
    nc <- ncdf4::nc_open(file_in)
    siteID <- ncdf4::ncvar_get(nc, "siteID")
    time <- ncdf4::ncvar_get(nc, "time")
    #tustr<-strsplit(ncdf4::ncatt_get(nc, varid = "time", "units")$value, " ")
    #time <-lubridate::as_date(time,origin=unlist(tustr)[3])
    t_string <- strsplit(ncdf4::ncatt_get(nc, varid = "time", "units")$value, " ")[[1]]
    if(t_string[1] == "days"){
      tustr<-strsplit(ncdf4::ncatt_get(nc, varid = "time", "units")$value, " ")
      time <-lubridate::as_date(time,origin=unlist(tustr)[3])
    }else{
      tustr <- lubridate::as_datetime(strsplit(ncdf4::ncatt_get(nc, varid = "time", "units")$value, " ")[[1]][3])
      time <- as.POSIXct.numeric(time, origin = tustr)
    } 
    
    targets <- names(nc$var)[which(names(nc$var) %in% target_variables)]
    combined_forecast <- NULL
    for(j in 1:length(targets)){
      forecast_targets <- ncdf4::ncvar_get(nc, targets[j])
      for(i in 1:length(siteID)){
        tmp <- forecast_targets[ ,i ,]
        d <- cbind(time, as.data.frame(tmp))
        names(d) <- c("time", seq(1,dim(tmp)[2]))
        d <- d %>%
          tidyr::pivot_longer(-time, names_to = reps_col, values_to = "value") %>%
          dplyr::mutate(siteID = siteID[i],
                        variable = targets[j])
        combined_forecast <- rbind(combined_forecast, d)
      }
    }
    ncdf4::nc_close(nc)
    combined_forecast <- combined_forecast %>%
      tidyr::pivot_wider(names_from = variable, values_from = value) %>% 
      mutate(filename = basename(file_in))
    
    out <- combined_forecast
  }else{
    out <- NA
    no_forecast <- TRUE
  }
  
  if(!is.na(no_forecast)){
    teams_tmp <- (str_split(basename(file_in), c("-"),simplify = TRUE))
    if(str_detect(teams_tmp[,1], "30min")){
      unique_dates <- sort(lubridate::as_datetime(unique(out$time)))
      dates <- lubridate::as_datetime(out$time)
    }else{
      unique_dates <- sort(lubridate::as_date(unique(out$time)))
      dates <- lubridate::as_date(out$time)
    }
    
    time_step <- unique_dates[2] - unique_dates[1]
    first_date <- unique_dates[1] - time_step
    horizon = as.numeric(dates - first_date) / as.numeric(time_step)
    team = tools::file_path_sans_ext(tools::file_path_sans_ext(last(teams_tmp[, ncol(teams_tmp)])))
    
    out <- out %>% 
      mutate(forest_start_time = first_date,
             horizon = horizon,
             team = team,
             theme = teams_tmp[,1])
  }
  out
}
## Generic scoring function.
crps_score <- function(forecast, 
                       target,
                       grouping_variables = c("siteID", "time"),
                       target_variables = c("oxygen", 
                                            "temperature", 
                                            "richness",
                                            "abundance", 
                                            "nee",
                                            "le", 
                                            "vswc",
                                            "gcc_90"),
                       reps_col = c("ensemble")){
  
  
  ## drop extraneous columns && make grouping vars into chr ids (i.e. not dates)
  
  if("ensemble" %in% colnames(forecast)){ 
    reps_col <- "ensemble"
    variables <- c(grouping_variables, target_variables, reps_col)
  }else  if("statistic" %in% colnames(forecast)){ 
    reps_col <- "statistic"
    variables <- c(grouping_variables, target_variables, reps_col) 
  }
  
  forecast <- forecast %>% dplyr::select(any_of(c(variables, "forest_start_time", "horizon", "team", "theme")))
  target <- target %>% select(any_of(variables))
  
  ## Teach crps to treat any NA observations as NA scores:
  scoring_fn_ensemble <- function(y, dat) {
    tryCatch(scoringRules::crps_sample(y, dat), error = function(e) NA_real_, finally = NA_real_)
  }
  
  scoring_fn_stat <- function(y, mean, sd) {
    tryCatch(scoringRules::crps_norm(y, mean = mean, sd = sd), error = function(e) NA_real_, finally = NA_real_)
  }
  
  ## Make tables into long format
  target_long <- target %>% 
    pivot_longer(any_of(target_variables), 
                 names_to = "target", 
                 values_to = "observed")
  forecast_long <- forecast %>% 
    pivot_longer(any_of(target_variables), 
                 names_to = "target", 
                 values_to = "predicted")
  
  if(reps_col == "ensemble"){
    
    inner_join(forecast_long, target_long, by = c(grouping_variables, "target"))  %>% 
      group_by(across(any_of(c(grouping_variables, "target", "horizon", "team", "forest_start_time", "theme")))) %>% 
      summarise(score = scoring_fn_ensemble(observed[[1]], predicted),
                .groups = "drop")
    
  } else {
    
    forecast_long %>%
      pivot_wider(names_from = statistic, values_from = predicted) %>%
      inner_join(target_long, by = c(grouping_variables, "target"))  %>% 
      group_by(across(any_of(c(grouping_variables, "target", "horizon", "team", "forest_start_time", "theme")))) %>% 
      summarise(score = scoring_fn_stat(observed[[1]], mean, sd),
                .groups = "drop")
    
  }
}


score_filenames <- function(forecast_files){
  f_name <- tools::file_path_sans_ext(paste0("scores-",
                                             basename(forecast_files)), compression = TRUE)
  file.path("scores", paste0(f_name, ".csv.gz"))
}


score_it <- function(targets_file, 
                     forecast_files, 
                     target_variables,
                     grouping_variables = c("time", "siteID"),
                     reps_col = c("ensemble"),
                     score_files = score_filenames(forecast_files)
){
  
  ## Read in data and compute scores!
  target <- read_forecast(targets_file)
  forecasts <- lapply(forecast_files, read_forecast)
  
  scores <- lapply(forecasts, 
                   crps_score, 
                   target = target,  
                   target_variables = target_variables, 
                   grouping_variables = c(grouping_variables),
                   reps_col = reps_col)
  
  ## write out score files
  purrr::walk2(scores, score_files, readr::write_csv)
  invisible(score_files)
}




