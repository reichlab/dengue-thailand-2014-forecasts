##### Functions for Formatting and Aggregating Thailand Dengue Counts and Forecasts
##### by Stephen Lauer, February 2015

#' Format counts
#' Loads in counts from a cntry.data.linelist object and converts them into a data frame
#'
#' @param den_data a cntry.data.linelist object with a count matrix and a linelist of counts
#'
#' @return counts a data frame with columns for pid (province), delivery_date, date_sick_year, date_sick_biweek, and number of counts
format_counts <- function(den_data){
    if("package:reshape" %in% search())
        detach("package:reshape", unload=TRUE)
    require(dengueThailand)
    require(dplyr)
    require(tidyr)
    require(lubridate)
    data(thai_prov_data)
    old_prov_ids <- data_frame(Province = c("Udon Thani", "Nong Bua Lam Phu",
                                            "Nakhon Phanom", "Mukdahan",
                                            "Ubon Ratchathani", "Amnat Charoen",
                                            "Prachin Buri", "Sa Kaeo"),
                               pid = c("TH76", "TH79",
                                       "TH73", "TH78",
                                       "TH75", "TH77",
                                       "TH74", "TH80"),
                               old_Province = c("Udon Thani", "Udon Thani",
                                                "Nakhon Phanom", "Nakhon Phanom",
                                                "Ubon Ratchathani", "Ubon Ratchathani",
                                                "Prachin Buri", "Prachin Buri"),
                               old_pid = c("TH19", "TH19", 
                                           "TH21", "TH21", 
                                           "TH71", "TH71",
                                           "TH45", "TH45"))
    
    ## find ISO numbers for Bueng Kan and Nong Khai for BK removal
    iso_BK <- thai_prov_data$ISO[which(thai_prov_data$Province == "Bueng Kan")]
    iso_NK <- thai_prov_data$ISO[which(thai_prov_data$Province == "Nong Khai")]
    
    ## use line listings to aggregate recent counts
    line_list <- den_data@line.list %>%
        filter(disease == 26, !is.na(date_sick)) %>% ## only use disease == 26 (DHF cases)
        mutate(delivery_date = as.Date(delivery_date),
               date_sick = as.Date(date_sick),
               province = ifelse(province == iso_BK, iso_NK, province), ## assign all of Bueng Kan's cases to Nong Khai
               pid = thai_prov_data$FIPS[match(province, thai_prov_data$ISO)],
               date_sick_year = year(date_sick),
               date_sick_biweek = date_to_biweek(date_sick))
    
    ## assign all counts a delivery date. For early counts, use 2011-04-09
    line_list$delivery_date[which(is.na(line_list$delivery_date))] <- as.Date("2011-04-09")
    
    ## sum all counts by pid, delivery_date, date_sick_biweek and date_sick_year
    agg_line_list <- line_list %>%
        group_by(pid, delivery_date, date_sick_year, date_sick_biweek) %>%
        summarise(count = n())
    
    ## extract prov_names and set all NA counts in den_data to 0
    prov_names <- row.names(den_data@.Data)
    den_df <- as.data.frame(den_data@.Data)
    den_df[is.na(den_df)] <- 0
    
    ## gather matrix of counts into a data frame
    ## only look at old counts and give delivery date.
    ## give old counts for split provinces old pids
    ## then bind with new counts
    counts <- den_df %>%
        mutate(province = prov_names,
               pid = old_prov_ids$old_pid[match(province,
                                                old_prov_ids$Province)],
               pid = ifelse(is.na(pid),
                            as.character(thai_prov_data$FIPS[match(province, 
                                                      thai_prov_data$Province)]),
                            pid)) %>%
        gather(biweek_col, count, starts_with("V")) %>%
        mutate(biweek_num = as.numeric(gsub("V", "", biweek_col)),
               date_sick_year = 1968 + floor((biweek_num - 1)/26),
               date_sick_biweek = ifelse(biweek_num %% 26 == 0, 26,
                                         biweek_num %% 26),
               delivery_date = as.Date("2011-04-09")) %>%
        group_by(pid, delivery_date, date_sick_year, date_sick_biweek) %>%
        summarise(count=sum(count, na.rm=T)) %>%
        filter(date_sick_year <= 1998) %>%
        bind_rows(agg_line_list)
    
    return(counts)
}


#' Outbreak level
#' The outbreak level is found by taking the average + 2 sd of the past eight years without the max year
#' @param x a numeric vector
#' @param stdevs the number of standard deviations above the mean to set the outbreak level
outbreak_level <- function(x, stdevs = 2){
    no_max <- x[-which.max(x)]
    floor(mean(no_max) + stdevs * sd(no_max)) + 1
}

#' Get epidemic thresholds
#' Calculate the outbreak thresholds for each province and biweek before a particular delivery date
#' @param counts a data frame of counts, like that outputted by format_counts()
#' @param deliv_date the most recent delivery_date to be used in epidemic threshold calculation
#' @param to_date_lag the number of biweeks in the past to be discarded in order to have only full data
#' @param trail_years number of years to use for threshold calculation
#' @param epi_thresh number of standard deviations above the average to use as the outbreak level
#'
#' @return epidemic_thresholds a number or numeric vector of outbreak thresholds
get_epidemic_thresholds <- function(counts,
                                    deliv_date,
                                    trail_years = 8,
                                    epi_thresh = 2) {
    require(lubridate)
    require(dplyr)
    require(tidyr)
    require(dengueThailand)
    
    ## if unspecified use the latest delivery date in counts file
    if(missing(deliv_date))
        deliv_date <- max(counts$delivery_date)
    deliv_date <- as.Date(deliv_date)
    
    ## the last_full_year is last year if the analysis date is after April, and two years ago otherwise
    last_full_year <- ifelse(month(deliv_date) > 4, year(deliv_date) - 1, year(deliv_date - 2))
    
    ## the from_year is ten years prior to the last_full_year
    from_year <- last_full_year - trail_years
    
    ## summarise counts (missing 0 count biweeks)
    recent_counts <- counts %>%
        filter(date_sick_year > from_year,
               date_sick_year <= last_full_year) %>%
        group_by(pid, date_sick_year, date_sick_biweek) %>%
        summarise(count = sum(count))
    
    ## add in all zero count biweeks
    all_provs_dates <- expand.grid(pid = unique(recent_counts$pid),
                                   date_sick_year = unique(recent_counts$date_sick_year),
                                   date_sick_biweek = unique(recent_counts$date_sick_biweek)) %>%
        left_join(recent_counts, by = c("pid", "date_sick_biweek", "date_sick_year"))
    all_provs_dates$count[is.na(all_provs_dates$count)] <- 0
    
    ## use outbreak_level to find the epi_thresholds; epi_thresh is # of std devs for threshold
    epidemic_thresholds <- all_provs_dates %>%
        group_by(pid, date_sick_biweek) %>%
        summarise(epidemic_threshold = outbreak_level(count, epi_thresh))
    
    return(epidemic_thresholds)
}

#' Format single forecast
#' Takes a cntry.data.linelist of forecasts and converts it into a data frame that can be used for forecast evaluation
#'
#' @param den_forecast a predicted.cntry.data object with an MC.sims slot
#' @param sim_results the MC.sims of a predicted.cntry.data object, can be used instead of den_forecast
#' @param analysis_date the date that analysis is run, if blank uses Sys.Date()
#' @param to_date_lag the number of biweeks in the past to be discarded in order to have only full data
#' @param counts a data frame of formatted counts, returned from format_counts()
#' @param den_data if missing(counts), used to derive counts object
#'
#' @return returns a data frame of formatted forecasts with columns for province (pid), biweek, year, analysis date, simulation number, predicted counts, epidemic threshold, and whether or not an outbreak is predicted
format_one_forecast <- function(den_forecast,
                                sim_results,
                                analysis_date,
                                to_date_lag,
                                counts,
                                den_data) {
    require(dengueThailand)
    require(spatialpred)
    require(lubridate)
    require(dplyr)
    require(tidyr)
    data(thai_prov_data)
    
    ## sim_results taken from either den_forecast or specified
    if(!missing(den_forecast))
        sim_results <- den_forecast@MC.sims
    
    if(missing(sim_results))
        stop("Current version requires a den_forecast object or sim_results object to be specified.")
    
    ## set if unspecified make analysis_date the current date
    if(missing(analysis_date))
        analysis_date <- Sys.Date()
    
    ## store and transform simulation results
    dimnames(sim_results) <- list(row.names(den_forecast@.Data), ## province names
                                  paste0(den_forecast@yr, "_", den_forecast@time.in.yr), ## timepoint of prediction
                                  paste0("sim_", 1:dim(sim_results)[3])) ## stochastic simulation number
    
    ## organize results into one row per province/date-sick-biweek/analysis_date combination
    ## multiple stochastic simulations per row
    prov_preds <- reshape2::melt(sim_results, varnames=c("pname", "predicted_time", "simulation"), value.name="pred_count") %>%
        separate(predicted_time,
                 into=c("date_sick_year", "date_sick_biweek"),
                 sep="_",
                 convert=TRUE) %>%
        mutate(pid = thai_prov_data$FIPS[match(pname, thai_prov_data$Province)],
               analysis_date = as.Date(analysis_date),
               analysis_biweek = date_to_biweek(analysis_date),
               analysis_year = year(analysis_date))
    
    ## get counts for epidemic threshold calculations
    if(missing(counts))
        counts <- format_counts(den_data = den_data)
    
    ## calculate epi thresholds
    epi_thresh <- get_epidemic_thresholds(counts,
                                          deliv_date = biweek_to_date(biweek = den_forecast@time.in.yr[to_date_lag],
                                                                      year = den_forecast@yr[to_date_lag]),
                                          epi_thresh = 2)
    
    ## add epi_thresholds and outbreak predictions to organized output
    prov_preds <- left_join(prov_preds, epi_thresh,
                            by = c("pid", "date_sick_biweek")) %>%
        mutate(pred_outbreak = ifelse(pred_count >= epidemic_threshold, 1, 0))
    return(prov_preds)
}

#' Find seasonal medians
#' Finds the seasonal medians of the past n years (or specified time frame) to use for evaluating forecasting models
#'
#' @param counts a data frame of formatted counts, returned from format_counts()
#' @param first_year the first year to use in finding seasonal medians
#' @param last_year the last year to use in finding seasonal medians
#' @param n_years if no first_year selected, can use this with last_year to
#'
#' @return a data frame with province (pid), biweek, and seasonal medians
find_seasonal_medians <- function(counts, first_year, last_year, n_years) {
    require(dplyr)
    require(tidyr)
    
    ## if first_year not specified, find using n_years
    if(missing(first_year))
        first_year <- last_year - n_years
    
    ## make a dataframe of all possible province, date_sick_biweek, date_sick_year combinations
    all_years <- seq(first_year, last_year)
    all_biweeks <- seq(1, 26)
    all_provs_all_years <- expand.grid(pid = unique(thai_prov_data$FIPS),
                                       date_sick_year = all_years,
                                       date_sick_biweek = all_biweeks)
    
    ## filter and summarise counts into recent biweeks and fill in 0s
    biweek_counts <- counts %>%
        filter(date_sick_year >= first_year, date_sick_year <= last_year) %>%
        group_by(pid, date_sick_year, date_sick_biweek) %>%
        summarise(count = sum(count))
    all_biweek_counts <- left_join(all_provs_all_years, biweek_counts,
                                   by = c("pid", "date_sick_year", "date_sick_biweek"))
    all_biweek_counts$count[is.na(all_biweek_counts$count)] <- 0
    
    ## find biweek_medians
    biweek_medians <- all_biweek_counts %>%
        group_by(pid, date_sick_biweek) %>%
        summarise(seasonal_median = median(count))
    return(biweek_medians)
}

#' Aggregate biweekly data
#' Combine forecasts, counts, seasonal medians, last_obs, and epidemic thresholds into one dataframe. If forecasts have columns for province (pid = FIPS code), date_sick_biweek, date_sick_year, delivery_date, and pred_count, then specify forecast_agg = "raw". If forecasts are from makeForecast (i.e. contain predicted_count, ub, and lb), then use forecast_agg = "simple"
#' @param forecasts a data frame of forecasts, generated by either run_one_year_forecasts() or from the makeForecast pipeline
#' @param counts a data frame of formatted counts, returned from format_counts()
#' @param seasonal_medians a data frame with province (pid), biweek, and seasonal medians from find_seasonal_medians()
#' @param forecast_agg use "raw" if using a data frame with all simulations (from run_one_year_forecasts()) or "simple, if using aggregated data frame
#' @param ana_date the analysis date for making the forecasts, if undefined is today's date
#' @param to_date_lag the number of biweeks in the past to be discarded in order to have only full data
#'
#' @return
aggregate_biweekly_data <- function(forecasts,
                                    counts,
                                    seasonal_medians,
                                    forecast_agg = c("raw", "simple"),
                                    analysis_year,
                                    ana_date,
                                    to_date_lag){
    require(dengueThailand)
    require(spatialpred)
    require(dplyr)
    require(tidyr)
    
    ## use today's date if ana_date is missing
    if(missing(ana_date))
        ana_date <- Sys.Date()
    
    ## use ana_date to find to_date
    ana_date <- as.Date(ana_date)
    
    ## use "raw" if forecast_agg undefined
    forecast_agg <- match.arg(forecast_agg)
    
    ## aggregate counts for all provinces, years, and biweeks and fill in 0s, and add last_obs value
    agg_counts <- counts %>%
        group_by(pid, date_sick_year, date_sick_biweek) %>%
        summarise(obs_count = sum(count, na.rm = TRUE))
    recent_counts <- left_join(expand.grid(pid = unique(counts$pid),
                                           date_sick_year = unique(counts$date_sick_year),
                                           date_sick_biweek = unique(counts$date_sick_biweek)),
                               agg_counts,
                               by = c("pid", "date_sick_year", "date_sick_biweek")) %>%
        mutate(obs_count = ifelse(is.na(obs_count), 0, obs_count))
    
    
    
    if(forecast_agg == "raw"){
        biweek_forecasts <- forecasts %>%
            filter(analysis_date <= ana_date,
                   date_sick_year == analysis_year) %>%
            left_join(recent_counts, by = c("pid", "date_sick_biweek", "date_sick_year")) %>%
            mutate(absolute_error = abs(pred_count - obs_count)) %>%
            group_by(pid, analysis_date, date_sick_year, date_sick_biweek, epidemic_threshold, obs_count) %>%
            summarise(pred_lb = quantile(pred_count, probs = .025, na.rm = TRUE),
                      pred_ub = quantile(pred_count, probs = .975, na.rm = TRUE),
                      pred_median = median(pred_count, na.rm = TRUE),
                      AE_lb = quantile(absolute_error, probs = .025, na.rm = TRUE),
                      AE_ub = quantile(absolute_error, probs = .975, na.rm = TRUE),
                      outbreak_prob = mean(pred_outbreak, na.rm = TRUE))
    }
    
    if(forecast_agg == "simple"){
        biweek_forecasts <- forecasts
        epi_thresholds <-
            get_epidemic_thresholds(counts = counts,
                                    deliv_date = min(biweek_forecasts$delivery_date),
                                    epi_thresh = 2, trail_years = 8)
        biweek_forecasts <- left_join(biweek_forecasts, epi_thresholds, by = c("pid", "date_sick_biweek"))
    }
    
    last_obs <- find_delivery_last_obs(biweek_forecasts$analysis_date, counts, to_date_lag)
    all_data <- biweek_forecasts %>%
        left_join(seasonal_medians, by = c("pid", "date_sick_biweek")) %>%
        left_join(last_obs, by = c("pid", "analysis_date")) %>%
        mutate(last_obs = ifelse(is.na(last_obs), 0, last_obs)) %>%
        mutate(obs_outbreak = as.numeric(obs_count >= epidemic_threshold),
               pred_covered = as.numeric(pred_lb <= obs_count & pred_ub >= obs_count),
               ## find the absolute errors for prediction, seasonal median, and last_obs
               AE_pred = abs(obs_count - pred_median),
               AE_seas = abs(obs_count - seasonal_median),
               AE_last_obs = abs(obs_count - last_obs),
               ## find how many steps ahead (or behind for now-casting) the predictions are
               step = date_sick_biweek + date_sick_year * 26 - date_to_biweek(analysis_date) - year(analysis_date) * 26)
    
    return(all_data)
}

#' Determine high season
#' Choose pct_cases to automatically each province's high season or choose season_start and season_end to set high season for all provinces yourself
#'
#' @param counts a data frame of formatted counts, returned from format_counts()
#' @param forecasts a data frame of formatted forecasts, generated by aggregate_biweekly_data()
#' @param analysis_year high season will be determined by years prior to analysis year
#' @param trail_years number of years to use when determining historical counts in high season
#' @param pct_cases the minimum percentage of cases to be enclosed by the high season, for automatic selection
#' @param season_start first biweek of high season, manually specified
#' @param season_end last biweek of high season, manually specified
#'
#' @return prov_high_season a data frame with rows for each province showing the start and end dates for the high season, the historical median, and total from the last full high season
determine_high_season <- function(counts,
                                  analysis_year,
                                  trail_years = 10,
                                  season_start,
                                  season_end,
                                  pct_cases = NULL){
    require(dplyr)
    require(tidyr)
    require(dengueThailand)
    require(lubridate)
    
    if(is.null(pct_cases) & missing(season_start))
        stop("Please set either pct_cases (to automatically determine high season) or season_start and season_end (to manually determine high season).")
    
    ## narrow counts to the ten years before each season's prediction
    last_full_year <- analysis_year - 1
    first_old_year <- last_full_year - trail_years
    
    ## get complete counts for trail_years
    historic_counts <- left_join(expand.grid(pid = unique(counts$pid),
                                             date_sick_year = unique(counts$date_sick_year),
                                             date_sick_biweek = unique(counts$date_sick_biweek)), counts,
                                 by = c("pid", "date_sick_year", "date_sick_biweek")) %>%
        filter(date_sick_year <= last_full_year, date_sick_year >= first_old_year) %>%
        group_by(pid, date_sick_year, date_sick_biweek) %>%
        summarise(obs_count = sum(count))
    historic_counts$obs_count[is.na(historic_counts$obs_count)] <- 0
    
    if(!is.null(pct_cases)){
        ## use means to determine peak, start and end of high season
        mean_annual_df <- historic_counts %>%
            group_by(pid, date_sick_biweek) %>%
            summarise(mean_count = mean(obs_count)) %>%
            group_by(pid) %>%
            mutate(historic_peak = date_sick_biweek[which.max(mean_count)],
                   annual_cases = cumsum(mean_count)) %>%
            group_by(pid, historic_peak) %>%
            summarise(start_biweek = date_sick_biweek[which(annual_cases >= ((1 - pct_cases) / 2) * annual_cases[length(annual_cases)])[1]],
                      end_biweek = date_sick_biweek[which(annual_cases >= (1 - (1 - pct_cases) / 2) * annual_cases[length(annual_cases)])[1]])
    }
    else {
        ## use means to determine peak date_sick_biweek then set start_biweek and end_biweek
        mean_annual_df <- historic_counts %>%
            group_by(pid, date_sick_biweek) %>%
            summarise(mean_count = mean(obs_count)) %>%
            group_by(pid) %>%
            summarise(historic_peak = date_sick_biweek[which.max(mean_count)],
                      start_biweek = season_start,
                      end_biweek = season_end)
    }
    ## find median high season counts and last high season's count
    prov_high_season <- left_join(historic_counts, mean_annual_df, by = "pid") %>%
        filter(date_sick_biweek >= start_biweek, date_sick_biweek <= end_biweek) %>%
        group_by(pid, date_sick_year, historic_peak, start_biweek, end_biweek) %>%
        summarise(annual_total = sum(obs_count)) %>%
        group_by(pid, historic_peak, start_biweek, end_biweek) %>%
        summarise(historic_median = median(annual_total),
                  last_season_count = annual_total[date_sick_year == last_full_year])
    
    return(prov_high_season)
}

#' Find delivery last_obs
#' Used within aggregate_high_season_data()
#'
#' @param all_data a data frame with counts and forecasts from one high season
#' @param counts a formatted counts data frame from format_counts()
#' @param to_date_lag the number of biweeks in the past to be discarded in order to have only full data
find_delivery_last_obs <- function(ana_dates, counts, to_date_lag = 4) {
    ana_date <- as.Date(unique(ana_dates))
    
    last_obs_data <- data.frame()
    for(i in 1:length(ana_date)) {
        last_obs_date <- biweek_to_date(date_to_biweek(ana_date[i]) -
                                            to_date_lag - 1,
                                        year(ana_date[i]))
        last_obs_biweek <- date_to_biweek(last_obs_date)
        last_obs_year <- year(last_obs_date)
        
        delivery_counts <- counts %>%
            filter(date_sick_biweek == last_obs_biweek, date_sick_year == last_obs_year, delivery_date <= ana_date[i]) %>%
            group_by(pid) %>%
            summarise(last_obs = sum(count)) %>%
            mutate(analysis_date = ana_date[i])
        
        last_obs_data <- bind_rows(last_obs_data, delivery_counts)
    }
    return(last_obs_data)
}

#' Aggregate high season data
#' Pulls together high season forecasts, historical counts, and observed counts for a given high season. If prov_high_season is specified, dont need to choose trail_years, pct_cases, season_start, or season_end (though specifying season_start and season_end can speed up process)
#'
#' @param forecasts forecasts a data frame of formatted forecasts, generated by aggregate_biweekly_data()
#' @param counts a data frame of formatted counts, returned from format_counts()
#' @param prov_high_season a data frame summarizing the high season from determine_high_seas
#' @param deliv_date the delivery date used when drawing the counts
#' @param ana_date the analysis date for making the forecasts, if undefined is today's date
#' @param to_date_lag the number of biweeks in the past to be discarded in order to have only full data
#' @param trail_years see determine_high_season()
#' @param pct_cases see determine_high_season()
#' @param season_start see determine_high_season()
#' @param season_end see determine_high_season()
#' @param analysis_year see determine_high_season()
#'
#' @return a data frame with high season forecasts, counts, historical medians, peak info, and MAEs
aggregate_high_season_data <- function(forecasts,
                                       counts,
                                       prov_high_season,
                                       to_date_lag,
                                       analysis_year,
                                       ana_date,
                                       trail_years = 10,
                                       season_start,
                                       season_end,
                                       pct_cases = NULL) {
    if(missing(prov_high_season)) {
        prov_high_season <- determine_high_season(counts,
                                                  forecasts,
                                                  analysis_year,
                                                  trail_years,
                                                  pct_cases,
                                                  season_start,
                                                  season_end)
    }
    if(missing(ana_date)){
        forecasts$analysis_date <- as.Date(forecasts$analysis_date)
        ana_date <- max(forecasts$analysis_date)
    }
    ana_date <- as.Date(ana_date)
    
    if(missing(to_date_lag))
        stop("Please specify to_date_lag")
    
    ## join forecasts and prov_high_season and find predicted peak and annual count for each simulation
    ## can save time if pct_cases not specified
    if(is.null(pct_cases) & !missing(season_start) & !missing(season_end)){
        high_season_forecasts <- forecasts %>%
            filter(date_sick_year == analysis_year,
                   analysis_biweek <= season_start + to_date_lag, ## only look at forecasts that occurred before the start of the season
                   date_sick_biweek >= season_start,
                   date_sick_biweek <= season_end) %>%
            left_join(prov_high_season, by = c("pid")) %>%
            group_by(pid, analysis_date, date_sick_year, start_biweek, end_biweek, historic_peak, historic_median, last_season_count, simulation) %>%
            summarise(pred_peak = date_sick_biweek[which.max(pred_count)],
                      pred_annual_count = sum(pred_count))
    } else{
        high_season_forecasts <- forecasts %>%
            left_join(prov_high_season, by = c("pid")) %>%
            filter(date_sick_year == analysis_year,
                   analysis_biweek <= start_biweek + to_date_lag, ## only look at forecasts that occurred before the start of the season
                   date_sick_biweek >= start_biweek,
                   date_sick_biweek <= end_biweek) %>%
            group_by(pid, analysis_date, date_sick_year, start_biweek, end_biweek, historic_peak, historic_median, last_season_count, simulation) %>%
            summarise(pred_peak = date_sick_biweek[which.max(pred_count)],
                      pred_annual_count = sum(pred_count))
    }
    
    ## join these with the number of counts we observed for the season
    annual_counts <- left_join(counts, prov_high_season, by = "pid") %>%
        filter(!is.na(start_biweek),
               date_sick_year == analysis_year,
               date_sick_biweek >= start_biweek,
               date_sick_biweek <= end_biweek) %>%
        group_by(pid, date_sick_year, date_sick_biweek) %>%
        summarise(obs_count = sum(count)) %>%
        group_by(pid, date_sick_year) %>%
        summarise(obs_peak = date_sick_biweek[which.max(obs_count)],
                  obs_annual_count = sum(obs_count))
    
    all_data <- left_join(high_season_forecasts, annual_counts, by = c("pid", "date_sick_year")) %>%
        mutate(AE_annual_count = abs(obs_annual_count - pred_annual_count)) %>%
        group_by(pid, analysis_date, date_sick_year, start_biweek, end_biweek, historic_peak, historic_median, last_season_count, obs_peak,
                 obs_annual_count) %>%
        summarise(pred_peak_median = median(pred_peak),
                  pred_peak_lb = quantile(pred_peak, probs = .025),
                  pred_peak_ub = quantile(pred_peak, probs = .975),
                  pred_total_median = median(pred_annual_count),
                  pred_total_lb = quantile(pred_annual_count, .025),
                  pred_total_ub = quantile(pred_annual_count, .975),
                  AE_lb = quantile(AE_annual_count, .025),
                  AE_ub = quantile(AE_annual_count, .975))
    
    ## add counts from last_obs from last_full_biweek before analysis_date
    last_obs_data <- find_delivery_last_obs(all_data$analysis_date, counts, to_date_lag)
    all_data <- all_data %>%
        left_join(last_obs_data, by = c("pid", "analysis_date")) %>%
        mutate(last_obs_pred = ifelse(is.na(last_obs), 0, last_obs) * (end_biweek - start_biweek + 1),
               pred_peak_covered = ifelse(pred_peak_lb <= obs_peak & pred_peak_ub >= obs_peak, 1,0),
               pred_total_covered = ifelse(pred_total_lb <= obs_annual_count & pred_total_ub >= obs_annual_count, 1,0),
               ## find the absolute errors for prediction, seasonal median, and last_obs
               AE_pred = abs(obs_annual_count - pred_total_median),
               AE_last_season = abs(obs_annual_count - last_season_count),
               AE_med = abs(obs_annual_count - historic_median),
               AE_last_obs = abs(obs_annual_count - last_obs_pred),
               biweeks_ahead = start_biweek + date_sick_year * 26 - date_to_biweek(analysis_date) - year(analysis_date) * 26)
    return(all_data)
}

#' Prepare biweek evaluations
#' A function that runs format_counts(), find_seasonal_medians(), and aggregate_biweekly_data, to prepare forecasts for evaluation
#'
#' @param den_data a cntry.data.linelist object with a count matrix and a linelist of counts
#' @param forecasts a data frame of forecasts, generated by either run_one_year_forecasts() or from the makeForecast pipeline
#' @param deliv_date the delivery date used when drawing the counts
#' @param to_date_lag the number of biweeks in the past to be discarded in order to have only full data
#' @param analysis_year year of interest
#' @param trail_years the number of past years to use for seasonal medians
#'
#' @return a list of objects with formatted counts, formatted forecasts, seasonal medians, and biweekly data
prepare_biweek_evaluations <- function(den_data,
                                       forecasts,
                                       deliv_date,
                                       to_date_lag,
                                       analysis_year,
                                       trail_years) {
    formatted_data <- format_counts(den_data)
    
    if(is.data.frame(forecasts)){
        formatted_forecasts <- forecasts
    } else{
        stop("forecasts needs to be a data frame")
    }
    
    seasonal_medians <- find_seasonal_medians(counts = formatted_data,
                                              first_year = analysis_year - 1 - trail_years,
                                              last_year = analysis_year - 1)
    
    biweekly_data <- aggregate_biweekly_data(forecasts = formatted_forecasts,
                                             counts = formatted_data,
                                             seasonal_medians = seasonal_medians,
                                             forecast_agg = "raw",
                                             analysis_year = analysis_year,
                                             to_date_lag = to_date_lag)
    return(list(formatted_data = formatted_data, formatted_forecasts = formatted_forecasts,
                seasonal_medians = seasonal_medians, biweekly_data = biweekly_data))
}

#' Subset den_data
#' Subset data by delivery date and MOPH regions of interest
#'
#' @param den_data a cntry.data.linelist object with a count matrix and a linelist of counts
#' @param deliv_date the delivery date used when drawing the counts
#' @param to_date_lag the number of biweeks in the past to be discarded in order to have only full data
#' @param moph_regions the MOPH regions of interest
#' @param complete_data whether or not to query based on complete data or real-time reported data
#'
#' @return a cntry.data.linelist object with a count matrix and a linelist of counts that has been subsetted
subset.den.data <- function(den_data,
                            deliv_date,
                            to_date_lag,
                            moph_regions,
                            complete_data=FALSE) {
    require(dengueThailand)
    require(lubridate)
    require(dplyr)
    
    if(missing(moph_regions))
        moph_regions <- seq(0, 12)
    if(missing(deliv_date))
        deliv_date <- Sys.Date()
    deliv_date <- as.Date(deliv_date)
    deliv_year <- year(deliv_date)
    deliv_biweek <- date_to_biweek(deliv_date)
    deliv_biweek_num <- deliv_year * 26 + deliv_biweek
    
    ## to_date is the first day of the t-l biweek,
    ## all cases delivered on or after this date should be excluded
    to_date <- biweek_to_date(deliv_biweek - to_date_lag, deliv_year)
    
    subset_pnames <- as.character(thai_prov_data$Province[which(thai_prov_data$MOPH_Admin_Code %in% moph_regions)])
    
    ## combine Bueng Kan with Nong Khai then remove Bueng Kan
    den_data@line.list$province <- ifelse(den_data@line.list$province == 38, 43, as.numeric(as.character(den_data@line.list$province)))
    if("Bueng Kan" %in% subset_pnames)
        subset_pnames <- subset_pnames[-which(subset_pnames =="Bueng Kan")]
    
    ## remove counts from data that have not been counted yet when:
    ## (1) more cases have arrived since analysis date
    ## AND
    ## (2) complete_data == FALSE
    ## ok to subset just by data prior to to_date when one of those conditions isn't met.
    if(complete_data==FALSE & length(which(as.Date(den_data@line.list$delivery_date) > deliv_date)) > 0) {
        ## isolate future counts
        future_line_list <- den_data@line.list %>%
            filter(disease == 26, as.Date(delivery_date) > deliv_date) %>%
            mutate(delivery_date = as.Date(delivery_date),
                   date_sick = as.Date(date_sick),
                   date_sick_year = year(date_sick),
                   date_sick_biweek = date_to_biweek(date_sick),
                   biweek_num = date_sick_year * 26 + date_sick_biweek,
                   pname = as.character(thai_prov_data$Province[match(province, thai_prov_data$ISO)])) %>%
            group_by(pname, date_sick_year, date_sick_biweek, biweek_num) %>%
            summarise(count = n()) %>%
            filter(biweek_num <= deliv_biweek_num)
        
        den_data_biweek_num <- den_data@yr * 26 + den_data@time.in.yr
        den_data_row_names <- row.names(den_data)
        deliv_data <- as.data.frame(den_data@.Data) %>%
            select(1:length(which(den_data_biweek_num <= deliv_biweek_num)))
        
        for(i in 1:dim(future_line_list)[1]){
            deliv_data[which(den_data_row_names == as.character(future_line_list[i, "pname"])),
                       which(den_data_biweek_num == as.numeric(as.character(future_line_list[i,"biweek_num"])))] <-
                deliv_data[which(den_data_row_names == as.character(future_line_list[i, "pname"])),
                           which(den_data_biweek_num == as.numeric(as.character(future_line_list[i,"biweek_num"])))] - future_line_list[i, "count"]
        }
        ## subset data to only include data prior to to_date
        den_data@.Data <- as.matrix(deliv_data[subset_pnames, 1:(dim(deliv_data)[2] - to_date_lag-1)])
        den_data@n.locs <- dim(den_data@.Data)[1]
        den_data@n.times <- dim(den_data@.Data)[2]
        den_data@t <- den_data@t[1:(dim(den_data@.Data)[2])]
        den_data@yr <- den_data@yr[1:(dim(den_data@.Data)[2])]
        den_data@time.in.yr <- den_data@time.in.yr[1:(dim(den_data@.Data)[2])]
        den_data@line.list$delivery_date[is.na(den_data@line.list$delivery_date)] <-
            min(as.Date(unique(den_data@line.list$delivery_date)), na.rm = TRUE)
        den_data@line.list <- filter(den_data@line.list, as.Date(delivery_date) <= deliv_date, as.Date(date_sick) < to_date)
    } else {
        ## subset data to only include data prior to to_date, explicit logic:
        ## (1) prior to year of to_date OR (2) <= year of to_date AND < biweek of to_date
        idx_to_keep <- which(den_data@yr<year(to_date) |
                                 (den_data@yr<=year(to_date) & den_data@time.in.yr<date_to_biweek(to_date)))
        ## den_data <- subset(den_data, t.criteria=idx_to_keep) ## This command wasn't making the right subset
        den_data@.Data <- den_data[subset_pnames, idx_to_keep]
        den_data@n.locs <- nrow(den_data@.Data)
        den_data@n.times <- ncol(den_data@.Data)
        den_data@t <- den_data@t[1:ncol(den_data@.Data)]
        den_data@yr <- den_data@yr[1:ncol(den_data@.Data)]
        den_data@time.in.yr <- den_data@time.in.yr[1:ncol(den_data@.Data)]
        #     den_data@line.list$delivery_date[is.na(den_data@line.list$delivery_date)] <-
        #         min(as.Date(unique(den_data@line.list$delivery_date)), na.rm = TRUE)
        #     den_data@line.list <- filter(den_data@line.list, as.Date(delivery_date) <= deliv_date, as.Date(date_sick) < to_date)
        ## not sure if it matters that linelist is not subset here.
    }
    return(den_data)
}

## run one full-year forecast
#' Run one full-year forecasts
#'
#' @param den_data a cntry.data.linelist object with a count matrix and a linelist of counts
#' @param analysis_date the date at which the analysis should be set to run at in YYYY-MM-DD format. Only data available before this date are included.
#' @param loc.fit.fn forecasting model to fit
#' @param to_date_lag the number of biweeks in the past to be discarded in order to have only full data
#' @param moph_regions the MOPH regions of interest
#' @param num.tops number of tops needed for spamd models
#' @param cor.lags number of lags needed for spamd models
#' @param stochastic whether to take stochastic steps
#' @param verbose should we print out verbose information?
#' @param MC_sims if we are taking stochastic steps, the number of MC chains, ignored if this is not a stochastic run.
#' @param predictions.only should we only return the predictions?
#' @param num.cores how many cores should be used to calculate the MC.sims
#' @param memory whether to include a memory term in the model fitting
#' @param fit.start.year first year to include in the analyzed data
#'
#' @return a data frame of all MC sims of all forecasts
run_one_full_year_forecast <- function(den_data,
                                       analysis_date,
                                       loc.fit.fn,
                                       to_date_lag,
                                       moph_regions,
                                       num.tops,
                                       cor.lags,
                                       stochastic,
                                       verbose,
                                       MC_sims,
                                       predictions.only,
                                       num.cores,
                                       memory=FALSE,
                                       fit.start.year){
    require(dplyr)
    require(tidyr)
    require(reshape2)
    require(dengueThailand)
    require(spatialpred)
    
    data(thai_prov_data)
    
    if(missing(moph_regions))
        moph_regions <- seq(0,12)
    
    ## run each forecast and bind them together
    all_MC_forecasts <- data.frame()
    
    ## get counts for epidemic threshold calculations
    counts <- format_counts(den_data)
    
    message(paste("starting forecasts for", analysis_date, "::", Sys.time()))
    den_subset <- subset.den.data(den_data,
                                  deliv_date = analysis_date,
                                  to_date_lag = to_date_lag,
                                  moph_regions = moph_regions)
    
    ## impute 26th biweek
    den_subset <- impute_final_biweek_26(den_subset)
    
    ## subset data to start with fit.start.year
    ## this removes the cntry.data.linelist
    if(!is.null(fit.start.year)) {
        den_subset <- subset(den_subset, t.criteria=(den_subset@t>=fit.start.year))
    }
    
    ## smooth subsetted data
    den_smooth <- smooth.cdata(den_subset)
    
    ## fit model
    den_mdl <- fit.cntry.pred.mdl.tasker(data = den_smooth,
                                         loc.fit.fn = loc.fit.fn,
                                         cor.lags = cor.lags,
                                         num.tops = num.tops,
                                         use.deltas = T,
                                         memory = memory,
                                         fit.start.year = fit.start.year,
                                         scale = "log",
                                         hold.outs=numeric(),
                                         family=poisson,
                                         ignore.self=F,
                                         verbose=verbose)
    
    ## calculate steps needed to reach end of year
    steps <- 27 - date_to_biweek(ymd(analysis_date)) + to_date_lag
    
    ## make forecasts
    den_forecast <- forecast(den_mdl,
                             as(den_subset, "cntry.data"), ## as() maybe not needed anymore
                             steps=steps,
                             stochastic=stochastic,
                             verbose=verbose,
                             MC.sims=MC_sims,
                             predictions.only=predictions.only,
                             num.cores=num.cores)
    
    ## store and transform simulation results
    sim_results <- den_forecast@MC.sims
    dimnames(sim_results) <- list(row.names(den_subset@.Data), ## province names
                                  paste0(den_forecast@yr, "_", den_forecast@time.in.yr), ## timepoint of prediction
                                  paste0("sim_", 1:MC_sims)) ## stochastic simulation number
    
    ## organize results into one row per province/date-sick-biweek/analysis-date combination
    ## multiple stochastic simulations per row
    prov_preds <- reshape2::melt(sim_results,
                                 varnames=c("pname", "predicted_time", "simulation"),
                                 value.name="pred_count") %>%
        spread(key=simulation, value=pred_count) %>%
        separate(predicted_time,
                 into=c("date_sick_year", "date_sick_biweek"),
                 sep="_",
                 convert=TRUE) %>%
        mutate(pid             = thai_prov_data$FIPS[match(pname, thai_prov_data$Province)],
               analysis_date   = as.Date(analysis_date, format = "%Y-%m-%d"),
               analysis_biweek = date_to_biweek(ymd(analysis_date)),
               analysis_year   = year(ymd(analysis_date)))
    
    ## calculate epi thresholds
    epi_thresh <- get_epidemic_thresholds(counts,
                                          deliv_date = as.Date(analysis_date, format = "%Y-%m-%d"),
                                          epi_thresh = 2)
    
    ## add epi_thresholds to organized output
    prov_preds <- left_join(prov_preds, epi_thresh,
                            by = c("pid", "date_sick_biweek"))
    
    ## combine into large data frame
    all_MC_forecasts <- bind_rows(all_MC_forecasts, prov_preds)
    
    ## tidy it into a dataframe by province, date_sick_biweek, date_sick_year, and simulation
    tidy_forecasts <- all_MC_forecasts %>%
        gather(key = simulation, value = "pred_count", starts_with("sim")) %>%
        mutate(pred_outbreak = ifelse(pred_count >= epidemic_threshold, 1, 0))
    return(forecasts = tidy_forecasts)
}


## run all forecasts
#' Run one year forecasts
#'
#' @param den_data a cntry.data.linelist object with a count matrix and a linelist of counts
#' @param loc.fit.fn forecasting model to fit
#' @param to_date_lag the number of biweeks in the past to be discarded in order to have only full data
#' @param moph_regions the MOPH regions of interest
#' @param num.tops number of tops needed for spamd models
#' @param cor.lags number of lags needed for spamd models
#' @param steps number of steps to forecast ahead
#' @param run_full_year if TRUE, steps will be determined by the number of biweeks to the end of the year
#' @param stochastic whether to take stochastic steps
#' @param verbose should we print out verbose information?
#' @param MC_sims if we are taking stochastic steps, the number of MC chains, ignored if this is not a stochastic run.
#' @param predictions.only should we only return the predictions?
#' @param num.cores how many cores should be used to calculate the MC.sims
#' @param first_deliv_date what is the earliest delivery date to consider running the model for
#' @param memory whether to include a memory term in the model fitting
#' @param fit.start.year first year to include in the analyzed data
#' @param max_steps the max number of steps to run
#' @param complete_data whether or not to run the forecasts on complete data
#'
#' @return a data frame of all MC sims of all forecasts
run_one_year_forecasts <- function(den_data,
                                   loc.fit.fn,
                                   to_date_lag,
                                   moph_regions,
                                   num.tops,
                                   cor.lags,
                                   analysis_year,
                                   in.year.analysis.dates,
                                   stochastic,
                                   verbose,
                                   MC_sims,
                                   predictions.only,
                                   num.cores,
                                   memory,
                                   fit.start.year,
                                   max_steps = 13,
                                   complete_data=FALSE,
                                   expand_counts,
                                   last_yr_complete_biweek=9,
                                   last_yr_active_biweek=2,
                                   smooth=c("spline", "loess", "none")
                                   ){
    if(complete_data & expand_counts)
        stop("Either expand counts or use complete data, but not both.")
    library(dplyr)
    library(tidyr)
    library(reshape2)
    library(dengueThailand)
    library(spatialpred)
    
    data(thai_prov_data)
    
    if(missing(moph_regions))
        moph_regions <- seq(0,12)
    
    ## define all analysis dates as first days of defined biweeks after a data delivery
    del_date_idx <- which(!is.na(den_data@line.list$delivery_date))
    all_deliv_dates <- sort(as.Date(unique(den_data@line.list[del_date_idx, "delivery_date"])))
    analysis_dates <- unique(biweek_to_date(date_to_biweek(all_deliv_dates) + 1,
                                            year(all_deliv_dates)))
    
    ## only keep delivery dates that make predictions within year of interest
    earliest_pred_date <- biweek_to_date(1+to_date_lag-max_steps, analysis_year)
    latest_pred_date <- biweek_to_date(26+to_date_lag-1,analysis_year)
    analysis_dates <- analysis_dates[analysis_dates <= latest_pred_date &
                                         analysis_dates >= earliest_pred_date]
    
    ## remove any analysis dates outside of the analysis year
    if(in.year.analysis.dates)
        analysis_dates <- analysis_dates[year(analysis_dates) == analysis_year]
    
    
    ## run each forecast and bind them together
    all_MC_forecasts <- data.frame()
    
    ## get counts for epidemic threshold calculations
    counts <- format_counts(den_data)
    
    for(i in 1:length(analysis_dates)){
        message(paste("starting forecasts for", analysis_dates[i], "::", Sys.time()))
        
        ## subset to counts available at analysis date
        if(complete_data) {
            ana_data <- subset.den.data(den_data,
                                        deliv_date = analysis_dates[i],
                                        to_date_lag = 0,
                                        moph_regions = moph_regions,
                                        complete_data = TRUE)
        } else {
            ana_data <- subset.den.data(den_data,
                                        deliv_date = analysis_dates[i],
                                        to_date_lag = 0,
                                        moph_regions = moph_regions,
                                        complete_data = FALSE)
        }
        
        analysis_biweek <- date_to_biweek(analysis_dates[i])
        
        ## expand recent case counts
        if(expand_counts){
            ana_data <- expand_cntry_data(cntry_data = ana_data,
                                          last_yr_complete_biweek = last_yr_complete_biweek,
                                          last_yr_active_biweek = last_yr_active_biweek,
                                          analysis_date = analysis_dates[i])
        }
        
        ## subset data according to to_date_lag
        if(complete_data) {
            den_subset <- subset.den.data(ana_data,
                                          deliv_date = analysis_dates[i],
                                          to_date_lag = to_date_lag,
                                          moph_regions = moph_regions,
                                          complete_data = TRUE)
        } else {
            den_subset <- subset.den.data(ana_data,
                                          deliv_date = analysis_dates[i],
                                          to_date_lag = to_date_lag,
                                          moph_regions = moph_regions,
                                          complete_data = FALSE)
        }
        ## impute 26th biweek
        den_subset <- impute_final_biweek_26(den_subset)
        
        ## subset data to start with fit.start.year
        ## this removes the cntry.data.linelist
        if(!is.null(fit.start.year)) {
            den_subset <- subset(den_subset, t.criteria=(den_subset@t>=fit.start.year))
        }
        
        ## smooth subsetted data
        if(smooth=="spline")
            den_smooth <- smooth.cdata(den_subset)
        # if(smooth=="loess")
        # add code for making a loess smooth for data
        if(smooth=="none")
            den_smooth <- den_subset
        
        ## fit model
        den_mdl <- fit.cntry.pred.mdl.tasker(data = den_smooth,
                                             loc.fit.fn = loc.fit.fn,
                                             cor.lags = cor.lags,
                                             num.tops = num.tops,
                                             use.deltas = T,
                                             memory = memory,
                                             fit.start.year = NULL,
                                             scale = "log",
                                             hold.outs=numeric(),
                                             family=poisson,
                                             ignore.self=F,
                                             verbose=verbose)
        
        ## calculate steps needed to reach end of year
        steps <- min(max_steps, 
                     analysis_year*26 + 1 - (date_to_biweek(analysis_dates[i]) + 
                                                 26*(year(analysis_dates[i])-1) -
                                                 to_date_lag))
        
        ## make forecasts
        den_forecast <- forecast(den_mdl,
                                 as(den_subset, "cntry.data"), ## as() maybe not needed anymore
                                 steps=steps,
                                 stochastic=stochastic,
                                 verbose=verbose,
                                 MC.sims=MC_sims,
                                 predictions.only=predictions.only,
                                 num.cores=num.cores)
        
        ## store and transform simulation results
        sim_results <- den_forecast@MC.sims
        dimnames(sim_results) <- list(row.names(den_subset@.Data), ## province names
                                      paste0(den_forecast@yr, "_", den_forecast@time.in.yr), ## timepoint of prediction
                                      paste0("sim_", 1:MC_sims)) ## stochastic simulation number
        
        ## organize results into one row per province/date-sick-biweek/analysis-date combination
        ## multiple stochastic simulations per row
        prov_preds <- reshape2::melt(sim_results,
                                     varnames=c("pname", "predicted_time", "simulation"),
                                     value.name="pred_count") %>%
            spread(key=simulation, value=pred_count) %>%
            separate(predicted_time,
                     into=c("date_sick_year", "date_sick_biweek"),
                     sep="_",
                     convert=TRUE) %>%
            mutate(pid             = thai_prov_data$FIPS[match(pname, thai_prov_data$Province)],
                   analysis_date   = as.Date(analysis_dates[i], format = "%Y%m%d"),
                   analysis_biweek = date_to_biweek(analysis_date),
                   analysis_year   = analysis_year)
        
        ## calculate epi thresholds
        epi_thresh <- get_epidemic_thresholds(counts,
                                              deliv_date = as.Date(analysis_dates[i], format = "%Y%m%d"),
                                              epi_thresh = 2)
        
        ## add epi_thresholds to organized output
        prov_preds <- left_join(prov_preds, epi_thresh,
                                by = c("pid", "date_sick_biweek"))
        
        ## combine into large data frame
        all_MC_forecasts <- bind_rows(all_MC_forecasts, prov_preds)
    }
    ## tidy it into a dataframe by province, date_sick_biweek, date_sick_year, and simulation
    tidy_forecasts <- all_MC_forecasts %>%
        gather(key = simulation, value = "pred_count", starts_with("sim")) %>%
        mutate(pred_outbreak = ifelse(pred_count >= epidemic_threshold, 1, 0))
    return(forecasts = tidy_forecasts)
}


#' quick imputation of biweek 26 for gnarly underreporting
#'
#' @param d a den.data object
#'
#' @return a den.data object with most recent biweek 26 imputed in a simple way.
impute_final_biweek_26 <- function(d) {
    require(spatialpred)
    last_biweek_26_idx <- max(which(d@time.in.yr == 26))
    
    ## is biweek 26 the last observation?
    if(last_biweek_26_idx == ncol(d)) {
        ## if yes, impute biweek 25's value (unless biweek 26 larger)
        d[,last_biweek_26_idx] <- pmax(d[,last_biweek_26_idx-1],
                                      d[,last_biweek_26_idx])
    } else {
        ## if no, impute from 25 - 26 - 1 (unless biweek 26 larger)
        d[,last_biweek_26_idx] <- pmax(round((d[,last_biweek_26_idx-1] + d[,last_biweek_26_idx+1])/2),
                                      d[,last_biweek_26_idx])
    }
    return(d)
}


#' Run forecast on latest observed data
#'
#' @param den_data a cntry.data.linelist object with a count matrix and a linelist of counts
#' @param loc.fit.fn forecasting model to fit
#' @param to_date_lag the number of biweeks in the past to be discarded in order to have only full data
#' @param moph_regions the MOPH regions of interest
#' @param num.tops number of tops needed for spamd models
#' @param cor.lags number of lags needed for spamd models
#' @param steps number of steps to forecast ahead
#' @param run_full_year if TRUE, steps will be determined by the number of biweeks to the end of the year
#' @param stochastic whether to take stochastic steps
#' @param verbose should we print out verbose information?
#' @param MC_sims if we are taking stochastic steps, the number of MC chains, ignored if this is not a stochastic run.
#' @param predictions.only should we only return the predictions?
#' @param num.cores how many cores should be used to calculate the MC.sims
#' @param memory whether to include a memory term in the model fitting
#' @param fit.start.year first year to include in the analyzed data
#' @param max_steps the max number of steps to run
#'
#' @return a data frame of all of the MC sims from this forecast
run_latest_forecast <- function(den_data,
                                loc.fit.fn,
                                to_date_lag,
                                moph_regions,
                                num.tops,
                                cor.lags,
                                steps,
                                run_full_year,
                                stochastic,
                                verbose,
                                MC_sims,
                                predictions.only,
                                num.cores,
                                memory,
                                fit.start.year,
                                expand_counts,
                                last_yr_complete,
                                last_yr_active,
                                last_yr_active_biweek,
                                smooth=c("spline", "loess", "none")){
    if("package:reshape" %in% search())
        detach("package:reshape", unload=TRUE)
    library(reshape2)
    library(dplyr)
    library(tidyr)
    library(dengueThailand)
    library(spatialpred)
    data(thai_prov_data)
    
    if(missing(moph_regions))
        moph_regions <- seq(0,12)
    
    ## define analysis date as first day of the defined biweek after the last data delivery
    last_deliv_date <- max(den_data@line.list$delivery_date, na.rm=TRUE)
    analysis_date <- biweek_to_date(date_to_biweek(last_deliv_date) + 1, year(last_deliv_date))
    
    ## expand recent case counts
    if(expand_counts){
        source("trunk/source/realtime-models/dynamic-expansion-for-cntry-data.R")
        den_data <- expand_cntry_data(cntry_data = den_data,
                                      last_yr_complete = last_yr_complete,
                                      last_yr_active = last_yr_active,
                                      last_yr_active_biweek = last_yr_active_biweek,
                                      analysis_date = analysis_date)
    }
    
    
    ## subset data to only the counts assumed to be full
    den_subset <- subset.den.data(den_data,
                                  deliv_date = analysis_date,
                                  to_date_lag = to_date_lag,
                                  moph_regions = moph_regions)
    
    ## impute 26th biweek
    den_subset <- impute_final_biweek_26(den_subset)
    
    ## smooth subsetted data
    if(smooth=="spline")
        den_smooth <- smooth.cdata(den_subset)
    # if(smooth=="loess")
        # add code for making a loess smooth for data
    if(smooth=="none")
        den_smooth <- den_subset
    
    ## fit model
    den_mdl <- fit.cntry.pred.mdl.tasker(data = den_smooth,
                                         loc.fit.fn = loc.fit.fn,
                                         cor.lags = cor.lags,
                                         num.tops = num.tops,
                                         use.deltas = T,
                                         memory = memory,
                                         fit.start.year = fit.start.year,
                                         scale = "log",
                                         hold.outs = numeric(),
                                         family = poisson,
                                         ignore.self = F,
                                         verbose = verbose)
    
    ## calculate steps needed to reach end of year
    if(run_full_year)
        steps <- 26 - date_to_biweek(last_deliv_date) + to_date_lag
    ## make forecasts
    
    den_forecast <- forecast(den_mdl,
                             as(den_subset, "cntry.data"),
                             steps=steps,
                             stochastic=stochastic,
                             verbose=verbose,
                             MC.sims=MC_sims,
                             predictions.only=predictions.only,
                             num.cores=num.cores)
    
    ## store and transform simulation results
    sim_results <- den_forecast@MC.sims
    dimnames(sim_results) <- list(row.names(den_smooth@.Data), ## province names
                                  paste0(den_forecast@yr, "_", den_forecast@time.in.yr), ## timepoint of prediction
                                  paste0("sim_", 1:MC_sims)) ## stochastic simulation number
    
    ## organize results into one row per province/date-sick-biweek combination
    ## multiple stochastic simulations per row
    prov_preds <- melt(sim_results, varnames=c("pname", "predicted_time", "simulation"), value.name="pred_count") %>%
        separate(predicted_time,
                 into=c("date_sick_year", "date_sick_biweek"),
                 sep="_",
                 convert=TRUE) %>%
        mutate(pid = thai_prov_data$FIPS[match(pname, thai_prov_data$Province)],
               analysis_date = analysis_date)
    
    ## get counts for epidemic threshold calculations
    counts <- format_counts(den_subset)
    
    ## calculate epi thresholds
    epi_thresh <- get_epidemic_thresholds(counts,
                                          deliv_date = as.Date(last_deliv_date, format = "%Y%m%d"),
                                          epi_thresh = 2)
    
    ## add epi_thresholds and outbreak predictions to organized output
    prov_preds <- left_join(prov_preds, epi_thresh, by = c("pid", "date_sick_biweek")) %>%
        mutate(pred_outbreak = ifelse(pred_count >= epidemic_threshold, 1, 0))
    
    return(list(forecasts = prov_preds, counts = counts))
}