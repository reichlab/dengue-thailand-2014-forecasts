##### Functions for Evaluating Thailand Dengue Forecasts
##### by Stephen Lauer, February 2015

#' Get model coefficients
#'
#' @param den_mdl a cntry.pred.mdl object
#' @param den_smooth a cntry.data.linelist, usually result of smooth.cdata()
#'
#' @return prov_models a data frame with model coefficients for each province
get_model_coefficients <- function(den_mdl,
                                   den_smooth){
  require(spatialpred)
  require(dplyr)
  prov_models <- c()
  for(i in 1:nrow(den_smooth@.Data)){
    prov_terms <- as.data.frame(t(round(den_mdl@loc.mdls[[i]]@mdl$coefficients, 1)))
    prov_terms$province <- row.names(den_smooth@.Data)[i]
    prov_models <- bind_rows(prov_models, prov_terms)
  }
  return(prov_models)
}

#' Summarize fit model
#' Find the model residuals of the cntry.pred.mdl object used for forecasting
#'
#' @param den_mdl a cntry.pred.mdl object
#' @param den_smooth a cntry.data.linelist, usually result of smooth.cdata()
#' @param deliv_date the delivery date used when drawing the counts
#' @param start_biweek the first biweek of the high season
#' @param end_biweek the last biweek of the high season
#'
#' @return a list of three data frames: for overall model residuals, low season residuals, and high season residuals
summarize_model <- function(den_mdl,
                            den_smooth,
                            deliv_date,
                            start_biweek,
                            end_biweek){
  require(spatialpred)
  require(dengueThailand)
  require(dplyr)
  
  if(missing(deliv_date))
    deliv_date <- max(as.Date(unique(den_smooth@line.list$delivery_date)), na.rm=TRUE)
  
  data(thai_prov_data)
  
  ## find model residuals for each province and bind them together
  model_data <- c()
  for(i in 1:length(den_mdl@loc.mdls)){
    mdl_dat <- den_mdl@loc.mdls[[i]]@data %>%
      rename_(obs_counts = "y", date_sick_biweek = "time.in.yr", date_sick_year = "yr")
    mdl_dat$fit_cases <- den_mdl@loc.mdls[[i]]@mdl$fitted.values
    mdl_dat$residuals <- mdl_dat$fit_cases - mdl_dat$obs_counts
    mdl_dat$province <- row.names(den_smooth@.Data)[i]
    mdl_dat$delivery_date <- deliv_date
    mdl_dat$std_res <- scale(mdl_dat$residuals)
    model_data <- bind_rows(model_data, mdl_dat)
  }
  
  ## NGR: this piped set of operations is the source of an error
  ## find low season residuals
  low_season_data <- model_data %>%
    filter(date_sick_biweek < start_biweek | date_sick_biweek > end_biweek) %>%
    mutate(date_sick_year = ifelse(date_sick_biweek > end_biweek, date_sick_year + 1, date_sick_year)) %>%
    group_by(province, date_sick_year, delivery_date) %>%
    summarise(annual_count = sum(obs_counts),
              fit_count = sum(fit_cases),
              residuals = fit_count - annual_count) %>%
    group_by(province) %>%
    mutate(std_res = residuals/sd(residuals),
           prov_quartile = ntile(annual_count,4))
  
  ## find high season residuals
  high_season_data <- model_data %>%
    filter(date_sick_biweek >= start_biweek, date_sick_biweek <= end_biweek) %>%
    group_by(province, date_sick_year, delivery_date) %>%
    summarise(annual_count = sum(obs_counts),
              fit_count = sum(fit_cases),
              residuals = fit_count - annual_count) %>%
    group_by(province) %>%
    mutate(std_res = residuals/sd(residuals),
           prov_quartile = ntile(annual_count,4))

  return(list(model_data = model_data,
              low_season_data = low_season_data,
              high_season_data = high_season_data))
}


#' High season outbreaks
#' Find the number of outbreaks that occurred in each province over the predicted high season
#'
#' @param biweekly_data aggregated biweekly forecasts and counts, output by aggregate_biweekly_data()
#' @param prov_high_season aggregated high season data, output by determine_high_season()
#' @param seasons season(s) of interest
#' @param to_date_lag the number of biweeks in the past to be discarded in order to have only full data
#'
#' @return a data frame of outbreaks by province
high_season_outbreaks <- function(biweekly_data,
                                  prov_high_season,
                                  seasons,
                                  to_date_lag) {
  require(tidyr)
  require(dplyr)
  
  number_of_outbreaks  <- left_join(biweekly_data, prov_high_season, by = "pid") %>%
    group_by(pid, start_biweek, end_biweek) %>%
    filter(date_sick_biweek >= start_biweek,
           date_sick_biweek <= end_biweek,
           step >= -to_date_lag,
           step < 26 - to_date_lag,
           date_sick_year %in% seasons,
           year(analysis_date) == date_sick_year) %>%
    group_by(pid, date_sick_year, date_sick_biweek, start_biweek, end_biweek, historic_peak) %>%
    summarise(outbreak = round(mean(obs_outbreak))) %>%
    group_by(pid, date_sick_year, start_biweek, end_biweek, historic_peak) %>%
    summarise(num_outbreaks = sum(outbreak))
  return(number_of_outbreaks)
}

#' Calculate outbreak sensitivity
#' find the sensitivity for provinces with outbreaks during the high season
#'
#' @param biweekly_data biweekly_data aggregated biweekly forecasts and counts, output by aggregate_biweekly_data()
#' @param number_of_outbreaks a data frame of outbreaks by province, as output by high_season_outbreaks()
#' @param steps_out the number of steps before the biweek that experienced outbreak
#' @param pred_outbreak_cutoff cutoff for outbreak_prob to make a binary outbreak prediction, default is 0.5
#'
#' @return data frame with the sensitivity for predicting outbreaks for each province in a given season
calc_sensitivity <- function(biweekly_data,
                             number_of_outbreaks,
                             steps_out = 0,
                             pred_outbreak_cutoff = 0.5) {
  require(tidyr)
  require(dplyr)
  
  sensitivity_df <- left_join(biweekly_data, number_of_outbreaks, by = c("pid", "date_sick_year")) %>%
    filter(date_sick_biweek >= start_biweek, ## only look at high season biweeks
           date_sick_biweek <= end_biweek,
           obs_outbreak == 1, ## only look at biweeks with an outbreak
           step %in% steps_out) %>%
    group_by(pid, date_sick_year, step) %>%
    summarise(n = n(),
              pred_outbreak_prob = mean(outbreak_prob),
              true_positive = sum(outbreak_prob >= pred_outbreak_cutoff),
              false_negative = sum(outbreak_prob < pred_outbreak_cutoff)) %>%
    mutate(sensitivity = true_positive / n)
  return(sensitivity_df)
}

## find the specificity for the provinces without high season outbreaks during their historic_peak
#' Calc specificity
#' Calculate the specificity for biweeks without outbreaks in each province
#'
#' @param biweekly_data biweekly_data aggregated biweekly forecasts and counts, output by aggregate_biweekly_data()
#' @param number_of_outbreaks a data frame of outbreaks by province, as output by high_season_outbreaks()
#' @param steps_out the number of steps before the biweek that experienced outbreak
#' @param pred_outbreak_cutoff cutoff for outbreak_prob to make a binary outbreak prediction, default is 0.5
#'
#' @return data frame with the sensitivity for outbreak-free biweeks for each province in a given season
calc_specificity <- function(biweekly_data,
                             number_of_outbreaks,
                             steps_out = 0,
                             pred_outbreak_cutoff = 0.5) {
  require(tidyr)
  require(dplyr)
  
  specificity_df <- left_join(biweekly_data, number_of_outbreaks, by = c("pid", "date_sick_year")) %>%
    filter(date_sick_biweek >= start_biweek, ## only look at high season biweeks
           date_sick_biweek <= end_biweek,
           obs_outbreak == 0, ## only look at provinces without outbreaks
           step %in% steps_out) %>%
    group_by(pid, date_sick_year, step) %>%
    summarise(n = n(),
              pred_outbreak_prob = mean(outbreak_prob),
              true_negative = sum(outbreak_prob < pred_outbreak_cutoff),
              false_positive = sum(outbreak_prob >= pred_outbreak_cutoff)) %>%
    mutate(specificity = true_negative / n)
  return(specificity_df)
}

#' Evaluate outbreak predictions
#'
#' @param biweekly_data biweekly_data aggregated biweekly forecasts and counts, output by aggregate_biweekly_data()
#' @param number_of_outbreaks a data frame of outbreaks by province, as output by high_season_outbreaks()
#' @param steps_out the number of steps before the biweek that experienced outbreak
#' @param pred_outbreak_cutoff cutoff for outbreak_prob to make a binary outbreak prediction, default is 0.5
#'
#' @return a data frame with the following columns:
#' \item{pid} province id
#' \item{date_sick_year} year that outbreaks occurred in
#' \item{step} the number of steps before the biweek that experienced outbreak
#' \item{n} the number of predictions
#' \item{num_outbreaks} the number of outbreaks
#' \item{mean_raw_sensitivity} the mean of the observation probabilities for biweeks that had an outbreak
#' \item{mean_raw_specificity} the mean of one minus the observation probabilities for biweeks that didn't have an outbreak
#' \item{true_positive} the number of biweeks with an outbreak in which the observation probability was greater than pred_outbreak_cutoff
#' \item{false_negative} the number of biweeks with an outbreak in which the observation probability was less than pred_outbreak_cutoff
#' \item{true_negative} the number of biweeks without an outbreak in which the observation probability was less than pred_outbreak_cutoff
#' \item{false_positive} the number of biweeks without an outbreak in which the observation probability was greater than pred_outbreak_cutoff
#' \item{sensitivity} the number of true positives predicted divided by the number of positives observed
#' \item{specificity} the number of true negatives predicted divided by the number of negatives observed
#' \item{pos_pred_value} the number of true positives predicted divided by the total number of positives predicted
#' \item{neg_pred_value} the number of true negatives predicted divided by the total number of negatives predicted
eval_outbreak_preds <- function(biweekly_data,
                                number_of_outbreaks,
                                steps_out = 0,
                                pred_outbreak_cutoff = 0.5) {
  require(tidyr)
  require(dplyr)
  
  outbreak_df <- left_join(biweekly_data, number_of_outbreaks, by = c("pid", "date_sick_year")) %>%
    filter(date_sick_biweek >= start_biweek, ## only look at high season biweeks
           date_sick_biweek <= end_biweek,
           step %in% steps_out) %>%
    group_by(pid, date_sick_year, step) %>%
    summarise(n = n(),
              num_outbreaks = sum(obs_outbreak == 1),
              mean_raw_sensitivity = mean(outbreak_prob[obs_outbreak == 1]),
              mean_raw_specificity = mean(1 - outbreak_prob[obs_outbreak == 0]),
              true_positive = sum(outbreak_prob >= pred_outbreak_cutoff & obs_outbreak == 1),
              false_negative = sum(outbreak_prob < pred_outbreak_cutoff & obs_outbreak == 1),
              true_negative = sum(outbreak_prob < pred_outbreak_cutoff & obs_outbreak == 0),
              false_positive = sum(outbreak_prob >= pred_outbreak_cutoff & obs_outbreak == 0)) %>%
    mutate(sensitivity = true_positive / (true_positive + false_negative),
           specificity = true_negative / (true_negative + false_positive),
           pos_pred_value = true_positive / (true_positive + false_positive),
           neg_pred_value = true_negative / (true_negative + false_negative))
  return(outbreak_df)
}

#' Calculate biweekly relative mean absolute error
#'
#' @param biweekly_data biweekly_data aggregated biweekly forecasts and counts, output by aggregate_biweekly_data()
#' @param denom either "seasonal" to compare to seasonal medians or "last_obs" to compare to the value from the last step
#' @param min_obs minimum number of observations for computing relative MAE (i.e., the rel_MAE of 1 observation may not be interesting)
#' @param vars variables to subset data by, which include "province", "step", "analysis_date", and "date_sick_biweek", choose two at most
#'
#' @return a data frame sorted by the specified vars showing the number of observations, the average observation, the average prediction (and upper and lower bounds for prediction) by the forecasting model, the average prediction by the specified denom, and the relative MAE as well as the upper and lower bounds for the MAE.
calculate_biweekly_rel_MAE <- function(biweekly_data,
                                       denom = c("seasonal", "last_obs"),
                                       min_obs = 10,
                                       vars = c("province", "step", "analysis_date", "date_sick_biweek")){
  ## if no denom is chosen, use "seasonal"
  denom <- match.arg(denom)
  
  if(denom == "seasonal")
    denom <- "seas"
  
  ## choose at most two variables to subset by
  if(length(vars) >= 3)
    stop('Choose at most two of the following for vars: "province", "step", "analysis_date", and "date_sick_biweek"')
  
  if("province" %in% vars)
    vars[which(vars == "province")] <- "pid"
  
  if("date_sick_biweek" %in% vars)
    vars <- c(vars, "date_sick_year")
  
  ## subset biweekly_data by the variables of interest, the observed counts, and the absolute errors
  AE_df <- biweekly_data[,c(vars, "obs_count", "AE_pred", "AE_lb", "AE_ub", "AE_seas", "AE_last_obs")]
  ## remove "AE_" from each column title
  colnames(AE_df) <- gsub("AE_", "", colnames(AE_df))

  if(length(vars) == 1){
    MAE_df <- AE_df %>%
      group_by_(vars) %>%
      summarise(n = n(),
                obs_count = mean(obs_count),
                pred = mean(pred),
                lb = mean(lb),
                ub = mean(ub),
                seas = mean(seas),
                last_obs = mean(last_obs)) %>%
      filter(n >= min_obs) %>%
      select_(vars, "obs_count", "n", "pred", "lb", "ub", denom)
  }

  if(length(vars) == 2){
    MAE_df <- AE_df %>%
      group_by_(vars[1], vars[2]) %>%
      summarise(n = n(),
                obs_count = mean(obs_count),
                pred = mean(pred),
                lb = mean(lb),
                ub = mean(ub),
                seas = mean(seas),
                last_obs = mean(last_obs)) %>%
      filter(n >= min_obs) %>%
      select_(vars[1], vars[2], "obs_count", "n", "pred", "lb", "ub", denom)
  }

  colnames(MAE_df)[dim(MAE_df)[2]] <- "denom"

  ## make new data frame with calculated rel_MAE
  rel_MAE_df <- MAE_df %>%
    mutate(rel_MAE = pred/denom,
           rel_MAE_lb = lb/denom,
           rel_MAE_ub = ub/denom)
  return(rel_MAE_df)
}

#' Calculate annual relative mean absolute error
#'
#' @param high_season_data a data frame with aggregated high season data, from aggregate_high_season_data()
#' @param denom the comparison metric, "median" compares to the median high season, "last" compares to last season, and "last_obs" uses the last observation before the forecast is made and projects it across the high season
#' @param steps_ahead the number of steps before the high season starts that the forecast is made
#'
#' @return a data frame with province, year, observed annual count, the mean absolute error for the predictions and comparison, and the relative mean absolute error.
calculate_annual_rel_MAE <- function(high_season_data,
                                     denom = c("median", "last_season", "last_obs"),
                                     steps_ahead) {
  if(denom == "median")
    denom <- "historic_median"
  if(denom == "last_obs")
    denom <- "last_obs_pred"
  if(denom == "last_season")
    denom <- "last_season_count"
  
  ## if steps_ahead is missing, then the week closest to the season start before the season starts
  if(missing(steps_ahead))
    steps_ahead <- which.min(high_season_data$biweeks_ahead[high_season_data$biweeks_ahead >= 0])

  colnames(high_season_data)[which(colnames(high_season_data) == denom)] <- "denominator"

  rel_MAE_df <- as.data.frame(high_season_data) %>%
    filter(biweeks_ahead %in% steps_ahead) %>%
    mutate(AE_pred = abs(pred_total_median - obs_annual_count),
           AE_denom = abs(denominator - obs_annual_count)) %>%
    group_by(pid, date_sick_year, obs_annual_count) %>%
    summarise(MAE_pred = mean(AE_pred),
              MAE_lb = mean(AE_lb),
              MAE_ub = mean(AE_ub),
              MAE_denom = mean(AE_denom)) %>%
    mutate(rel_MAE = MAE_pred / MAE_denom,
           rel_MAE_lb = MAE_lb / MAE_denom,
           rel_MAE_ub = MAE_ub / MAE_denom)

  return(rel_MAE_df)
}

#' Peak forecasts
#'
#' @param forecasts a data frame of forecasts, generated by run_one_year_forecasts()
#' @param high_season_data a data frame with aggregated high season data, from aggregate_high_season_data()
#' @param k number of weeks around the predicted peak week that the peak may fall into
#' @param steps_ahead the number of steps before the high season starts that the forecast is made
#'
#' @return a data frame with analysis date, province, year, observed peak, the percentage of simulations that correctly predicted the peak and the percentage of simulations that were within k weeks of the peak
peak_forecasts <- function(forecasts,
                           high_season_data,
                           k,
                           steps_ahead){
  peak_df <- left_join(forecasts, high_season_data) %>%
    filter(date_sick_biweek >= start_biweek,
           date_sick_biweek <= end_biweek,
           biweeks_ahead == steps_ahead) %>%
    group_by(pid, analysis_date, date_sick_year, obs_peak, simulation) %>%
    summarise(pred_peak = date_sick_biweek[which.max(pred_count)]) %>%
    mutate(correct_peak = ifelse(pred_peak == obs_peak, 1, 0),
           peak_wi_k = ifelse(abs(pred_peak - obs_peak) <= k, 1, 0)) %>%
    group_by(pid, analysis_date, date_sick_year, obs_peak) %>%
    summarise(peak_pct = mean(correct_peak),
              peak_wi_k_pct = mean(peak_wi_k))
  return(peak_df)
}


#' Forecast evaluation
#' A function that runs several specific functions to get an overview of how well the forecasts performed
#'
#' @param forecasts forecasts a data frame of formatted forecasts, generated by aggregate_biweekly_data()
#' @param counts a data frame of formatted counts, returned from format_counts()
#' @param season_start first biweek of high season
#' @param season_end last biweek of high season
#' @param to_date_lag the number of biweeks in the past to be discarded in order to have only full data
#' @param biweekly_data aggregated biweekly forecasts and counts, output by aggregate_biweekly_data()
#' @param seasons season(s) of interest
#' @param k number of weeks around the predicted peak week that the peak may fall into
#' @param trail_years number of years to use when determining historical counts in high season
#' @param pct_cases the minimum percentage of cases to be enclosed by the high season, for automatic selection
#' @param pred_outbreak_cutoff cutoff for outbreak_prob to make a binary outbreak prediction, default is 0.5
#'
#' @return
forecast_evaluation <- function(counts,
                             forecasts,
                             season_start,
                             season_end,
                             seasons,
                             trail_years,
                             pct_cases = NULL,
                             to_date_lag,
                             biweekly_data,
                             pred_outbreak_cutoff = 0.5,
                             k) {
  require(dengueThailand)
  
  data("thai_prov_data")
  ## set high season
  prov_high_season <- determine_high_season(counts = counts,
                                            season_start = season_start,
                                            season_end = season_end,
                                            analysis_year = seasons,
                                            trail_years = trail_years,
                                            pct_cases = pct_cases)

  ## find number of outbreaks in 2014
  number_of_outbreaks <- high_season_outbreaks(biweekly_data = biweekly_data,
                                               prov_high_season = prov_high_season,
                                               seasons = seasons,
                                               to_date_lag = to_date_lag)

  ## find sensitivity and specificity of outbreak predictions from 2014 for 0, 2, and 4 steps_out
  outbreak_preds <- eval_outbreak_preds(biweekly_data = biweekly_data,
                                        number_of_outbreaks = number_of_outbreaks,
                                        steps_out = c(0,2,4),
                                        pred_outbreak_cutoff = pred_outbreak_cutoff)

  ## prepare data to be evaluated at the seasonal (as opposed to biweekly) level
  high_season_data <- aggregate_high_season_data(forecasts = forecasts,
                                                 counts = counts,
                                                 prov_high_season = prov_high_season,
                                                 to_date_lag = to_date_lag,
                                                 analysis_year = seasons,
                                                 trail_years = trail_years,
                                                 season_start = season_start,
                                                 season_end = season_end,
                                                 pct_cases = pct_cases)

  pre_season_data <- filter(high_season_data, biweeks_ahead == 0)

  ## compare predicted counts with those of the biweekly medians
  rel_MAE_annual_median <- calculate_annual_rel_MAE(high_season_data = high_season_data,
                                                    denom = "median",
                                                    steps_ahead = 0)

  ## compare predicted counts with those of last season
  rel_MAE_annual_last_season <- calculate_annual_rel_MAE(high_season_data = high_season_data,
                                                  denom = "last_season",
                                                  steps_ahead = 0)

  ## compare predicted counts with an last_obs model (takes last full observation before delivery date and extrapolates it throughout high season)
  rel_MAE_annual_last_obs <- calculate_annual_rel_MAE(high_season_data = high_season_data,
                                                 denom = "last_obs",
                                                 steps_ahead = 0)

  ## compare peak predictions to peak observations
  peak_predictions <- peak_forecasts(forecasts = forecasts,
                                     high_season_data = high_season_data,
                                     k = k,
                                     steps_ahead = 0)

  summary_table <- data_frame(pid = number_of_outbreaks$pid,
                              pname = thai_prov_data$Province[match(pid, thai_prov_data$FIPS)],
                              num_outbreaks = number_of_outbreaks$num_outbreaks,
                              obs_annual_count = pre_season_data$obs_annual_count[match(pid, pre_season_data$pid)],
                              pred_annual_count = pre_season_data$pred_total_median[match(pid, pre_season_data$pid)],
                              rel_MAE_median = rel_MAE_annual_median$rel_MAE[match(pid, rel_MAE_annual_median$pid)],
                              rel_MAE_last_obs = rel_MAE_annual_last_obs$rel_MAE[match(pid, rel_MAE_annual_last_obs$pid)],
                              rel_MAE_last_season = rel_MAE_annual_last_season$rel_MAE[match(pid, rel_MAE_annual_last_season$pid)],
                              obs_peak = pre_season_data$obs_peak[match(pid, pre_season_data$pid)],
                              pred_peak_pct = peak_predictions$peak_pct[match(pid, peak_predictions$pid)],
                              pred_peak_wi_k_pct = peak_predictions$peak_wi_k_pct[match(pid, peak_predictions$pid)]
  )
  return(list(annual_pred_summary = summary_table,
              outbreak_preds = outbreak_preds,
              rel_MAE_annual_last_obs = rel_MAE_annual_last_obs,
              rel_MAE_annual_last_season = rel_MAE_annual_last_season,
              rel_MAE_annual_median = rel_MAE_annual_median,
              prov_high_season = prov_high_season,
              number_of_outbreaks = number_of_outbreaks,
              high_season_data = high_season_data))
}

