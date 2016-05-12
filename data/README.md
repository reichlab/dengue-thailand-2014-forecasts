# Data codebook

field name           | description
---------------------|--------------
`pid`                | province ID
`analysis_date`      | date forecasts were generated
`date_sick_year`     | year of the date for which cases were reported
`date_sick_biweek`   | biweek of the date for which cases were reported
`epidemic_threshold` | epidemic threshold for the `date_sick` 
`obs_count`          | the observed case count for that `date_sick`
`pred_lb`            | the lower bound of the 95% prediction interval from the model
`pred_ub`            | the upper bound of the 95% prediction interval from the model
`pred_median`        | the median of the predictive distribution
`AE_lb`              | the absolute error of the lower bound
`AE_ub`              | the absolute error of the upper bound
`outbreak_prob`      | the outbreak probability
`seasonal_median`    | the seasonal median for the particular `date_sick`
`last_obs`           | the last observed count
`obs_outbreak`       | binary, whether incidence was above the epidemic threshold
`pred_covered`       | binary, whether the prediction interval covered the truth
`AE_pred`            | absolute error of the predicted median
`AE_seas`            | absolute error of the seasonal prediction
`AE_last_obs`        | absolute error of the last observation 
`step`               | the prediction horizon, i.e. the difference in biweeks between `date_sick` and `analysis_date`