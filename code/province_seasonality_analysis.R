## find amplitude of seasonal time series

prov_seasonalities <- function(data, plot=FALSE) {  
    require(mgcv)
    require(dplyr)
    require(dengueThailand)
    nprov <- length(unique(data$prov_num))
    max_seas <- data_frame(prov_num = unique(data$prov_num), seasonality=NA, resid_sd_seas=NA)
    for (i in 1:nprov) {
        prov_dat <- filter(data, prov_num == unique(data$prov_num)[i], date_sick_year >= 1995) %>%
            mutate(date_sick = date_sick_year + (date_sick_biweek-1)/26)
        ## fit seasonal model with secular trend data
        fm1 <- gam(case_count~s(date_sick, k=3) + s(date_sick_biweek, bs="cc"), data=prov_dat, family="poisson")
        fm1_output <- plot(fm1, select=2)
        max_seas[i, "seasonality"] <- max(fm1_output[[2]]$fit)
        max_seas[i, "resid_sd_seas"] <- sd(fm1$residuals)
    }
    return(max_seas)
}

library(dplyr)
load("~/Documents/code_versioned/spamd/trunk/source/realtime-models/20150806-spamd-nomem-allyrs-10step-6lag-casedata.RData")
dat <- dhf_data %>% 
    group_by(prov_num, FIPS, prov_name, date_sick_year, date_sick_biweek) %>%
    summarize(case_count = sum(case_count, na.rm=TRUE))

prov_seasonalities <- prov_seasonalities(dat)

save(prov_seasonalities, file="data/20151027-province-seasonalities.RData")
