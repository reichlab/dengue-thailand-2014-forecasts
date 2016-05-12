

#' reformat a cntry.data.linelist object into a long-format dataset with delivery dates
#'
#' @param d a cntry.data.linelist object
#'
#' @return a tbl_df of biweekly case count data. 
#' Columns for province, date_sick_year, date_sick_biweek, delivery_biweek, 
#' and case_counts.
#' 
make_long_aggregated_counts_through_2014 <- function(d) {
    require(dengueThailand)
    require(spatialpred)
    require(dplyr)
    require(tidyr)
    require(lubridate)
    require(reshape2)
    
    data(thai_prov_data)
    thai_prov_narrow <- select(thai_prov_data, Province, province, FIPS) %>%
        dplyr::rename(prov_name = Province)
    
    ## put wide-form counts into long-form
    colnames(d) <- paste0("BW.", d@time.in.yr, ".", den_data@yr) 
    elongated_counts <- tbl_df(reshape2::melt(d@.Data, 
                                              varnames=c("prov_name", "char.time"))) %>%
        dplyr::rename(case_count = value) %>%
        separate(char.time, 
                 into = c("bw_label", "date_sick_biweek", "date_sick_year"),
                 convert=TRUE) %>%
        filter(date_sick_year < 2014) %>% ## 2014 cases to be added with deliv_dates
        select(-bw_label) %>%                       ## drop unneeded label column
        mutate(case_count = round(case_count)) %>%  ## round interpolated cases
        left_join(thai_prov_narrow) %>%
        mutate(prov_name = as.factor(prov_name)) 
    ## total records = 61256 = 76 provinces * 26 biweeks * 46 years
        
    ## aggregate 2014 linelist cases
    aggregated_linelist_2014 <- tbl_df(d@line.list) %>%
        filter(!is.na(date_sick),     ## remove NA date sick 
               year(date_sick)==2014, ## only 2014
               province != 38,        ## remove Bueng Kan, already added to Nong Khai 
               disease==26) %>%       ## subset to DHF
        mutate(date_sick_biweek = date_to_biweek(date_sick),
               date_sick_year = year(date_sick),
               province = as.numeric(province)) %>%
        group_by(province, delivery_date, date_sick_year, date_sick_biweek) %>%
        summarize(case_count = n()) %>%
        ungroup() %>%
        arrange(province, date_sick_year, date_sick_biweek, delivery_date) 
    
    ## for 2014, make table with one row per province-biweek-year-deliv_date
    ## total rows = 76 provinces * 26 biweeks * 29 delivery_dates * 1 year
    expanded_linelist_cols <- aggregated_linelist_2014 %>%
        filter(date_sick_year == 2014) %>%
        tidyr::expand(province, date_sick_year, date_sick_biweek, delivery_date)
        ## could remove the rows with delivery_date < date_sick
    
    ## merge expanded grid with observed data
    ## total rows = 148954
    full_aggregated_linelist <- left_join(expanded_linelist_cols, aggregated_linelist_2014) %>%
        left_join(thai_prov_narrow)
    
    ## merge with old data
    all_dhf_data <- bind_rows(full_aggregated_linelist, elongated_counts) %>%
        dplyr::rename(prov_num = province) %>%
        mutate(prov_name = factor(prov_name))
    
    return(all_dhf_data)
}



#' Plot all case data from an aggregated counts file
#'
#' @param counts results from a function like make_long_aggregated_counts_through_2014()
#' @param type scale of data, either "incidence" or "raw"
#'
#' @return NULL
make_all_data_raster <- function(counts, type=c("incidence", "raw")) {
    type <- match.arg(type)
    require(dengueThailand)
    require(dplyr)
    require(scales)
    require(grid)
    
    data(thai_prov_data)
    
    prov_data <- thai_prov_data %>% 
        transmute(FIPS=FIPS, pop = Population, moph_region = MOPH_Admin_Code)
    
    ## aggregate 2014 counts across delivery dates
     counts_aggregated <- counts %>% 
        ## for 2014 only, add sum counts with na.rm=TRUE 
        filter(date_sick_year == 2014) %>%
        group_by(date_sick_year, date_sick_biweek, FIPS, prov_name) %>%
        summarize(case_count = sum(case_count, na.rm=TRUE)) %>% 
        ungroup() %>%
        ## add back other years
        bind_rows(filter(counts, date_sick_year != 2014))
    
    ## merge couts and province data
    counts_with_prov <- left_join(counts_aggregated, prov_data) %>%
        mutate(prov_name  = reorder(prov_name, pop))
    
    ## rescale to avoid -Inf
    zero_idx <- which(counts_with_prov$case_count == 0)
    counts_with_prov[zero_idx, "case_count"] <- 1

    ## color schemes from 9-class OrRd: http://colorbrewer2.org/
    if(type=="incidence") {
        ggplot(counts_with_prov) + 
            geom_raster(aes(x=date_sick_year+(date_sick_biweek-1)/26, 
                            y=prov_name, 
                            fill=log10((case_count)/pop))) +
            #facet_grid(moph_region~., space="free", scales="free") +
            scale_fill_gradientn(colours=c("#fff7ec", "#fdbb84", "#7f0000"), #c("#ffffcc", "#fd8d3c", "#800026"),
                                 values = rescale(c(-7,-4,-2.2)),
                                 guide = "colorbar", limits=c(-7,-2.2),
                                 name="cases per 100,000 residents",
                                 na.value="lightgrey",
                                 labels=c(1, 10, 100, 500), breaks=c(-5, -4, -3, log10(5/1000)))+
            scale_x_continuous(expand=c(0,0)) +
            theme_bw() + xlab("") + ylab("") +
            theme(panel.margin = unit(0, "lines"))
            
    } else if (type=="raw") {
        ggplot(counts_with_prov) + 
            geom_raster(aes(x=date_sick_year+(date_sick_biweek-1)/26, 
                            y=prov_name, 
                            fill=log10(case_count))) + 
            #facet_grid(moph_region~., space="free", scales="free") +
            scale_fill_gradientn(colours=c("#fff7ec", "#fc8d59", "#7f0000"),
                                 values = rescale(c(0, 4)),
                                 guide = "colorbar", limits=c(0,4),
                                 name="number of cases",
                                 na.value="lightgrey",
                                 labels=c(0, 10, 100, 1000), breaks=log10(c(1, 11, 101, 1001))) +
            scale_x_continuous(expand=c(0,0)) +
            theme_bw() + xlab("") + ylab("") +
            theme(panel.margin = unit(0, "lines"))
    } 
}