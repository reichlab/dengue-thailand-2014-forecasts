## function to make outbreak probability plot for paper
## Nicholas Reich 
## June 2015

## adapted from dengueThailand function new_plot_forecast_map.R

#'@param forecast_data forecast file dataframe
#'@param cdata cntry.data object with spatial polygons needed for plotting
#'@param biweek_to_plot biweek to plot
#'@param analysis_date analysis date to plot predictions from
#'@param include_legend logical, whether to include legend
#'@param plot_type one of either "incidence" or "outbreak"

plot_outbreak_map <- function(forecast_data,
                              wld.shape = "../../data/updated_adm/ne_10m_admin_1_states_provinces.shp",
                              biweek_to_plot, 
                              analysis_date,
                              include_legend=TRUE,
                              plot_type=c("incidence", "outbreak")) {
    require(ggplot2)
    require(dplyr)
    require(rgeos)
    require(mapproj)
    require(maptools)
    
    
    data(thai_prov_data)
    
    ## load in spatial polygons
    world <- readShapeSpatial(wld.shape)
    thai_poly <- world[world$admin=="Thailand",]
    
    ## adjust for end of year biweeks
    if(biweek_to_plot>26)
        biweek_to_plot <- biweek_to_plot - 26
    
    if(!(biweek_to_plot %in% unique(forecast_data$date_sick_biweek)))
        stop("date_sick_biweek must be in forecast_data.")
    
    if(nrow(filter(forecast_data, 
                   analysis_date == analysis_date, 
                   date_sick_biweek=biweek_to_plot)) < 1)
        stop("no match for plotting biweek and analysis date.")
    
    ## merge thai_prov_data with forecasts to get population
    forecast_data_merged <- left_join(forecast_data, thai_prov_data, by = c("pid" = "FIPS")) %>%
        mutate(incidence = pred_median/Population)
    
    ## foritfy polygon info
    thai_poly <- thai_poly[which(thai_poly@data$gns_adm1 != "TH81"),]
    thai_poly$ID_1 <- dense_rank(thai_poly$woe_name)
    thai_locs <- fortify(thai_poly, region="ID_1")
    thai_locs[['region']] <- thai_locs[['id']]
    
    ## store loc.info@data
    loc_info <- thai_poly@data
    
    ## combine polygon info with thai data
    data_to_plot <- left_join(loc_info, thai_prov_data, by = c("gns_adm1" = "FIPS")) %>%
        mutate(id = ID_1,
               pid = as.character(gns_adm1))
    
    
    ## plotting choices based on type
    if(plot_type=="incidence") {
        fill_var <- "log10(incidence)"
        plot_lims <- range(forecast_data_merged$incidence)
        plot_lims <- c(floor(log10(plot_lims)[1]),
                       ceiling(log10(plot_lims)[2]))
        plot_breaks <- seq(plot_lims[1], plot_lims[2])
        plot_midpoint <- mean(plot_lims)
        legend_title <- "incidence"
        plot_labels <- paste0("1e", plot_breaks)
        
    } else {
        fill_var <- "outbreak_prob"
        plot_lims <- c(0,1)
        plot_midpoint <- .5
        plot_breaks <- c(0, .5, 1)
        plot_labels <- c(0, .5, 1)
        legend_title <- "outbreak probability"
    }
    
    ## set legend position, if any
    legend_pos <- ifelse(include_legend, "right", "none")
    
    ## text for map label
    forecast_data_subset <- subset(forecast_data_merged, 
                                   date_sick_biweek == biweek_to_plot,
                                   analysis_date == analysis_date)
    map_date <- format(as.Date(biweek_to_date(biweek_to_plot, 
                                              forecast_data_subset$date_sick_year[1])), "%d %b %Y")
    
    ## merge forecast data with data_to_plot
    new_dp <- left_join(data_to_plot, forecast_data_subset)
    
    sp_map <- ggplot(new_dp, aes(map_id=id)) + 
        geom_map(aes_string(fill=fill_var), map=thai_locs) + 
        expand_limits(x = thai_locs$long, y = thai_locs$lat) +
        ## use color blind friendly colors
        scale_fill_gradient2(low = "#053061", mid="#CCCCCC", high = "#FF2C19", 
                             name=legend_title,
                             limits=plot_lims, 
                             midpoint=plot_midpoint, 
                             breaks=plot_breaks,
                             labels=plot_labels) +
        theme_bw() +
        theme(axis.ticks = element_blank(), 
              axis.text = element_blank(),
              panel.background = element_rect(fill = "transparent",colour = NA),
              legend.position = legend_pos) + # or element_blank()
        coord_map(projection="mercator") + ## to keep scaling right
        annotate("text", x = 103.5, y = 20.1, label = map_date, size=4) +
        xlab("") + ylab("")        
    return(sp_map)
}