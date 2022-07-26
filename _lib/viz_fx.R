# ----HEADER------------------------------------------------------------------------------------------------------------
# Author: JF
# Date: 07/26/2022
# Purpose: Store functions used for mapping and visualizations
# source("/homes/jfrostad/_code/ehd_mapsense/_lib/viz_fx.R", echo=T)
#***********************************************************************************************************************

# ----FUNCTIONS---------------------------------------------------------------------------------------------------------
#function to create custom maps 
cartographeR <- function(shapefile=tract_sf, dt,
                         map_varname, map_label=NA, map_title=NA,
                         subset_var=NA, subset_val=F,
                         facet_var=F,
                         filter_geocodes=drop_geocodes,
                         scale_type='cont',
                         scale_vals=NULL) {
  
  #if necessary, subset data
  if(subset_var %>% is.character) {
    dt <- dt[get(subset_var)==subset_val]
    
    #also cleanup the name for file
    subset_name <- str_replace_all(subset_val, ' ', "000") %>% 
      tolower %>% 
      str_replace_all("[^[:alnum:]]", "") %>% 
      str_replace_all('000', "_")
    
    message(subset_name)
    
  }
  
  #make variable to plot
  if(scale_type!='cont') dt[, map_var := get(map_varname) %>% as.factor]
  else dt[, map_var := get(map_varname)]
  
  #cleanup geocodes if needed
  if(filter_geocodes %>% is.character) dt <- dt[!(GEOID %in% filter_geocodes)]
  
  #merge dropouts to shapefile and plot
  shp <- shapefile %>% 
    merge(dt, by='GEOID', allow.cartesian=T) 
  
  #make plot
  plot <- ggplot() + 
    geom_sf(data = shp, aes(fill = map_var), lwd=0) + 
    geom_sf(data = water_sf, lwd=0, color = 'gray70', fill = 'gray95') +
    #TODO seems to make the map tilted even though this is the NAD83 proj?
    #coord_sf(crs=4269) + 
    theme_void()
  
  #facet by theme if needed
  if(facet_var %>% is.character) plot <- plot + facet_wrap(reformulate(facet_var))
  
  #add on the colorscale
  #TODO could probably just make this reflexive based on input data
  if(scale_type=='div') plot <- plot + scale_fill_brewer(map_label, palette='RdBu', na.value = "grey75", direction=-1)
  else if(scale_type=='div_man') plot <- plot + scale_fill_manual(map_label, values=scale_vals, na.value = "grey75")
  else if(scale_type=='drops') plot <- plot + scale_fill_brewer(map_label, palette='PuOr', na.value = "grey75", direction=-1)
  else if(scale_type=='bin') plot <- plot + scale_fill_viridis_d(map_label, option = "plasma", na.value = "grey75")
  else if(scale_type=='cont') plot <- plot + scale_fill_viridis_c(map_label, option='magma', na.value = "grey75")
  
  #title 
  #TODO add more functionality
  if(map_title %>% is.character & subset_var %>% is.na) plot <- plot + ggtitle(map_title)
  else if(subset_var %>% is.character  & map_title %>% is.na) plot <- plot + ggtitle(subset_name %>% to_upper_camel_case)
  else if(subset_var %>% is.character & map_title %>% is.character) plot <- plot + 
    ggtitle(paste0(map_title, ': \n', subset_name %>% to_upper_camel_case))
  
  #save the plot
  file.path(viz.dir, paste0(map_varname, '_',
                            ifelse(subset_var %>% is.character,
                                   subset_name,
                                   ''),
                            ifelse(facet_var %>% is.character,
                                   paste0('_by_', facet_var),
                                   ''),
                            '_map.png')
  ) %>% ggsave(height=8, width=12)
  
}

#***********************************************************************************************************************
#
# ---SCRAP -------------------------------------------------------------------------------------------------------------

#***********************************************************************************************************************