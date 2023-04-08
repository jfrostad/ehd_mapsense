# ----HEADER------------------------------------------------------------------------------------------------------------
# Author: JF
# Date: 07/26/2022
# Purpose: Store functions used for mapping and visualizations
# source("/homes/jfrostad/_code/ehd_mapsense/_lib/viz_fx.R", echo=T)
#***********************************************************************************************************************

# ----FUNCTIONS---------------------------------------------------------------------------------------------------------
#function to create custom maps 
cartographeR <- function(shapefile=tract_sf, dt,
                         lvl=1, #default is overall rank (2=theme, 3=measure)
                         map_varname, map_label=NA, map_title=NA, tag=NA, #use tag to label files
                         subset_var=NA, subset_val=F,
                         facet_var=F,
                         filter_geocodes=drop_geocodes,
                         scale_type='cont',
                         scale_vals=NULL,
                         get_plot=F) {

  #filter long dataset to appropriate level
  dt <- dt[level==lvl]

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
  if(scale_type%like%'cont'&!(scale_type%like%'man')) dt[, map_var := get(map_varname)] 
  else dt[, map_var := get(map_varname) %>% as.factor]
  
  #cleanup missing data
  #remove the -1s which represent nonmissing NAs (plot as NA)
  #but note that there can be negatives for the change vars
  if(map_varname %like% 'measure') dt <- dt[map_var>=0]
  
  #cap variables for manual scales
  # if(scale_type%like%'man') {
  #   
  #   browser()
  #   
  #   message('manual scale specified, capping values at extremes')
  #   scale_extremes <-
  #     scale_vals %>%
  #     names %>% 
  #     as.numeric %>% 
  #     summary %>% 
  #     as.matrix %>% 
  #     .[c(1,6)]
  #   
  #   #cap max
  #   message(max(dt$map_var %>% as.numeric, na.rm=T), ' capped to ', scale_extremes[2])
  #   dt[map_var > scale_extremes[2], map_var := scale_extremes[2]]
  #   
  #   #cap min
  #   message(min(dt$map_var, na.rm=T), ' capped to ', scale_extremes[1])
  #   dt[map_var < scale_extremes[1], map_var := scale_extremes[1]]
  #   
  # }
  
  #cleanup geocodes if needed
  if(filter_geocodes %>% is.character) dt <- dt[!(GEOID %in% filter_geocodes)]

  #merge dropouts to shapefile and plot
  shp <- shapefile %>% 
    merge(dt, by='GEOID', allow.cartesian=T) 

  #make plot
  plot <- ggplot() + 
    geom_sf(data = shp, aes(fill = map_var), lwd=0) + 
    geom_sf(data = road_sf, lwd=.25, color = 'gray20', fill = 'gray95', show.legend = 'line') +
    geom_sf(data = water_sf, lwd=0, color = 'gray70', fill = 'gray95') +
    #TODO seems to make the map tilted even though this is the NAD83 proj?
    #coord_sf(crs=4269) + 
    theme_void()
  
  #facet by theme if needed
  if(facet_var %>% is.character) plot <- plot + facet_wrap(reformulate(facet_var))

  #add on the colorscale
  #TODO could probably just make this reflexive based on input data
  if(scale_type=='div') plot <- plot + scale_fill_distiller(map_label, palette='RdBu', na.value = "grey75", direction=-1)
  else if(scale_type=='div_man') plot <- plot + scale_fill_manual(map_label, values=scale_vals, na.value = "grey75")
  else if(scale_type=='drops') plot <- plot + scale_fill_brewer(map_label, palette='PuOr', na.value = "grey75", direction=-1)
  else if(scale_type=='bin') plot <- plot + scale_fill_viridis_d(map_label, option = "plasma", na.value = "grey75")
  else if(scale_type=='cont_man') plot <- plot + scale_fill_manual(map_label, values=scale_vals, na.value = "grey75")
  else if(scale_type=='cont') plot <- plot + scale_fill_viridis_c(map_label, option='magma', na.value = "grey75")
  else if(scale_type=='cont_grad') plot <- plot + scale_fill_gradient2(map_label, na.value = "grey75")
  else if(scale_type=='cont_vir') plot <- plot + scale_fill_viridis_c(map_label, option='viridis', na.value = "grey75")
  else if(scale_type=='bivar') plot <- plot + bi_scale_fill(pal='GrPink', guide='none') + bi_theme()
  else if(scale_type=='identity') { plot <- plot + scale_fill_identity(map_label, labels = scale_vals %>% names, 
                                                                       scale_vals,
                                                                       guide='legend') }
  
  
  #title 
  #TODO add more functionality
  if(map_title %>% is.character & subset_var %>% is.na) plot <- plot + ggtitle(map_title)
  else if(subset_var %>% is.character  & map_title %>% is.na) plot <- plot + ggtitle(subset_name %>% to_upper_camel_case)
  else if(subset_var %>% is.character & map_title %>% is.character) plot <- plot + 
    ggtitle(paste0(map_title, ': \n', subset_name %>% to_upper_camel_case))
  
  #save the plot
  #TODO make it so you can change the output dir?
  file.path(viz.dir, paste0(map_varname, '_',
                            ifelse(tag %>% is.character,
                                   tag,
                                   ''),
                            ifelse(subset_var %>% is.character,
                                   subset_name,
                                   ''),
                            ifelse(facet_var %>% is.character,
                                   paste0('_by_', facet_var),
                                   ''),
                            '_map.png')
  ) %>% ggsave(height=8, width=12)
  
  if(get_plot) return(plot)
  
}

#***********************************************************************************************************************
#
# ---SCRAP -------------------------------------------------------------------------------------------------------------

#***********************************************************************************************************************