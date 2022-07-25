# ----HEADER------------------------------------------------------------------------------------------------------------
# Author: JF
# Date: 04/08/2020
# Purpose: Produce HAP paper results
# source("/homes/jfrostad/_code/lbd/hap/cooking/rover/opportunity.R", echo=T)
#***********************************************************************************************************************

# ----CONFIG------------------------------------------------------------------------------------------------------------
# clear memory
rm(list=ls())

# runtime configuration
if (Sys.info()["sysname"] == "Linux") {
  j_root <- "/home/j/"
  h_root <- "/homes/jfrostad/"
  arg <- commandArgs()[-(1:3)] # First args are for unix use only
  
  if (length(arg)==0) {
    # arg <- c("IND", #current project iteration
    #          "8", #output version
    #          1) #number of cores provided to multicore functions
  }
  
} else {
  j_root <- "J:"
  h_root <- "H:"
  # arg <- c("IND", #current project iteration
  #          "4", #output version
  #          1) #number of cores provided to multicore functions
}

#use cairo to render instead of quartz (quartz causes big slowdowns with geom_sf)
if(!identical(getOption("bitmapType"), "cairo") && isTRUE(capabilities()[["cairo"]])){
  options(bitmapType = "cairo")
}

#set opts
options(scipen=999) #readability

## Set core_repo location and indicator group
user            <- Sys.info()['user']
my_repo         <- '/homes/jfrostad/_code/ehd_mapping'

#load packages
package_lib    <- sprintf('%s_code/_lib/pkg_R',h_root)
## Load libraries and  MBG project functions.
.libPaths(package_lib)

#--------------------------------------------------------
# Purpose: re-create IBL ranks & QA observed vs expected ranks
# Author: Christopher Ahmed (christopher.ahmed@doh.wa.gov)
# Date: January 2022
#
# Note: this script reads in a data extract created from the backend
#   of WTN. The extract includes all the data IBL uses to create its ranks
#   the extract is a long data-set as opposed to wide.
#--------------------------------------------------------

pacman::p_load(tidyverse, readxl, snakecase, janitor, data.table, naniar, visdat,
               magrittr, scales, ggplot2, ggpubr, ggridges, ggrepel, gridExtra, isoband, RColorBrewer, 
               sf, viridis, farver, reldist, ggnewscale, ggallin,
               daff, waldo, tigris, tidycensus, ggcorrplot,
               stargazer)

#***********************************************************************************************************************

# ----IN/OUT------------------------------------------------------------------------------------------------------------
###Input###
#raw data
main.dir <- file.path('/mnt/share/scratch/users/jfrostad/ehd_mapping/')
data.dir <- file.path(main.dir, 'data')
local.dir <- file.path('C:/Users/jfrostad/Downloads/Datafiles_EHDMapv2.0_AllCTs/Data All CTs')
data.dir <- file.path(main.dir, 'data')
data_extract_EHDv1 <- 'ehd_data_v1.xlsx' #TODO rename
data_extract_EHDv2 <- 'ehd_data_v3.xlsx' #TODO rename

###Output###
out.dir <- file.path('C:/Users/jfrostad/Desktop/shiny_app/data')
viz.dir  <- file.path(local.dir, 'viz')
#***********************************************************************************************************************

# ---FUNCTIONS----------------------------------------------------------------------------------------------------------
#test space
cleanAndRank <- function(dir, path,
                         nranks=10, 
                         clean_names=F,
                         debug=F) {
  
  if(debug) browser()
  
  # File Details:
  # Weights: Measures: All 1's = No Special Weights
  #          Themes: "Environmental Effects" = .5, Other Themes Are 1's = No Special Weights
  # Rank Method: Themes: Average
  #              Index: (("Environmental Effects"+"Environmental Exposures")/2)*
  #                     (("Socioeconomic Factors"+"Sensitive Populations")/2)
  # Read Data Extract (3-Sheets)
  message('reading')
  extract_measure <- read_xlsx(file.path(dir,path), sheet = "Measure",
                               col_types=c('numeric', 'numeric', 'numeric',
                                           'guess', 'guess', 'guess',
                                           'guess', 'guess', 'guess', 
                                           'guess')) %>% as.data.table
  
  extract_theme <- read_xlsx(file.path(dir,path), sheet = "Theme",
                             col_types=c('numeric', 'numeric', 'numeric',
                                         'guess', 'guess', 'numeric',
                                         'guess', 'guess')) %>% as.data.table
  
  extract_index <- read_xlsx(file.path(dir,path), sheet = "Index",
                             col_types=c('numeric', 'numeric', 'numeric',
                                         'guess', 'guess', 'guess',
                                         'guess', 'guess')) %>% as.data.table
  
  # Clean Data Extract Data
  # Apply Weights If Not = 1
  message('cleaning')
  measure_dt <- extract_measure %>% 
    setnames(., names(.), names(.) %>% tolower) %>% 
    .[rank==-1, rank := NA] %>% 
    .[rankcalculatedvalue==-1, rankcalculatedvalue := NA] %>% 
    .[, .(GEOID=as.character(geocode), item=itemname, theme=themename, measure_rank=rank, 
          measure_rank_val=rankcalculatedvalue)]
  
  index_dt <- extract_index %>% 
    setnames(., names(.), names(.) %>% tolower) %>% 
    .[, .(GEOID=as.character(geocode), index_rank=rank, index_rank_val=rankcalculatedvalue)]
  
  theme_dt <- extract_theme %>% 
    setnames(., names(.), names(.) %>% tolower) %>% 
    .[, theme_weights := 1] %>% 
    .[themename=="Environmental Effects", theme_weights := .5] %>% 
    .[, .(GEOID=as.character(geocode), theme=themename, theme_weights,
          theme_rank_val=rankcalculatedvalue, theme_rank=rank)] %>% 
    merge(index_dt, by='GEOID') #merge on index information
  
  # Remap item names
  if(clean_names %>% is.data.table) {
    
    measure_dt <- setnames(measure_dt, 'item', 'item_old') %>% 
      merge(item_map, by='item_old')
    
  }
  
  # Calculate Ranks: Data Measures
  message('calculating measure ranks')
  measure_ranks <- measure_dt[, item := factor(item)] %>% 
    .[, rank_order := frank(measure_rank_val, 
                            na.last='keep', 
                            ties.method = 'min'), by=item] %>% 
    .[, bin_size := floor(sum(!is.na(measure_rank_val))/nranks), by=item] %>% 
    #bin the ranks and number them by rounding up
    .[, measure_rank_integer := (rank_order / bin_size) %>% ceiling, by=item] %>% 
    # Non-Missing Tract's WIth NA's Are Given A Measure Rank = -1
    .[measure_rank_integer>nranks, measure_rank_integer := nranks] %>% 
    setnafill(type='const', fill=-1, cols=c('measure_rank_val',
                                            'measure_rank_integer',
                                            'measure_rank')) %>% 
    .[, qa_measures := measure_rank == measure_rank_integer] #verify that calcs match published
  
  
  # Calculate Ranks: Themes
  #crossjoin to create all combos of geocode and theme+name
  #thene fill with 0s
  message('calculating theme ranks')
  theme_ranks <- measure_dt[, CJ(GEOID, item, unique=T), by=theme] %>% 
    merge(measure_dt, by=c('GEOID', 'item', 'theme'), all.x=T) %>% 
    merge(theme_dt, by=c('GEOID', 'theme')) %>% 
    # Missing Tracts Are Added To Measures With A Measure Rank = 0
    setnafill(., type='const', fill=0, cols=names(.) %>% .[.%like%'rank']) %>% 
    .[, measure_num := .N, by=.(GEOID, theme)] %>% 
    .[, measure_rank_sum := sum(measure_rank_integer), by=.(GEOID, theme)] %>% 
    .[, measure_rank_sum := sum(measure_rank_integer), by=.(GEOID, theme)] %>% 
    unique(., by=c('GEOID', 'theme')) %>% 
    .[, item := ''] %>%  #item no longer relevant after collapse
    #TODO ranking could be better as a function??
    .[, theme_avg_rank := round_half_up(measure_rank_sum/measure_num, 2), by=theme] %>% 
    .[, theme_rank_order := frank(theme_avg_rank, 
                                  na.last='keep', 
                                  ties.method = 'min'), by=theme] %>% 
    .[, theme_bin_size := floor(sum(!is.na(theme_avg_rank))/nranks), by=theme] %>% 
    #bin the ranks and number them by rounding up
    .[, theme_rank_integer := (theme_rank_order / theme_bin_size) %>% ceiling, by=theme] %>% 
    # Non-Missing Tract's WIth NA's Are Given A Measure Rank = -1
    .[theme_rank_integer>nranks, theme_rank_integer := nranks] %>% 
    .[, qa_themes := theme_rank == theme_rank_integer] #verify that calcs match published
  
  
  # Calculate Ranks: Index
  # File Details Above Indicates Theme "Environmental Effects" is Weighted At 0.5 And
  # The Rank Method = (("Environmental Effects"+"Environmental Exposures")/2)*
  #                     (("Socioeconomic Factors"+"Sensitive Populations")/2)
  # The First Sum In This Method = The Pollution Burden, The Second = Population Characteristics
  # The Theme Field We Start With = calc_theme_avg_measure_ranks
  message('calculating index ranks')
  index_ranks <- setorder(theme_ranks, GEOID, theme) %>% 
    .[, weighted_ranks := theme_avg_rank * theme_weights] %>% 
    .[theme%like%'Environmental', burden_sum := sum(weighted_ranks), by=GEOID] %>% 
    setnafill(type='locf', cols='burden_sum') %>%  #env ordered first, so carry fwd
    .[theme%like%'Environmental', burden_avg := mean(weighted_ranks), by=GEOID] %>% 
    setnafill(type='locf', cols='burden_avg') %>%  #env ordered first, so carry fwd
    .[!(theme%like%'Environmental'), pop_chars := sum(weighted_ranks), by=GEOID] %>% 
    setnafill(type='nocb', cols='pop_chars') %>%  #env ordered first, so carry back
    .[!(theme%like%'Environmental'), pop_avg := mean(weighted_ranks), by=GEOID] %>% 
    setnafill(type='nocb', cols='pop_avg') %>%  #env ordered first, so carry back
    unique(., by=c('GEOID')) %>% 
    .[, item := ''] %>%  #theme and item no longer relevant after collapse
    .[, theme := ''] %>%  #theme and item no longer relevant after collapse
    #scale both per cal enviroscreen standard
    .[, burden_scaled := burden_avg / max(burden_avg) * 10] %>% 
    .[, pop_scaled := pop_avg / max(pop_avg) * 10] %>% 
    #TODO ranking could be better as a function??
    .[, index_avg_rank := round_half_up(burden_sum/2*pop_chars/2, 
                                        2)] %>% 
    .[, index_rank_cal := burden_avg * pop_avg] %>% 
    #.[order(index_avg_rank), index_rank_int := floor( 1 + nranks * (.I-1)/.N)] %>% 
    .[, index_rank_order := frank(index_avg_rank, 
                                  na.last='keep', 
                                  ties.method = 'min')] %>% 
    .[, index_bin_size := floor(sum(!is.na(index_avg_rank))/nranks)] %>% 
    #bin the ranks and number them by rounding up
    .[, index_rank_integer := (index_rank_order / index_bin_size) %>% ceiling] %>% 
    .[index_rank_integer>nranks, index_rank_integer := nranks] %>% 
    #redo everything with the scaled cal method
    .[, index_rank_order := frank(index_rank_cal, 
                                  na.last='keep', 
                                  ties.method = 'min')] %>% 
    .[, index_bin_size := floor(sum(!is.na(index_rank_cal))/nranks)] %>% 
    #bin the ranks and number them by rounding up
    .[, index_rank_integer_cal := (index_rank_order / index_bin_size) %>% ceiling] %>% 
    # Non-Missing Tract's WIth NA's Are Given A Measure Rank = -1
    .[index_rank_integer>nranks, index_rank_integer := nranks] %>% 
    .[index_rank_integer_cal>nranks, index_rank_integer_cal := nranks] %>% 
    .[, qa_index := index_rank == index_rank_integer] #verify that calcs match published
  
  list('measure'=measure_ranks,
       'theme'=theme_ranks,
       'index'=index_ranks,
       'measure_raw'=measure_dt) %>% 
    return
  
}

#helper function to copy things out of R
writeExcel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
}

#label outliers statistically
is_outlier <- function(x) {
  return(x < quantile(x, 0.25, na.rm=T) - 1.5 * IQR(x, na.rm=T) | x > quantile(x, 0.75, na.rm=T) + 1.5 * IQR(x, na.rm=T))
}

#***********************************************************************************************************************

# ---PREP DATA----------------------------------------------------------------------------------------------------------
##read in and prep datasets for analysis##
#names of themes have changed, make a map
item_map <- file.path(local.dir, 'ehd_map_theme_names.csv') %>% fread

#life expectancy data
le_dt <-  file.path(local.dir, 'le_at_birth_2015_2019.csv') %>% fread
setnames(le_dt, names(le_dt), 
         c('county', 'geocode', 'le', 'lower', 'upper'))
le_dt[, c('le', 'lower', 'upper') := lapply(.SD, as.numeric), .SDcols=c('le', 'lower', 'upper')]

#also bring in the census tracts shapefile in order to do some cartography
#can be downloaded from the census website using tigris
tract_sf <- tracts('WA', cb=T) #%>% 
  #st_transform(32148) #%>% 
  #erase_water(area_threshold = 0.9) #intersect with water overlay and remove

# tract_sf <- file.path(local.dir, 'shapefile', 'tl_2010_53_tract10.shp') %>% 
#   st_read

#use the water shapefile as an overlay
counties_list <- counties('WA', cb=T)
water_sf <- area_water('WA', counties_list$COUNTYFP %>% unique) %>% 
  st_simplify(preserveTopology = TRUE, dTolerance = 100)

#also overlay roads
road_sf <- primary_secondary_roads(state='WA') %>% 
  st_simplify(preserveTopology = TRUE, dTolerance = 100)

#also overlay places
  places_sf <- places(state = 'WA', cb = T) 

#first read in and calculate all the ranks
ranks_old <- cleanAndRank(dir=local.dir, path=data_extract_EHDv1, clean_names=item_map, debug=F)
ranks_new <- cleanAndRank(dir=local.dir, path=data_extract_EHDv2)

#merge measures (old v. new) to compare
measure_ranks <- merge(ranks_old$measure[, .(GEOID, item, theme, 
                                          measure_old=measure_rank_integer)],
                    ranks_new$measure[, .(GEOID, item, theme, 
                                          measure_new=measure_rank_integer),],
                    by=c('GEOID', 'item', 'theme'), all = T) %>% 
  .[, measure_shift := measure_new-measure_old] %>% 
  #.[, avg_shift := mean(measure_shift, na.rm=T), by='item'] %>% 
  #.[, measure_old_neg := measure_old * -1] %>% 
  .[, measure_shift_capped := measure_shift] %>% 
  .[measure_shift>=5, measure_shift_capped := 5] %>% 
  .[measure_shift<=-5, measure_shift_capped := -5]  %>% 
  .[, measure_shift := measure_shift_capped] %>% 
  .[, measure_shift_capped := NULL]
  

#merge raw measures (old v. new) to compare
measure_raw <- merge(ranks_old$measure_raw[, .(GEOID, item, theme, 
                                             measure_old=measure_rank_val)],
                       ranks_new$measure[, .(GEOID, item, theme, 
                                             measure_new=measure_rank_val)],
                       by=c('GEOID', 'item', 'theme'), all=T) %>% 
  .[, measure_shift_raw := measure_new-measure_old] %>% 
  #.[, avg_shift := mean(measure_shift_raw, na.rm=T), by='item'] %>% 
  #.[, measure_old_neg := measure_old * -1] %>% 
  setnames(., c('measure_new', 'measure_old'), c('measure_new_raw', 'measure_old_raw'))

#merge both measures datasets
measure_dt <- merge(measure_ranks,
                    measure_raw[, .(GEOID, item, measure_old_raw, measure_new_raw, measure_shift_raw)],
                    by=c('GEOID', 'item'))

#merge themes
theme_dt <- merge(ranks_old$theme[, .(GEOID, item, theme, 
                                      index_old=theme_rank_integer)],
                  ranks_new$theme[, .(GEOID, item, theme, 
                                      index_new=theme_rank_integer)],
                  by=c('GEOID', 'item', 'theme')) %>% 
  .[, index_shift := index_new-index_old] %>% 
  .[, index_shift_capped := index_shift] %>% 
  .[index_shift>=5, index_shift_capped := 5] %>% 
  .[index_shift<=-5, index_shift_capped := -5]  %>% 
  .[, index_shift := index_shift_capped] %>% 
  .[, index_shift_capped := NULL]

#merge indexes (old v. new) to compare
index_dt <- merge(ranks_old$index[, .(GEOID, item, theme, 
                                          index_old=index_rank_integer)],
                    ranks_new$index[, .(GEOID, item, theme, 
                                          index_new=index_rank_integer)],
                    by=c('GEOID', 'item', 'theme')) %>% 
  .[, index_shift := index_new-index_old] %>% 
  .[, index_shift_capped := index_shift] %>% 
  .[index_shift>=5, index_shift_capped := 5] %>% 
  .[index_shift<=-5, index_shift_capped := -5] %>% 
  .[, index_shift := index_shift_capped] %>% 
  .[, index_shift_capped := NULL]

#identify dropouts
threshold_val <- 9
index_dt[, impacted_new := 0]
index_dt[, impacted_old := 0]
index_dt[index_new>=threshold_val, impacted_new := 1]
index_dt[index_old>=threshold_val, impacted_old := 1]
index_dt[, dropout := as.factor(impacted_old - impacted_new)]

#add life expectancy data
index_dt <- merge(index_dt, le_dt, by.x='GEOID', by.y='geocode')
index_dt[, le_state_average := mean(le, na.rm=T)]
index_dt[, life_expectancy := paste0(le, ' years (', lower, '-', upper, ')')]

#merge the dropouts back onto the other tables for graphing
measure_dt <- merge(measure_dt, index_dt[, .(GEOID, dropout, index=index_new, life_expectancy)], by='GEOID')
theme_dt <- merge(theme_dt, index_dt[, .(GEOID, dropout, index=index_new, life_expectancy)], by='GEOID')

#drop water codes with weird data from maps
drop_geocodes <- c('53057990100' #san juan water area with lots of big changes
                   )

#label the outliers
#index_dt[, outlier_lab := ifelse(is_outlier(le), name, NA_character_), by=index_new]

#save all data
out <- list(
  'index'=index_dt,
  'ranks'=measure_ranks,
  'ranks_new'=ranks_new,
  'ranks_old'=ranks_old,
  'measures'=measure_raw,
  'themes'=theme_dt,
  'tracts'=tract_sf,
  'water'=water_sf
)
saveRDS(out, file=file.path(out.dir, 'all_data.RDS'))

#save a lite version of the data for the online mapping tool
out <- list(
  'index'=index_dt[, -c('lower', 'upper', 'le_state_average'), with=F],
  'measures'=measure_dt,
  'themes'=theme_dt,
  'tracts'=tract_sf,
  'water'=water_sf,
  'roads'=road_sf,
  'places'=places_sf
)
saveRDS(out, file=file.path(out.dir, 'viz_data.RDS'))
#***********************************************************************************************************************

# ---MAP----------------------------------------------------------------------------------------------------------------

#plot dropouts
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


#create a manual diverging color scale to make sure that the index shifts are uniformly depicted
div_colors <- RColorBrewer::brewer.pal(11, 'RdBu') %>% rev
names(div_colors) <- unique(measure_ranks$measure_shift_capped) %>% sort

#create a manual discrete color scale for the continuous index and make sure the color scale legend has all integers
#cont_colors <- viridis::magma(n = 10)
#update colors with allies scale
cont_colors <- c(
  rgb(248, 250, 251, maxColorValue = 255), # 1:     Red: 248, Green: 250, Blue: 251
  rgb(236, 242, 246, maxColorValue = 255), # 2:     Red: 236, Green: 242, Blue: 246
  rgb(220, 230, 239, maxColorValue = 255), # 3:     Red: 220, Green: 230, Blue: 239
  rgb(203, 218, 233, maxColorValue = 255), # 4:     Red: 203, Green: 218, Blue: 233
  rgb(194, 199, 223, maxColorValue = 255), # 5:     Red: 194, Green: 199, Blue: 223
  rgb(194, 178, 213, maxColorValue = 255), # 6:     Red: 194, Green: 178, Blue: 213
  rgb(192, 157, 203, maxColorValue = 255), # 7:     Red: 192, Green: 157, Blue: 203
  rgb(189, 132, 186, maxColorValue = 255), # 8:     Red: 189, Green: 132, Blue: 186
  rgb(161, 124, 177, maxColorValue = 255), # 9:     Red: 161, Green: 124, Blue: 177
  rgb(163, 124, 162, maxColorValue = 255) # 10:   Red: 163, Green: 124, Blue: 162
)

names(cont_colors) <- 1:10

#create a series of colored maps for the different vars in the dataset
  #denote tracts that dropped in or out of high impact
  cartographeR(dt=index_dt, map_varname = 'dropout', map_label = 'Dropout Direction',
               map_title = 'Tracts Newly Dropped/Added From High Impact',
               scale_type='drops')
  #show overall change in the index ranking
  cartographeR(dt=index_dt, map_varname = 'index_shift_capped', map_label = 'Index Change',
               map_title = 'Shift in Overall Index from V1 to V2',
               scale_type='div_man', scale_vals = div_colors)
  #show change in the in the theme specific ranking
  cartographeR(dt=theme_dt, map_varname = 'index_shift_capped', map_label = 'Index Change',
               map_title = 'Shift in Overall Index from V1 to V2',
               facet_var='theme', scale_type='div_man', scale_vals = div_colors)
  #show the current version of the overall index
  cartographeR(dt=index_dt, map_varname = 'index_new', map_label = 'Index',
               map_title = 'Overall Index (V2)', scale_type='div_man', scale_vals = cont_colors)
  #show the old version of the overall index
  cartographeR(dt=index_dt, map_varname = 'index_old', map_label = 'Index', 
               map_title = 'Overall Index (V1)', scale_type='cont')
  #show which tracts are highly impacted (current version)
  cartographeR(dt=index_dt, map_varname = 'impacted_new', map_label = 'Impacted', 
               map_title = 'Highly Impacted Tracts (V2)', scale_type='bin')
  #show which tracts are highly impacted (old version)
  cartographeR(dt=index_dt, map_varname = 'impacted_old', map_label = 'Impacted', 
               map_title = 'Highly Impacted Tracts (V1)', scale_type='bin')
  #see if the scaling method has any impact
  cartographeR(dt=index_dt, map_varname = 'scaling_effect', map_label = 'Effect of Scaling Method', 
               map_title = 'Scaling Method Impact on Final Index Ranking', scale_type='cont')
  
  #loop through changes in the measures binned by theme
  lapply(unique(measure_ranks$theme),
         cartographeR,
         shapefile=tract_sf,
         dt=measure_ranks[!(item%like%'Trans')], map_varname = 'measure_shift_capped', map_label= 'Index Change',
         map_title = 'Shift in Measure Ranking from V1 to V2',
         subset_var='theme',
         facet_var='item', scale_type='div_man', scale_vals = div_colors)
  
  #loop through changes in the measures
  lapply(unique(measure_raw$item),
         cartographeR,
         shapefile=tract_sf,
         dt=measure_raw, map_varname = 'measure_shift_raw', 
         subset_var='item', scale_type='cont')
  
  #loop through the measures themselves (both rank and raw)
  #rank
  lapply(unique(measure_raw$item),
         cartographeR,
         shapefile=tract_sf,
         dt=measure_ranks, map_varname = 'measure_new', map_label= 'Tract Ranking',
         map_title = 'Measure Ranking (V2)', 
         subset_var='item', scale_type='div_man', scale_vals = cont_colors)
  #raw
  lapply(unique(measure_raw$item),
         cartographeR,
         shapefile=tract_sf,
         dt=measure_raw, map_varname = 'measure_new_raw', map_label= 'Raw Measure',
         map_title = 'Measure Ranking (V2)', 
         subset_var='item', scale_type='cont')
  
  #loop through changes in the measure ranks
  lapply(unique(measure_ranks$item),
         cartographeR,
         shapefile=tract_sf,
         dt=measure_ranks, map_varname = 'measure_shift_capped', map_label= 'Index Change',
         map_title = 'Shift in Measure Ranking from V1 to V2',
         subset_var='item', scale_type='div_man', scale_vals = div_colors)

###explore changes
#reshape wide in order to run cor
corr_dt <- measure_ranks[, .(GEOID, item=paste0(theme, ': ', item), measure_new)] %>% 
  dcast(geocode~item, value.var=c('measure_new')) %>% 
  na.omit %>% 
  .[, -c('GEOID'), with=F] #drop the the geocode variable from this figure
  

#generate correlations
corr_dt <- corr_dt %>% cor

p_mat <- corr_dt %>% cor_pmat

ggcorrplot(corr_dt, type = "upper",
           colors = c("#6D9EC1", "white", "#E46726")) +
  theme(axis.text.x=element_text(angle=25,hjust=1))  

#save the plot
file.path(viz.dir, 'corrplot_v2.png') %>% ggsave(height=12, width=20)

#scaling effect basic scatterplot
ggplot(index_dt, aes(x=index_new, y=index_new_cal, color=scaling_effect) ) +
  geom_point(position='jitter') +
  #geom_hex(bins = 70) +
  scale_color_viridis(option='magma') +
  theme_bw() 

#life expectancy basic scatterplot
ggplot(index_dt, aes(x=index_new, y=le, color=impacted_new %>% as.factor) ) +
  geom_point(position='jitter') +
  #geom_hex(bins = 70) +
  scale_color_brewer('Dropouts', palette='Pastel1') +
  theme_bw() 

#create plot of life expectancy distributions
ggplot(index_dt, aes(index_new %>% as.factor, le)) + 
  geom_violin(aes(col = index_new %>% as.factor, fill = index_new %>% as.factor), alpha = 0.25)+
  labs(x = "Index", y = "Life Expectancy",
       title = "Relationship between EHD Index and Life Expectancy",
       subtitle = "Average Life Expectancy Between 2015-2019 (dotted line represents state average of 80.3 years)"
  ) +
  geom_vline(xintercept = 8.5, linetype='dashed', color='dark red') +
  geom_hline(yintercept = index_dt[1, le_state_average], linetype='dotted', color='dark blue') +
  scale_x_discrete() +
  scale_color_brewer('Index', palette='PuOr', direction=-1) +
  scale_fill_brewer('Index', palette='PuOr', direction=-1) +
  geom_boxplot(color = "gray20", width = 0.15, coef = 1.5) +
  annotate("text", x = 9.25, y = 99, label = "Highly Impacted Tracts", vjust = -0.5, color='dark red') +
  theme_minimal() +
  theme(legend.position = c(1, .99), legend.justification = c(1, 1),
        plot.margin = unit(c(0, 0, 0, 0), "in"))

#save the plot
file.path(viz.dir, 'le_violin_minimal.png') %>% ggsave(height=8, width=12)

#scatter the measure ranks to compare V1:v2
ggplot(measure_ranks, aes(x=measure_old, measure_new, color=dropout) ) +
  geom_point(position='jitter') +
  #geom_hex(bins = 70) +
  facet_wrap(~item) + 
  scale_x_continuous("Measure Ranking, V1") +
  scale_y_continuous("Measure Ranking, V2") +
  scale_color_brewer('Highly Impacted: \nDropouts', palette='Pastel1') +
  theme_minimal() 
#save the plot
file.path(viz.dir, 'measure_shift_scatters.png') %>% ggsave(height=8, width=12)

#***********************************************************************************************************************

# ---ANALYZE------------------------------------------------------------------------------------------------------------
####regression analysis####
#use logistic regression to find most influential vars
reg_dt <- index_dt[, .(GEOID, impacted_new, impacted_old, index_new, index_shift, dropout)] %>% 
  merge(measure_ranks[, .(GEOID, item, theme, measure_new, measure_old, measure_shift)], by='GEOID') %>% 
  .[, dropout_int := dropout==1] %>% 
  .[dropout_int %>% is.na, dropout_int:=0] #only regress on dropouts (not dropins)

mod1 <- glm(dropout_int ~ item:measure_shift, data=reg_dt, family='binomial')
stargazer(mod1, type='html') %>% 
  capture.output(file=file.path(out.dir, 'table_1.html'))

mod2 <- glm(index_shift ~ item:measure_shift, data=reg_dt, family='gaussian')
stargazer(mod1, mod2, 
          type='html',
          title="Regression Results", dep.var.labels=c("Dropouts","Index Movement"), 
          keep.stat="n", ci=TRUE, ci.level=0.90, single.row=TRUE) %>% 
  capture.output(file=file.path(out.dir, 'table_2.html'))

mod1 <- glm(dropout_int ~ item:measure_new, data=reg_dt, family='binomial')
mod2 <- glm(index_shift ~ item:measure_new, data=reg_dt, family='gaussian')
mod3 <- glm(index_new ~ item:measure_new, data=reg_dt, family='gaussian')
stargazer(mod1, mod2, mod3,
          type='html',
          title="Regression Results", dep.var.labels=c("Dropouts","Index Movement", 'New Index'), 
          keep.stat="n", ci=TRUE, ci.level=0.90, single.row=TRUE) %>% 
  capture.output(file=file.path(out.dir, 'table_3.html'))

mod1 <- glm(impacted_new ~ item:measure_shift, data=reg_dt, family='binomial')
mod2 <- glm(impacted_new ~ item:measure_new, data=reg_dt, family='binomial')
mod3 <- glm(impacted_old ~ item:measure_old, data=reg_dt, family='binomial')
stargazer(mod1, mod2, mod3,
          type='html',
          title="Regression Results", 
          keep.stat="n", ci=TRUE, ci.level=0.90, single.row=TRUE) %>% 
  capture.output(file=file.path(out.dir, 'table_4.html'))

mod1 <- glm(impacted_new ~ item:measure_shift, data=reg_dt, family='binomial')
stargazer(mod1, type='html') %>% 
  capture.output(file=file.path(out.dir, 'table_5.html'))


#save the plot
file.path(viz.dir, 'measure_scatters.png') %>% ggsave(height=8, width=12)

#
ggplot(reg_dt, aes(x=index_shift, measure_shift, color=dropout) ) +
  geom_point(position='jitter') +
  #geom_hex(bins = 70) +
  facet_wrap(~item) + 
  scale_color_brewer('Dropouts', palette='Pastel1') +
  theme_bw() 

#save the plot
file.path(viz.dir, 'measure_shift_scatters.png') %>% ggsave(height=8, width=12)
#***********************************************************************************************************************

# ---SCRAP -------------------------------------------------------------------------------------------------------------
##explore problems with wastewater data
cartographeR(dt=ranks_old$measure_raw, map_varname = 'rank_order', scale_type='cont',
             subset_var='item', subset_val='Wastewater Discharge')

cartographeR(dt=ranks_new$measure_raw, map_varname = 'rank_order', scale_type='cont',
             subset_var='item', subset_val='Wastewater Discharge')


#more complicated version of the violin plot (with labeled outliers)
ggplot(index_dt, aes(index_new %>% as.factor, le)) + 
  geom_violin(aes(col = index_new %>% as.factor, fill = index_new %>% as.factor), alpha = 0.25)+
  labs(x = "Index", y = "Life Expectancy (Years; 2015-2019)",
       # title = "Temperature rise in India",
       # subtitle = "Monthly temperature fluctuations between 1961-2019",
       # caption = "Data source: Kaggle.com"
  ) +
  geom_vline(xintercept = 2.5, linetype='dashed', color='dark red') +
  scale_x_discrete() +
  scale_color_brewer('Index', palette='PuOr', direction=-1) +
  scale_fill_brewer('Index', palette='PuOr', direction=-1) +
  geom_boxplot(color = "gray20", width = 0.15, coef = 1.5) +
  geom_text_repel(aes(label = outlier_lab), na.rm = TRUE, show.legend = F) +
  annotate("text", x = 2.1, y = 95, label = "Highly Impacted", vjust = -0.5, color='dark red') +
  theme_minimal() +
  theme(legend.position = c(1, 1), legend.justification = c(0, 0),
        #legend.text=element_text(size=10),
        plot.title = element_text(hjust=0.5), plot.margin=unit(c(0, 0, 0, 0), "in"))

#save the plot
file.path(viz.dir, 'le_violin.png') %>% ggsave(height=8, width=12)


ggplot(index_dt, aes(x = index_new %>% as.factor, y = le)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA) + 
  geom_boxplot(
    width = .25, 
    outlier.shape = NA
  ) +
  geom_point(
    size = 1.3,
    alpha = .3,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) + 
  coord_cartesian(clip = "off") +
  scale_x_discrete('Index') +
  scale_y_continuous('Life Expectancy (Years; 2015-2019)') +
  theme_minimal()

#save the plot
file.path(viz.dir, 'le_distribution.png') %>% ggsave(height=8, width=12)

# Bin size control + color palette
ggplot(measure_ranks, aes(x=item, fill=avg_shift) ) +
  geom_bar(aes(y=measure_shift), stat = "summary", fun = "mean") +
  #geom_hex(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  scale_y_continuous(trans = pseudolog10_trans) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=20,hjust=1))

# Bin size control + color palette
ggplot(index_dt, aes(x=GEOID, index_shift, color=index_shift) ) +
  geom_point(position='jitter') +
  #geom_hex(bins = 70) +
  scale_color_continuous(type = "viridis") +
  theme_bw() 

# Bin size control + color palette
ggplot(measure_raw[item=='Lead Risk From Housing'], aes(x=measure_old, measure_new, color=measure_shift) ) +
  geom_point(position='jitter') +
  #geom_hex(bins = 70) +
  scale_color_continuous(type = "viridis") +
  theme_bw() 

ggplot(measure_raw[item=='PM 2.5 -Diesel Emissions (Annual Tons/Km2)'], aes(x=measure_old, measure_new, color=measure_shift) ) +
  geom_point(position='jitter') +
  #geom_hex(bins = 70) +
  scale_color_continuous(type = "viridis") +
  theme_bw() 

# Clean Data Extract Data
# Apply Weights If Not = 1
# extract_measure_clean <- extract_measure %>%
#   mutate(IBLmeasure_rank = ifelse(Rank == "-1",NA,as.numeric(Rank)),
#          IBLmeasure_RankCalculatedValue = ifelse(RankCalculatedValue == "-1",NA,
#                                                  as.numeric(RankCalculatedValue)))%>%
#   select(GeoCode,ItemName,ThemeName,IBLmeasure_rank,IBLmeasure_RankCalculatedValue)
# 
# extract_theme_and_index_clean <- extract_theme %>% 
#   mutate(IBLtheme_rank = as.numeric(Rank),
#          IBLtheme_RankCalculatedValue = as.numeric(RankCalculatedValue),
#          IBLtheme_ThemeRank = as.numeric(ThemeRank),
#          IBLtheme_weights = ifelse(ThemeName == "Environmental Effects",.5,1)) %>%
#   select(GeoCode,ThemeName,IBLtheme_weights,IBLtheme_rank,IBLtheme_RankCalculatedValue,IBLtheme_ThemeRank) %>%
#   left_join((extract_index %>% 
#                mutate(IBLindex_rank = as.numeric(Rank),
#                       IBLindex_RankCalculatedValue = as.numeric(RankCalculatedValue)) %>%
#                select(GeoCode,IBLindex_rank,IBLindex_RankCalculatedValue,IndexName)),
#             by = "GeoCode")
# 
# # Calculate Ranks: Data Measures
# # Non-Missing Tract's WIth NA's Are Given A Measure Rank = -1
# ranks_measures <- extract_measure_clean %>%
#   mutate(ItemName = as.factor(ItemName)) %>%
#   group_by(ItemName)%>%
#   mutate(calc_measure_ordered_rank = rank(IBLmeasure_RankCalculatedValue,
#                                           na.last = "keep",
#                                           ties.method = "min"),
#          calc_measure_desired_bin_size = floor(sum(is.na(IBLmeasure_RankCalculatedValue)==F)/10),
#          calc_measure_rank_decimal = calc_measure_ordered_rank/calc_measure_desired_bin_size,
#          calc_measure_rank_integer = ifelse(calc_measure_rank_decimal > 10,10,
#                                             ceiling(calc_measure_rank_decimal)),
#          IBLmeasure_rank = ifelse(is.na(IBLmeasure_rank) == T, -1,IBLmeasure_rank),
#          IBLmeasure_RankCalculatedValue = ifelse(
#            is.na(IBLmeasure_RankCalculatedValue) == T, -1,
#            IBLmeasure_RankCalculatedValue),
#          calc_measure_rank_integer = ifelse(
#            is.na(calc_measure_rank_integer == T), -1, 
#            calc_measure_rank_integer),
#          calc_measure_rank_integer = ifelse(
#            is.na(calc_measure_rank_integer == T), -1, 
#            calc_measure_rank_integer),
#          QA_measure_ranks = ifelse(IBLmeasure_rank == calc_measure_rank_integer,
#                                    "pass","fail")) %>% ungroup()
# 
# # Calculate Ranks: Themes
# # Missing Tracts Are Added To Measures With A Measure Rank = 0
# ranks_theme <- (merge(
#   unique(extract_theme_and_index_clean$GeoCode),
#   unique(paste0(ranks_measures$ItemName,"@",ranks_measures$ThemeName))) %>%
#     rename(GeoCode = x, ItemName =y) %>%
#     mutate(ThemeName = as.factor(gsub(".*@","",ItemName)),
#            ItemName = gsub("@.*","",ItemName))) %>%
#   left_join(ranks_measures,
#             by = c("GeoCode","ItemName","ThemeName")) %>%
#   left_join(extract_theme_and_index_clean,
#             by = c("GeoCode","ThemeName")) %>%
#   mutate_all(~replace(., is.na(.), 0)) %>%
#   group_by(ThemeName) %>%
#   mutate(calc_theme_how_many_measures = sum(GeoCode == 53000000000)) %>%
#   ungroup() %>%
#   group_by(ThemeName,GeoCode) %>%
#   mutate(calc_theme_sum_measure_ranks = sum(calc_measure_rank_integer)) %>%
#   ungroup() %>%
#   distinct(ThemeName,GeoCode, .keep_all = T) %>%
#   group_by(ThemeName) %>%
#   mutate(calc_theme_avg_measure_ranks = round_half_up(calc_theme_sum_measure_ranks/calc_theme_how_many_measures,2),
#          calc_theme_ordered_rank = rank(calc_theme_avg_measure_ranks,
#                                         na.last = "keep",
#                                         ties.method = "min"),
#          calc_theme_desired_bin_size = floor(sum(is.na(calc_theme_avg_measure_ranks)==F)/10),
#          calc_theme_rank_decimal = calc_theme_ordered_rank/calc_theme_desired_bin_size,
#          calc_theme_rank_integer = ifelse(calc_theme_rank_decimal > 10,10,
#                                           ceiling(calc_theme_rank_decimal)),
#          QA_theme_ranks = ifelse(IBLtheme_rank == calc_theme_rank_integer,
#                                  "pass","fail")) %>% ungroup()
# 
# # Calculate Ranks: Index
# # File Details Above Indicates Theme "Environmental Effects" is Weighted At 0.5 And
# # The Rank Method = (("Environmental Effects"+"Environmental Exposures")/2)*
# #                     (("Socioeconomic Factors"+"Sensitive Populations")/2)
# # The First Sum In This Method = The Pollution Burden, The Second = Population Characteristics
# # The Theme Field We Start With = calc_theme_avg_measure_ranks
# ranks_index <- ranks_theme %>%
#   left_join((ranks_theme %>%
#                filter(ThemeName == "Environmental Effects" | ThemeName == "Environmental Exposures") %>%
#                mutate(calc_index_theme_avg_measure_ranks_weighted = calc_theme_avg_measure_ranks*
#                         IBLtheme_weights) %>%
#                group_by(GeoCode) %>%
#                summarize(calc_index_theme_sum_pollution_burden = sum(
#                  calc_index_theme_avg_measure_ranks_weighted))),
#             by="GeoCode") %>%
#   left_join((ranks_theme %>%
#                filter(ThemeName == "Sensitive Populations" | ThemeName == "Socioeconomic Factors") %>%
#                mutate(calc_index_theme_avg_measure_ranks_weighted = calc_theme_avg_measure_ranks*
#                         IBLtheme_weights) %>%
#                group_by(GeoCode) %>%
#                summarize(calc_index_theme_sum_population_character = sum(
#                  calc_index_theme_avg_measure_ranks_weighted))),by="GeoCode") %>%
#   distinct(GeoCode, .keep_all = T) %>%
#   mutate(calc_index_value_to_be_ranked = round_half_up((calc_index_theme_sum_pollution_burden/2) *
#                                                          (calc_index_theme_sum_population_character/2),2),
#          calc_index_ordered_rank = rank(calc_index_value_to_be_ranked,
#                                         na.last = "keep",
#                                         ties.method = "min"),
#          calc_index_desired_bin_size = floor(sum(is.na(calc_index_value_to_be_ranked)==F)/10),
#          calc_index_rank_decimal = calc_index_ordered_rank/calc_index_desired_bin_size,
#          calc_index_rank_integer = ifelse(calc_index_rank_decimal > 10,10,
#                                           ceiling(calc_index_rank_decimal)),
#          QA_index_ranks = ifelse(IBLindex_rank == calc_index_rank_integer,
#                                  "pass","fail")) %>%
#   relocate(IndexName, .after = last_col()) %>%
#   relocate(IBLindex_rank, .after = last_col()) %>%
#   relocate(IBLindex_RankCalculatedValue, .after = last_col())
# 
# # QA Messages
# QA_EHD <- paste("Measures - QA Results:\n # tracts with expected ranks = ",
#                 sum(ranks_measures$QA_measure_ranks == "pass"), " (",
#                 round(100*(sum(ranks_measures$QA_measure_ranks == "pass")/
#                              length(ranks_measures$QA_measure_ranks)),0),"%)",
#                 "\n # tracts with unexpected ranks = ",
#                 sum(ranks_measures$QA_measure_ranks == "fail"), " (",
#                 round(100*(sum(ranks_measures$QA_measure_ranks == "fail")/
#                              length(ranks_measures$QA_measure_ranks)),0),"%)",
#                 "\n\nThemes - QA Results:\n # tracts with expected ranks = ",
#                 sum(ranks_theme$QA_theme_ranks == "pass"), " (",
#                 round(100*(sum(ranks_theme$QA_theme_ranks == "pass")/
#                              length(ranks_theme$QA_theme_ranks)),0),"%)",
#                 "\n # tracts with unexpected ranks = ",
#                 sum(ranks_theme$QA_theme_ranks == "fail"), " (",
#                 round(100*(sum(ranks_theme$QA_theme_ranks == "fail")/
#                              length(ranks_theme$QA_theme_ranks)),0),"%)",
#                 "\n\nIndex - QA Results:\n # tracts with expected ranks = ",
#                 sum(ranks_index$QA_index_ranks == "pass"), " (",
#                 round(100*(sum(ranks_index$QA_index_ranks == "pass")/
#                              length(ranks_index$QA_index_ranks)),0),"%)",
#                 "\n # tracts with unexpected ranks = ",
#                 sum(ranks_index$QA_index_ranks == "fail"), " (",
#                 round(100*(sum(ranks_index$QA_index_ranks == "fail")/
#                              length(ranks_index$QA_index_ranks)),0),"%)")
# message(QA_EHD)

# For More QA Uncomment The Below And Run:
# table(ranks_measures$ItemName, ranks_measures$QA_measure_ranks)
# table(ranks_theme$ThemeName, ranks_theme$QA_theme_ranks)
# table(ranks_index$QA_index_ranks)

# For More QA Uncomment The Below And Run:
# Look At Each Data Measure Individually As It's Own Data Frame
# ranks_measures <- ranks_measures %>% split(f = ranks_measures$ItemName) %>% list2env(envir=.GlobalEnv)

#***********************************************************************************************************************