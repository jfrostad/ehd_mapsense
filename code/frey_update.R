# ----HEADER------------------------------------------------------------------------------------------------------------
# Author: JF
# Date: 04/08/2022
# Purpose: Prepare EHD rankings from raw data, explore and analyze sensitivity
# source("/homes/jfrostad/_code/ehd_mapsense/repo/code/frey.R", echo=T)
#***********************************************************************************************************************

# ----CONFIG------------------------------------------------------------------------------------------------------------
# clear memory
rm(list=ls())

#set opts
set.seed(98118)
options(scipen=999) #readability
#use cairo to render instead of quartz (quartz causes big slowdowns with geom_sf)
if(!identical(getOption("bitmapType"), "cairo") && isTRUE(capabilities()[["cairo"]])){
  options(bitmapType = "cairo")
}

#set control flow params
reload <- F #set true if you want to reprep all the data
remap <- F #set true if you want to redraw all the data maps
rerun_gsa <- F #set true if you are rerunning the global sensitivity analysis
gsa_cores <- 1

#version
gsa_version <- 1 #first draft with params=
gsa_version <- 2 #reduced noise
gsa_version <- 3 #switch to MARC instead of just avg. rerunning with fixed PCA sign

## Set core_repo location
user            <- Sys.info()['user']
main.dir         <- ifelse(Sys.info()["sysname"] == "Linux",
                           file.path('/homes', user, '_code/ehd_mapsense/'),
                           file.path('/Users', user, 'Documents/work/ehd_mapsense/'))
my_repo <- file.path(main.dir, 'repo')
setwd(my_repo)

#load packages
#TODO only relevant to running in linux on shared cluster
# package_lib    <- ifelse(.Platform$GUI=='RStudio',
#                          '/mnt/share/homes/jfrostad/_code/_lib/pkg_R_ehd',
#                          '/mnt/share/homes/jfrostad/_code/_lib/pkg_R_int')
# ## Load libraries and  MBG project functions.
# .libPaths(package_lib)

#TODO cleanup old packages
pacman::p_load(readxl, snakecase, janitor, data.table, naniar, stringr, magrittr, scales, Hmisc,
               ggtern, ggplot2, ggpubr, ggridges, ggrepel, ggdist, grid, gridExtra, RColorBrewer, #viz pkgs
               sf, viridis, farver, reldist, ggnewscale, ggallin, biscale, cowplot,
               tigris, tidycensus, ggcorrplot,
               broom.mixed, ggstance, jtools, factoextra, scam,
               COINr, randtoolbox, sensobol, #sens packages
               stargazer,
               parallel, pbmcapply,
               cluster, ggdendro, #HCA packages
               #caret, mlbench, randomForest, pls,
               zoo)

#***********************************************************************************************************************

# ----IN/OUT------------------------------------------------------------------------------------------------------------
###Input###
#raw data
code.dir <- file.path(my_repo, 'code')
data.dir <- file.path(my_repo, 'data')
data_extract_EHDv1 <- 'ehd_data_v1.xlsx' #TODO rename
data_extract_EHDv2 <- 'ehd_data_v3.xlsx' #TODO rename

###Output###
out.dir <- file.path(main.dir, 'output')
scratch.dir <- '/mnt/share/scratch/users/jfrostad/ehd_mapsense/output'
viz.dir  <- file.path(main.dir, 'viz')
vizdata.dir  <- file.path(my_repo, 'code/baldR/data')
#***********************************************************************************************************************

# ---FUNCTIONS----------------------------------------------------------------------------------------------------------
#source custom functions that are relevant to this module
file.path(code.dir, '_lib', 'mod_fx.R') %>% source
file.path(code.dir, '_lib', 'prep_fx.R') %>% source
file.path(code.dir, '_lib', 'sens_fx.R') %>% source
file.path(code.dir, '_lib', 'viz_fx.R') %>% source

##custom utilities##
#helper function to copy things out of R
writeExcel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
}

#label outliers statistically
isOutlier <- function(x) {
  quantile(x, 0.25, na.rm=T) - 1.5 * IQR(x, na.rm=T) | x > quantile(x, 0.75, na.rm=T) + 1.5 * IQR(x, na.rm=T)
}

#***********************************************************************************************************************

# ---PREP DATA----------------------------------------------------------------------------------------------------------
##read in and prep datasets for analysis##
if(reload) {
#names of themes have changed, make a map
item_map <- file.path(data.dir, 'ehd_map_theme_names.csv') %>% fread
gsa_sample_map <- file.path(data.dir, 'gsa_sample_map.xlsx') %>% 
  read_excel(sheet='Sheet1') %>% 
  as.data.table

#life expectancy data
le_dt <-  file.path(data.dir, 'le_at_birth_2015_2019.csv') %>% fread
setnames(le_dt, names(le_dt), 
         c('county', 'geocode', 'le', 'lower', 'upper'))
le_dt[, c('le', 'lower', 'upper') := lapply(.SD, as.numeric), .SDcols=c('le', 'lower', 'upper')]

#ihme life expectancy data
# regexp <- "[[:digit:]]+"
# king_le_dt <-  file.path(data.dir, 'ihme_kingco_le_1990_2014.csv') %>% 
#   fread %>% 
#   .[year_id==2010 & sex_id==3, .(location_name, val, lower, upper)] %>%  #TODO re-evaluate year
#   .[, tract_id := stringr::str_replace(location_name, 'King County Census Tract ', '')] %>% 
#   .[tract_id!='King County'] %>%  #drop the all county vals %>% 
#   .[, location_name := NULL] %>% 
#   setnames(., names(.), c('tract_le', 'tract_le_lower', 'tract_le_upper', 'tract_id'))

#bring in the CJEST data
cjest_dt <- file.path(data.dir, 'cjest_data_v1.csv') %>% 
  fread %>% 
  .[, c(1,2,3,16,20,21,23,108)] %>% 
  setnames(.,
           names(.),
           c('GEOID',
             'county_name',
             'state',
             'cjest_cat_exceeded',
             'cjest_impacted',
             'cjest_pct_imp',
             'cjest_pop',
             'cjest_le')) %>% 
  .[, GEOID := as.character(GEOID)]

eji_dt <- file.path(data.dir, 'eji_data.csv') %>% 
  fread %>% 
  .[, c(6,13,15,16)] %>% 
  setnames(.,
           names(.),
           c('GEOID',
             'eji_pop',
             'eji_spl',
             'eji_rpl'))%>% 
  .[, GEOID := as.character(GEOID)]

#bring in the urban rural classifications data
ruca_dt <- read_xlsx(file.path(data.dir, 'ruca_2010.xlsx'), skip=1, sheet = "Data") %>% 
  as.data.table %>% 
  .[, c(2,4:9)] %>% 
  setnames(.,
           names(.),
           c('state_abr',
             'GEOID',
             'ruca_code',
             'ruca_code_2',
             'ruca_pop',
             'ruca_sqm',
             'ruca_pop_dens')) %>%
  .[state_abr=='WA']

#use scheme 2 ruca classifications
ruca_dt[ruca_code==1, ruca_level := 1]
ruca_dt[ruca_code>1 & ruca_code<4, ruca_level := 2]
ruca_dt[ruca_code>=4 & ruca_code<7, ruca_level := 3]
ruca_dt[ruca_code>=7, ruca_level := 4]
ruca_dt[, ruca_level := factor(ruca_level,
                             labels=c('Urban',
                                      'Suburban',
                                      'Large Rural',
                                      'Small rural'))]

#also bring in the census tracts shapefile in order to do some cartography
#can be downloaded from the census website using tigris
tract_sf <- tracts('WA', year=2010, cb=T) %>% 
  st_transform(32148) %>% 
  erase_water(area_threshold = 0.9) %>% #intersect with water overlay and remove
  mutate('GEOID'=substring(GEO_ID, 10)) #remove the excess first 9 chr and rename GEOID

# #merge the GEOIDs onto the kingco le data
# king_le_dt <- tract_sf %>% 
#   as.data.table %>% 
#   .[COUNTY=='033', .(NAME, GEOID)] %>% 
#   merge(., king_le_dt,
#         by.x='NAME',
#         by.y='tract_id')

#use the water shapefile as an overlay
counties_list <- counties('WA', cb=T)
water_sf <- area_water('WA', counties_list$COUNTYFP %>% unique) %>% 
  st_simplify(preserveTopology = TRUE, dTolerance = 100)

#also overlay roads
road_sf <- primary_secondary_roads(state='WA') %>% 
  st_simplify(preserveTopology = TRUE, dTolerance = 100)

#also overlay places
places_sf <- places(state = 'WA', cb = T) 

#first read in and calculate all the ranks using custom function
ranks_old <- rankeR(dir=data.dir, path=data_extract_EHDv1, clean_names=item_map, debug=F)
ranks_new <- rankeR(dir=data.dir, path=data_extract_EHDv2, debug=F) 

##create comparisons##
#merge measures (old v. new) to compare
measure_ranks <- merge(ranks_old$measure[, .(GEOID, item, item_short, theme, level,
                                          rank_v1=measure_rank_integer)],
                    ranks_new$measure[, .(GEOID, item, theme, level,
                                          rank=measure_rank_integer),],
                    by=c('GEOID', 'item', 'theme', 'level'), all = T) %>% 
  .[, rank_shift := rank-rank_v1] %>% 
  .[, rank_shift_capped := rank_shift] %>% 
  .[rank_shift>=5, rank_shift_capped := 5] %>% #cap shift to max for plotting
  .[rank_shift<=-5, rank_shift_capped := -5]  %>% 
  .[, rank_shift := rank_shift_capped] %>% 
  .[, rank_shift_capped := NULL]
  

#merge raw measures (old v. new) to compare
measure_raw <- merge(ranks_old$measure_raw[, .(GEOID, item, item_short, theme, level,
                                             measure_v1=measure_rank_val)],
                       ranks_new$measure[, .(GEOID, item, theme, level,
                                             measure=measure_rank_val)],
                       by=c('GEOID', 'item', 'theme', 'level'), all=T) %>% 
  .[, measure_shift := measure-measure_v1] %>% 
  .[, measure_ratio := measure/measure_v1] %>% 
  .[measure_ratio %>% is.infinite, measure_ratio := NA] #0s create issues here
  
#merge both measures datasets
measure_dt <- merge(measure_ranks,
                    measure_raw[, .(GEOID, item, measure, measure_v1, measure_shift, measure_ratio)],
                    by=c('GEOID', 'item'))

#merge themes
theme_dt <- merge(ranks_old$theme[, .(GEOID, item, item_short, theme, level,
                                      rank_v1=theme_rank_integer)],
                  ranks_new$theme[, .(GEOID, item, theme, level,
                                      rank=theme_rank_integer)],
                  by=c('GEOID', 'item', 'theme', 'level')) %>% 
  .[, rank_shift := rank-rank_v1] %>% 
  .[, rank_shift_capped := rank_shift] %>% 
  .[rank_shift>=5, rank_shift_capped := 5] %>% #cap shift to max for plotting
  .[rank_shift<=-5, rank_shift_capped := -5]  %>% 
  .[, rank_shift := rank_shift_capped] %>% 
  .[, rank_shift_capped := NULL]

#merge indexes (old v. new) to compare
index_dt <- merge(ranks_old$index[, .(GEOID, item, item_short, theme, level,
                                      rank_v1=index_rank_integer)],
                    ranks_new$index[, .(GEOID, item, theme, level, 
                                        rank=index_rank_integer)],
                    by=c('GEOID', 'item', 'theme', 'level')) %>% 
  .[, rank_shift := rank-rank_v1] %>% 
  .[, rank_shift_capped := rank_shift] %>% 
  .[rank_shift>=5, rank_shift_capped := 5] %>% #cap shift to max for plotting
  .[rank_shift<=-5, rank_shift_capped := -5]  %>% 
  .[, rank_shift := rank_shift_capped] %>% 
  .[, rank_shift_capped := NULL]

##post estimations##
#identify dropout units based on impacted threshold of >8
threshold_val <- 9
index_dt[, impacted := 0]
index_dt[, impacted_v1 := 0]
index_dt[rank>=threshold_val, impacted := 1]
index_dt[rank_v1>=threshold_val, impacted_v1 := 1]
index_dt[, dropout := factor(impacted_v1 - impacted,
                                levels=c(-1, 0, 1),
                                labels=c('Addition',
                                         'Steady',
                                         'Dropout'))]
index_dt[, impacted_hierarchy := 0]
index_dt[impacted==1 & impacted==impacted_v1, impacted_hierarchy := 4]
index_dt[impacted==1 & impacted!=impacted_v1, impacted_hierarchy := 3]
index_dt[impacted==0 & impacted!=impacted_v1, impacted_hierarchy := 2]
index_dt[impacted==0 & impacted==impacted_v1, impacted_hierarchy := 1]
index_dt[, impacted_hierarchy := factor(impacted_hierarchy,
                                        levels=c(1, 2, 3, 4),
                                        labels=c('No v2 + No v1.1\n Certain Unimpacted',
                                                 'No v2 + Yes v1.1\n Uncertain Unimpacted',
                                                 'Yes v2 + No v1.1\n Uncertain Impacted',
                                                 'Yes v2 + Yes v1.1\n Certain Impacted'))]

#add life expectancy data
index_dt <- merge(index_dt, le_dt, by.x='GEOID', by.y='geocode')
index_dt[, le_state_average := mean(le, na.rm=T)]
index_dt[, county_le := paste0(le, ' years (', lower, '-', upper, ')')]
#index_dt <- merge(index_dt, king_le_dt, by='GEOID', all.x=T)

#add in the cjest and eji data
#TODO ID the cjest merge issue
index_dt <- merge(index_dt, cjest_dt, by='GEOID', all.x=T)
index_dt <- merge(index_dt, eji_dt, by='GEOID', all.x=T) %>% 
  unique(by='GEOID') #TODO for some reason duplicates GEOID==53057950100??

#add in the RUCA data
index_dt <- merge(index_dt, ruca_dt, by='GEOID', all.x=T)

#merge the dropouts back onto the other tables for graphing
measure_dt <- merge(measure_dt, index_dt[, .(GEOID, dropout, index=rank, impacted_hierarchy,
                                             impacted, impacted_v1, county_le, cjest_le,
                                             ruca_pop, ruca_sqm, ruca_pop_dens, ruca_level)], by='GEOID', all.x=T)
theme_dt <- merge(theme_dt, index_dt[, .(GEOID, dropout, index=rank, impacted_hierarchy,
                                         impacted, impacted_v1, county_le, cjest_le,
                                         ruca_pop, ruca_sqm, ruca_pop_dens, ruca_level)], by='GEOID', all.x=T)

#label the outliers
#index_dt[, outlier_lab := ifelse(isOutlier(le), name, NA_character_), by=index_new]

##output##
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
save(out, file=file.path(out.dir, 'all_data.RData'))

#reformat the data long to simplify for the mapping tool
#TODO - eventually i think it makes sense to adapt all future code to use this long version
dt <- list(index_dt[, -c('le', 'lower', 'upper', 
                         'tract_le_Lower', 'tract_le_upper',
                         'le_state_average', 'county'), with=F],
           theme_dt,
           measure_dt) %>% 
  rbindlist(use.names=T, fill=T)
dt[item=='Aggregated', item_short := 'Agg'] #TODO fix earlier
dt <- dt[!(is.na(item_short))] #rows that didn't have data for v1

#make it possible to use a logged scale on measure, which tends skewed
#dt[, measure_trans := log(measure + 0.001)]
#dt[item %like% '%' & measure>0, measure_trans := car::logit(measure+0.001)]
#dt[, measure_v1_trans := log(measure_v1 + 0.001)]
#dt[item %like% '%' & measure_v1>0, measure_v1_trans := car::logit(measure_v1+0.001)]

#save a  version of the data for the online mapping tool
out <- list(
  'dt'=dt,
  'tract_sf'=tract_sf,
  'water_sf'=water_sf,
  'road_sf'=road_sf,
  'places_sf'=places_sf
)

saveRDS(out, file=file.path(vizdata.dir, 'viz_data.RDS'))

#also save a csv of the lite data for edmund
out <- list(
  'dt'=dt,
  'tracts'=tract_sf
)
saveRDS(out, file=file.path(out.dir, 'lite_data.RDS'))
write.csv(dt, file=file.path(out.dir, 'lite_data.csv'))

} else file.path(vizdata.dir, 'viz_data.RDS') %>% readRDS %>% list2env(., globalenv())

#***********************************************************************************************************************

# ---MAP----------------------------------------------------------------------------------------------------------------
#create a manual diverging color scale to make sure that the index shifts are uniformly depicted
div_colors <- RColorBrewer::brewer.pal(11, 'RdBu') %>% rev
names(div_colors) <- dt[, unique(rank_shift) %>% sort]

#create a manual discrete color scale for the continuous index and make sure the color scale legend has all integers
cont_colors <- viridis::magma(n = 10)
names(cont_colors) <- 1:10

#update colors with DOH scale
doh_cont_colors <- c(
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
names(doh_cont_colors) <- 1:10

if (remap) {
#create a series of colored maps for the different vars in the dataset
  #denote tracts that dropped in or out of high impact
  cartographeR(dt=dt, map_varname = 'dropout', map_label = 'Dropout Direction',
               map_title = 'Tracts Newly Dropped/Added From High Impact (top 20% of ranks)',
               scale_type='drops')
  #show overall change in the index ranking
  cartographeR(dt=dt, map_varname = 'rank_shift', map_label = 'Rank Change',
               map_title = 'Shift in Overall Rank from V1.1 to V2.0',
               scale_type='div_man', scale_vals = div_colors)
  #show change in the in the theme specific ranking
  cartographeR(dt=dt, map_varname = 'rank_shift', map_label = 'Rank Change',
               lvl=2, #theme level
               map_title = 'Shift in Overall Rank from V1.1 to V2.0',
               facet_var='theme', scale_type='div_man', scale_vals = div_colors)
  #show the current version of the overall index
  cartographeR(dt=dt, map_varname = 'rank', map_label = 'Rank',
               map_title = 'Overall Rank (V2.0)', scale_type='cont_grad', scale_vals = cont_colors)
  #show the old version of the overall index
  cartographeR(dt=dt, map_varname = 'rank_v1', map_label = 'Rank',
               map_title = 'Overall Rank (V1.1)', scale_type='cont_man')
  #show which tracts are highly impacted (current version)
  cartographeR(dt=dt, map_varname = 'impacted', map_label = 'Impacted',
               map_title = 'Highly Impacted Tracts (V2.0)', scale_type='bin')
  #show which tracts are highly impacted (old version)
  cartographeR(dt=dt, map_varname = 'impacted_v1', map_label = 'Impacted',
               map_title = 'Highly Impacted Tracts (V1.1)', scale_type='bin')
  #see if the scaling method has any impact
  # cartographeR(dt=dt, map_varname = 'scaling_effect', map_label = 'Effect of Scaling Method', 
  #              map_title = 'Scaling Method Impact on Final Ranking', scale_type='cont')

  #loop through changes in the measures binned by theme
  lapply(unique(measure_ranks$theme),
         cartographeR,
         shapefile=tract_sf,
         dt=dt[!(item%like%'Trans')], map_varname = 'rank_shift', map_label= 'Rank Change',
         lvl=3, #most granular level
         map_title = 'Shift in Measure Ranking from V1.1 to V2.0',
         subset_var='theme',
         facet_var='item_short', scale_type='div_man', scale_vals = div_colors)

  #loop through the measures themselves (both rank and raw)
  #rank
  lapply(unique(dt$item),
         cartographeR,
         shapefile=tract_sf,
         dt=dt, map_varname = 'rank', map_label= 'Tract Ranking',
         lvl=3, #most granular level
         map_title = 'Measure Ranking (V2.0)', 
         subset_var='item', scale_type='cont_man', scale_vals = cont_colors)
    #raw
  lapply(unique(measure_raw$item),
         cartographeR,
         shapefile=tract_sf,
         dt=measure_raw, map_varname = 'measure_new_raw', map_label= 'Raw Measure',
         map_title = 'Measure Ranking (V2.0)', 
         subset_var='item', scale_type='cont_man')
  
  #loop through changes in the measure ranks
  lapply(unique(measure_ranks$item),
         cartographeR,
         shapefile=tract_sf,
         dt=measure_ranks, map_varname = 'measure_shift_capped', map_label= 'Index Change',
         map_title = 'Shift in Measure Ranking from V1 to V2',
         subset_var='item', scale_type='div_man', scale_vals = div_colors)
  
  
  #some plots requested by esther for some community meetings
  #show the current version of the overall index
  impacted_colors <- viridis::plasma(n = 10)[9:10]
  names(impacted_colors) <- 9:10
  
  cartographeR(dt=dt[rank>8], map_varname = 'rank', map_label = 'Rank',
               map_title = 'Overall Rank (V2.0)', scale_type='cont_man', scale_vals = impacted_colors,
               tag='community_meeting_910_v2')
  #show the old version of the overall index
  cartographeR(dt=dt[rank_v1>8], map_varname = 'rank_v1', map_label = 'Rank',
               map_title = 'Overall Rank (V1.1)', scale_type='cont_man', scale_vals = impacted_colors,
               tag='community_meeting_910_v1')
  #same with 8s included
  impacted_colors <- viridis::plasma(n = 10)[8:10]
  names(impacted_colors) <- 8:10
  cartographeR(dt=dt[rank>7], map_varname = 'rank', map_label = 'Rank',
               map_title = 'Overall Rank (V2.0)', scale_type='cont_man', scale_vals = impacted_colors,
               tag='community_meeting_8910_v2')
  cartographeR(dt=dt[rank_v1>7], map_varname = 'rank_v1', map_label = 'Rank',
               map_title = 'Overall Rank (V1.1)', scale_type='cont_man', scale_vals = impacted_colors,
               tag='community_meeting_8910_v1')
  
}
#***********************************************************************************************************************

# ---COINr---------------------------------------------------------------------------------------------------------------
#setup COINr to run a sensitivity analysis across your param samples
#prepare the data and metadata to coin standards
#data
ehd_coin_dt <- dt[level==3, .(uCode=GEOID, item_short, measure)] %>% 
  .[item_short %like% '%', measure := measure /100] %>% #convert to pct
  .[, item_short := str_to_lower(item_short) %>% #cleanup so we can reshape wide
      str_replace_all(., ' \\(%\\)', '') %>% 
      str_replace_all(., ' ', '_') %>% 
      str_replace_all(., '\\(rsei\\)', 'rsei') %>%
      str_replace_all(., '2.5', '25')] %>% 
  .[measure==-1, measure := NA] %>%  #-1s are code for missing val
  dcast(uCode~..., value.var='measure') #shape wide

#add in the new indicators for sens
#pesticides
pest_dt<- file.path(data.dir, 'indicators', 'pesticides.csv') %>% 
  fread %>% 
  .[, .(uCode=CT_GEOID %>% as.character, pesticides=pesticide_lbs_mile2)]

#asthma
asthma_dt <- file.path(data.dir, 'indicators', 'asthma.rds') %>% 
  readRDS() %>% 
  as.data.table %>% 
  .[, .(uCode=TRACT %>% as.character, asthma=Modeled_Age_Adj_Rate)]

#wildfire smoke
smoke_dt<- file.path(data.dir, 'indicators', 'wildfire_smoke.csv') %>% 
  fread %>% 
  .[, .(uCode=County %>% as.character, wildfire_smoke=`Cumulative Smoke Score (2016-2022)`)]

#merge
ehd_coin_dt <- ehd_coin_dt %>% 
  merge(pest_dt, by='uCode') %>% 
  merge(asthma_dt, by='uCode') %>% 
  merge(smoke_dt, by='uCode')

#new metadata
ehd_coin_meta_new <- data.table(
  iCode=c('pesticides', 'asthma', 'wildfire_smoke'),
  iName=c('Pesticide Use', 'Asthma Rates', 'Wildfire Smoke Score'),
  Direction=c(1,1,1),
  Weight=c(1,1,1),
  Parent=c('environmental_exposures', 'sensitive_populations', 'environmental_exposures'),
  Level=c(1,1,1),
  Type='Indicator'
)

#meta
ehd_coin_meta <- dt[, .(level, item_short, item, theme, 
                        Direction=1,
                        Weight=1,
                        Parent=theme)] %>% 
  unique %>% 
  .[level==2, item_short := theme] %>% 
  .[level==2, item := theme] %>% 
  .[level==1, item_short := 'ehd_rank'] %>%
  .[level==1, item := 'EHD Ranking'] %>%
  #COINr has the levels structure inverted
  .[level==1, Level := 4] %>% 
  .[level==2, Level := 2] %>% 
  .[level==3, Level := 1] %>% 
  #cleanup names to match the wide dt
  .[, item_short := str_to_lower(item_short) %>% 
      str_replace_all(., ' \\(%\\)', '') %>% 
      str_replace_all(., ' ', '_') %>% 
      str_replace_all(., '\\(rsei\\)', 'rsei') %>%
      str_replace_all(., '2.5', '25')] %>% 
  #cleanup names to match the wide dt
  .[, Parent := str_to_lower(Parent) %>% 
      str_replace_all(., ' ', '_')] %>% 
  #need to set up some of the hierarchy manually
  .[level==1, Parent := NA] %>% 
  .[level==2, Parent := 'population_chars'] %>% 
  .[level==2 & item_short %like% 'environmental', Parent := 'pollution_burden'] %>% 
  .[level==2 & item_short %like%'effects', Weight := .5] %>% 
  list(.,
       data.table(Level=c(3,3),
                  item_short=c('population_chars', 'pollution_burden'),
                  item=c('Population Characteristics', 'Pollution Burden'),
                  Parent='ehd_rank',
                  Direction=1,
                  Weight=1)) %>% 
  rbindlist(fill=T) %>% 
  #COINr typing (base level only is indicator)
  .[, Type := 'Aggregate'] %>% 
  .[level==3, Type := 'Indicator'] %>% 
  #cleanup
  setnames(., 
           c('item_short', 'item'),
           c('iCode', 'iName')) %>% 
  .[, iName := as.character(iName)] %>%  #required for meta type
  .[, `:=` (theme=NULL, level=NULL)]

#add the new inds
ehd_coin_meta <- ehd_coin_meta %>% 
  list(., ehd_coin_meta_new) %>% 
  rbindlist

#build coin
coin <- new_coin(iData = ehd_coin_dt,
                 iMeta = ehd_coin_meta,
                 level_names = c("Indicator", "Pillar", "Sub-index", "Index"))

#treat data 
#TODO should build a version that doesn't treat at all to keep default
#needed to use this in order to analyze sens
coin <- Treat(coin, dset = "Raw", global_specs = list(f1 = 'winsorise'))

#normalize the data distribution by percentile ranking
#TODO note that this method is by percentiles instead of deciles
#could round in a custom function
coin <- Normalise(coin, dset = "Treated",
                  global_specs=list(f_n='n_brank'))

#aggregate using base function to compare to our results
coin <- Aggregate(coin, dset = "Normalised", 
                  w='Original',
                  flatten_hierarchy = F,
                  #note that we take the arith mean until the last aggregation, where we multiply for risk score
                  f_ag = c("a_amean", "a_amean", "prod")) 

#compare our results
coin_results_dt <- coin$Data$Aggregated %>% as.data.table

#plot scatter
plot_dt <- coin_results_dt[, .(uCode, ehd_rank)] %>% 
  merge(., dt[level==1, .(GEOID, rank)], by.x='uCode', by.y='GEOID') %>% 
  .[, index_rank_order := frank(ehd_rank, 
                                na.last='keep', 
                                ties.method = 'min')] %>% 
  .[, index_bin_size := floor(sum(!is.na(ehd_rank))/10)] %>% 
  #bin the ranks and number them by rounding up
  .[, coin_brank := (index_rank_order / index_bin_size) %>% ceiling] %>% 
  .[coin_brank>10, coin_brank := 10] %>% 
  #also calculate equal intervals
  .[, int_brank := cut(ehd_rank, breaks = 10, labels = 1:10)] %>% 
  .[, coin_int_rank := int_brank %>% as.integer] %>% 
  #also calculate z scores
  .[, coin_zrank := (rank - mean(rank)) / sd(rank)] %>% 
  .[, impacted := 0] %>% 
  .[, zimpacted := 0] %>%
  .[coin_brank>8, impacted := 1] %>% 
  .[coin_zrank>quantile(coin_zrank, p=zscore_p80), zimpacted:=1] %>% 
  .[, diff := coin_brank - rank] %>% 
  .[, class_diff := coin_brank - coin_int_rank] %>% 
  .[, impacted_diff := impacted==zimpacted] %>% 
  .[, GEOID := uCode] %>% 
  .[, level := 1] #for mapping

#investigate where the points are off, generally agreement seems good
ggplot(plot_dt, aes(x=coin_brank, y=rank, color=diff)) +
  geom_point(position='jitter') +
  scale_color_viridis() +
  scale_x_continuous('EHD Rank v3') +
  scale_y_continuous('EHD Rank v2') +
  theme_minimal()

#investigate where the points are off, generally agreement seems good
ggplot(plot_dt, aes(x=forcats::fct_reorder(coin_brank %>% as.factor, index_rank_order), 
                    y=coin_zrank, color=impacted_diff)) +
  geom_point(position='jitter') +
  #scale_color_viridis() +
  scale_x_discrete('Binned Ranking') +
  scale_y_continuous('Equal Interval Ranking') +
  theme_minimal()

#coin brank vs coin intervals mapping
cartographeR(dt=plot_dt, map_varname = 'coin_brank', map_label = 'EHD Index \n(Deciles)',
             map_title = '',
             tag = 'coin_results', scale_type='cont_man', scale_vals = cont_colors,
             get_plot=T)

cartographeR(dt=plot_dt, map_varname = 'diff', map_label = 'EHD Index \n(Deciles)',
             map_title = '',
             tag = 'coin_results', scale_type='div_man', scale_vals = div_colors,
             get_plot=T)

cartographeR(dt=plot_dt, map_varname = 'coin_int_rank', map_label = 'EHD Index \n(Equal Intervals)',
             map_title = '',
             tag = 'coin_results', scale_type='cont_man', scale_vals = cont_colors)

cartographeR(dt=plot_dt, map_varname = 'coin_zrank', map_label = 'EHD Index \n(Equal Intervals)',
             map_title = '',
             tag = 'coin_results', scale_type='cont')

#***********************************************************************************************************************

# ---ALT COINs----------------------------------------------------------------------------------------------------------
if(FALSE) {

#adjustments/tests
#remove the new indicators
# remove new indicators and regenerate the coin
coin_old <- change_ind(coin, drop = c("asthma", 'wildfire_smoke', 'pesticides'), regen = TRUE)
#add/drop inds to test
coin_asthma <- change_ind(coin_old, add = 'asthma', regen=T)
coin_smoke <- change_ind(coin_old, add = 'wildfire_smoke', regen=T)
coin_pest <- change_ind(coin_old, add = 'pesticides', regen=T)
coin_poc <- change_ind(coin_old, drop = c("people_of_color"), regen = T) 

# remove two indicators and regenerate the coin
coin_list <- list(
  'coins'=list(coin_old, coin_asthma, coin_smoke, coin_pest, coin_poc),
  'names'=list('ehd_rank', 'ehd_asthma', 'ehd_smoke', 'ehd_pest', 'ehd_white')
)


extractCoin <- function(i, coin_list) {
  
  this_coin <- coin_list$coins[[i]]
  this_name <- coin_list$names[[i]]
  
  dt <- this_coin$Data$Aggregated %>% 
    as.data.table %>% 
    .[, .(uCode, ehd_rank, var=this_name)] %>% 
    .[, ehd_rank := n_brank(ehd_rank)] %>% 
    return
  
}

sens_dt <- lapply(1:length(coin_list$coins), extractCoin, coin_list=coin_list) %>% 
  rbindlist %>% 
  dcast(uCode~var, value.var = 'ehd_rank') %>% 
  .[, asthma_diff := ehd_asthma - ehd_rank] %>%
  .[, smoke_diff := ehd_smoke - ehd_rank] %>%
  .[, pest_diff := ehd_pest - ehd_rank] %>%
  .[, white_diff := ehd_white - ehd_rank] %>%
  .[, level := 1] %>% 
  setnames('uCode', 'GEOID')
  
cartographeR(dt=sens_dt, map_varname = 'asthma_diff', map_label = 'Difference w/ \nAsthma',
             map_title = '',
             tag = 'coin_results', scale_type='div_man', scale_vals = div_colors,
             get_plot=T)

cartographeR(dt=sens_dt, map_varname = 'smoke_diff', map_label = 'Difference w/ \nWildfire Smoke',
             map_title = '',
             tag = 'coin_results', scale_type='div_man', scale_vals = div_colors,
             get_plot=T)

cartographeR(dt=sens_dt, map_varname = 'pest_diff', map_label = 'Difference w/ \nPesticide Use',
             map_title = '',
             tag = 'coin_results', scale_type='div_man', scale_vals = div_colors,
             get_plot=T)

cartographeR(dt=sens_dt, map_varname = 'white_diff', map_label = 'Difference w/o \n% PoC',
             map_title = '',
             tag = 'coin_results', scale_type='div_man', scale_vals = div_colors,
             get_plot=T)

#get some stats!!
testChange <- function(dt, this_col, thresh=8) {
  dt[, test_col := get(this_col)]
  message('modifying ', this_col, ' changes the ranking for:')
  dt[, sum(ehd_rank!=test_col)]  %>% print
  dt[, sum(ehd_rank!=test_col)/.N]  %>% print
  message('modifying ', this_col, ' changes the impact for:')
  dt[, sum((ehd_rank>=thresh)!=(test_col>=thresh))] %>% print
  dt[, sum((ehd_rank>=thresh)!=(test_col>=thresh))/.N]  %>% print
  dt[, test_col := NULL]
}

testChange(dt=sens_dt, this_col='ehd_white')
testChange(dt=sens_dt, this_col='ehd_smoke')
testChange(dt=sens_dt, this_col='ehd_pest')
testChange(dt=sens_dt, this_col='ehd_asthma')

#output the PoC index so it can be used by the agency
file.path(out.dir, 'ehd_ranks_without_poc.csv') %>% 
  write.csv(sens_dt[, .(GEOID, ehd_rank_wo_poc=ehd_white)], file=.)

#TODO compare it to LBW to show how some indicators are much more sensitive
  
#TODO
#impute to fix issues in the wastewater data

#treat data 
#TODO should build a version that doesn't treat at all to keep default
#needed to use this in order to analyze sens
#coin_alt <- Treat(coin, dset = "Raw", global_specs = list(f1 = 'winsorise'))

#normalize to test the zscore method
norm_alts <- list(
  minmax=list(f_n = "n_minmax", f_n_para = list(c(0.001,1))),
  z=list(f_n = "n_zscore_log", f_n_para = list(c(0,.5))),
  pct=list(f_n = 'n_prank_log'),
  dec=list(f_n = 'n_brank')
)

coin_alt <- Normalise(coin, dset = "Raw",
                  global_specs=norm_alts[['n_brank']])

#see how different norm methods change the data distributions
#note how ranking introduces nonlinearity
plot_dist(coin, dset='Normalised', iCodes = unique(coin$Data$Meta$iCode), type='Dot')

#compare inds
# compare one of the indicators
COINr::plot_scatter(coin_alt, dset = c("Raw", "Normalised"),
                    iCodes = "ozone_concentration")

#aggregate using base function to compare to our results
coin_alt <- Aggregate(coin_alt, dset = "Normalised", 
                  w='Original',
                  #note that we take the arith mean until the last aggregation, where we multiply for risk score
                  #f_ag = c("a_amean", "a_amean", "prod"),
                  f_ag='pca') 

#compare our results
#coin_results_dt <- get_results(coin, dset = "Aggregated", tab_type = "Aggs") %>% as.data.table
alt_coin_results_dt <- coin_alt$Data$Aggregated %>% 
  as.data.table %>% 
  .[, level:=1] %>% 
  .[, GEOID := uCode]

cartographeR(dt=alt_coin_results_dt, map_varname = 'ehd_rank', map_label = 'EHD Index \n(Equal Intervals)',
             map_title = '',
             tag = 'coin_results_pca', scale_type='cont')

#plot scatter
plot_dt <- alt_coin_results_dt[, .(uCode, ehd_rank)] %>% 
  merge(., dt[level==1, .(GEOID, rank)], by.x='uCode', by.y='GEOID') %>% 
  .[, coin_rank := n_brank(ehd_rank)] %>% 
  .[, diff := coin_rank - rank]

#investigate where the points are off, generally agreement seems good
ggplot(plot_dt, aes(x=coin_rank, y=rank, color=diff)) +
  geom_jitter(height=.5, width=.5) +
  scale_color_gradientn(colors=RColorBrewer::brewer.pal(10, 'PuOr')) +
  theme_minimal()

#investigative plots
#understand the agreement within groups
getCronbach(ASEM, dset = "Aggregated", icodes = "Conn", aglev = 2)


#remove elements
l_res <- remove_elements(coin, Level = 1, dset = "Aggregated", iCode = "ehd_rank")

# #consider another difference plot with more granularity
# ggplot(SA_res$diffs_dt[param=='Norm'], aes(x=sample, y=average_diff)) +
#   geom_bar(stat='identity')

}
#***********************************************************************************************************************

# ---GSA SPECS----------------------------------------------------------------------------------------------------------
#now use our COIN to run a global sensitivity analysis across our parameters
#first we build out our parameters for each step
# # TREATMENT OPTIONS
# l_winmax <- list(Address = "$Log$Treat$global_specs$f1_para$winmax",
#                  Distribution = 1:5,
#                  Type = "discrete")
# 
# treat_alts <- list(
#   #list(f_t = "log_CT", f_t_para = list(na.rm = TRUE)),
#   list(f_t = "winsorise", f_t_para = list(na.rm = TRUE,
#                                           winmax = 5,
#                                           skew_thresh = 2,
#                                           kurt_thresh = 3.5,
#                                           force_win = FALSE)),
#   list(f_t = "winsorise", f_t_para = list(na.rm = TRUE,
#                                           winmax = 3,
#                                           skew_thresh = 2,
#                                           kurt_thresh = 3.5,
#                                           force_win = FALSE)),
#   list(f_t = "winsorise", f_t_para = list(na.rm = TRUE,
#                                           winmax = 1,
#                                           skew_thresh = 2,
#                                           kurt_thresh = 3.5,
#                                           force_win = FALSE))
# )
# 
# l_treatfx <- list(Address = "$Log$Treat$global_specs",
#                   Distribution = treat_alts,
#                   Type = "discrete")

#Normalization OPTIONS:
#default we define the two alternatives: minmax or zscore (along with respective parameters)
norm_alts <- list(
  list(f_n = "n_minmax", f_n_para = list(c(0.001,1))),
  list(f_n = "n_zscore_log", f_n_para = list(c(0,.5))),
  list(f_n = 'n_prank_log'),
  list(f_n = 'n_brank')
)

l_norm <- list(Address = "$Log$Normalise$global_specs",
               Distribution = norm_alts,
               Type = "discrete")

# get nominal weights
w_nom <- coin$Meta$Weights$Original

#build an equally weighted version
# copy original weights
w_equal <- copy(w_nom)
w_equal$Weight <- 1
coin$Meta$Weights$EqualWeights <- w_equal

#build a version weighted by the first comp of PCA
pca_results <- get_PCA(coin, dset = "Raw", Level = 1, out2 = "list", by_groups=F, weights_to='PCAComp1', impute=T)
coin$Meta$Weights$PComp1 <- pca_results$Weights
coin$Meta$Weights$PComp1$Weight <- abs(pca_results$Weights$Weight) #TODO what to do about negative loadings? arbitrary

#build a version that uses the inverted correlations of the ranked data in order to downweight highly correlated inds
#see spiegelhalter 2012
corr_results <- get_corr(coin, dset='Normalised', Levels = 1)  %>% 
  as.data.table() %>% 
  .[Var1!=Var2] %>% #drop the equivalent ones
  .[, Weight := 1/abs(sum(Correlation, na.rm = T)), by=Var1] %>% #TODO how to account for negative correlations?
  unique(by='Var1') %>% 
  .[, .(iCode=Var1, Level=1, Weight)]

coin$Meta$Weights$InvCorr <- list(
  corr_results,
  w_nom[Level>1] #use the same weights for higher levels
) %>% rbindlist

#drop indicators to proxy indicator selection
# get 100 replications
drop_cols <- names(coin$Data$Raw)[-1]
l_ind<- list(Address = "$Log$new_coin$exclude",
                Distribution = drop_cols,
                Type = "discrete")

#simulate random noise into the raw var for measurement error
# get 100 replications
noisy_dts <- get_noisy_dt(dt = coin$Data$Raw, noise_factor = .25, Nrep = 100)
l_noise <- list(Address = "$Log$new_coin$iData",
                Distribution = noisy_dts,
                Type = "discrete")

# component of SA_specs for weights
l_weights <- list(Address = "$Log$Aggregate$w",
                  Distribution = c('Original', 'EqualWeights', 'PComp1', 'InvCorr'),
                  Type = "discrete")

## AGG OPTIONS
#different types of aggregation formula functions
agg_alts <- list(
  f_ag = c('pca'),
  f_ag = c("a_amean", "a_amean", "prod"),
  f_ag = c("a_amean", "a_amean", "a_amean"),
  f_ag = c("a_gmean", "a_gmean", "prod"),
  f_ag = c("a_gmean", "a_gmean", "a_gmean")
)
l_agg <- list(Address = "$Log$Aggregate$f_ag",
              Distribution = agg_alts,
              Type = "discrete")

#flatten the hierarchy or keep it hierarchical
l_hier <- list(Address = '$Log$Aggregate$flatten_hierarchy',
               Distribution = c(T,F),
               Type= 'discrete')


#test different classification methods
l_class <- list(Address = '$classification',
               Distribution = c('equal_int', 'deciles', 'zscore'),
               Type= 'discrete')

#build sampling framework
SA_specs <- list(
  Selection = l_ind,
  Measurement_Error = l_noise,
  Norm = l_norm,
  Weighting_Scheme = l_weights,
  Aggregation_Formula = l_agg,
  Classification = l_class,
  Hierarchy = l_hier
)
  
##Run Sens/Unc Analysis
if(rerun_gsa) {
 
  SA_res <- get_sensitivity(coin, SA_specs = SA_specs, N = 5000, SA_type = "SA", use_branks=F,
                            dset = "Aggregated", iCode = "ehd_rank", Nboot = 100, ncores = gsa_cores, sock=T)
  saveRDS(SA_res, file=file.path(scratch.dir, paste0('gsa_output_v', gsa_version, '.RDS')))
  
  browser()
  
} else SA_res <- file.path(scratch.dir, paste0('gsa_output_v', gsa_version, '.RDS')) %>% readRDS
#***********************************************************************************************************************

# ---GSA RESULTS--------------------------------------------------------------------------------------------------------
#make some plot and stat DTs from the GSA outputs
#for ranks
rank_stats_dt <- SA_res$RankStats %>% 
  as.data.table %>%
  setnames('uCode', 'GEOID') %>% 
  merge(., dt[level==1, .(GEOID, impacted_hierarchy, county_name, ruca_level, ruca_pop, ehd_rank=rank)], by='GEOID') %>% 
  .[, range := Q95-Q5] %>% 
  .[, deviation := Median-Nominal_rank] %>% 
  .[, level := 1] 

#for binranks
brank_stats_dt <- SA_res$BRankStats %>% 
  as.data.table %>% 
  .[, range := Q95-Q5] %>% 
  .[, deviation := Median-Nominal_rank] %>% 
  .[, level := 1] %>% 
  setnames('uCode', 'GEOID') %>% 
  merge(., dt[level==1, .(GEOID, impacted_hierarchy, county_name, ruca_level, ruca_pop)], by='GEOID')

#generate summary stats
#do some birdwatchin
rank_stats_dt[, weighted.mean(deviation %>% abs, ruca_pop)] #avg deviation in ranking
rank_stats_dt[, median(range), by=ehd_rank][order(ehd_rank)] #avg deviation across bin groups from base EHD

rank_stats_dt[, weighted.mean(deviation %>% abs, ruca_pop), by=ehd_rank][order(ehd_rank)] %>% 
  ggplot(., aes(ehd_rank, V1, fill=V1)) + 
  geom_bar(stat='identity') + 
  scale_y_continuous('Average Rank Deviation') + 
  scale_x_continuous('Baseline EHD Ranking', breaks=pretty_breaks()) + 
  scale_fill_viridis() +
  theme_minimal()

file.path(viz.dir, 'mean_rank_deviation_binned.png') %>% ggsave(height=8, width=12)

#***********************************************************************************************************************

# ---GSA FIGURES--------------------------------------------------------------------------------------------------------
#make plots for the manuscript
# Figure 2a ----------------------------------------------------------
#make some combined plots for figure 2
#first make a general layout for this plot type
lay <- rbind(c(2,2,1,1,1,1,1,1,1),
             c(2,2,1,1,1,1,1,1,1),
             c(2,2,1,1,1,1,1,1,1),
             c(2,2,1,1,1,1,1,1,1),
             c(2,2,1,1,1,1,1,1,1),
             c(2,2,1,1,1,1,1,1,1))

#figure 2a
#first make a histogram of county level deviations
dev_hist <-
  rank_stats_dt[, .(V1=weighted.mean(deviation, ruca_pop),
                    V2=weighted.mean(deviation, ruca_pop) %>% abs), by=county_name][order(county_name)] %>% 
  .[, county_short := str_remove(county_name, ' County')] %>% 
  ggplot(., aes(forcats::fct_reorder(county_short %>% as.factor, V1), V2, fill=V1)) + 
  geom_bar(stat='identity') + 
  scale_x_discrete('') +
  scale_y_continuous('Average Rank Change') + 
  scale_fill_gradient2(guide='none') +
  coord_flip() +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
file.path(viz.dir, 'mean_rank_deviation_county.png') %>% ggsave(height=8, width=12)

#then map the tract level deviations
dev_map <- 
  cartographeR(dt=rank_stats_dt, map_varname = 'deviation', map_label = 'Median Change',
               map_title = '',
               tag = 'gsa_results',
               scale_type='cont_grad',
               get_plot=T)

#now combine into a single map
all_grobs <- list(dev_map+
                    theme(plot.margin = unit(c(0,0,0,-1), units = "cm"),
                          legend.position = "none"), 
                  dev_hist)
plot <- arrangeGrob(grobs=all_grobs, layout_matrix=lay, 
                    top=textGrob("Baseline EHD Rank - Median Simulated Rank", 
                                 gp = gpar(fontsize=17))
) %>% 
  grid.arrange

ggsave(plot=plot, filename=file.path(viz.dir, 'fig_2a.png'),
       width=12, height=8, units='in', dpi=900)

# Figure 2b ----------------------------------------------------------
#first make a histogram of county level uncertainty
range_hist <-
  rank_stats_dt[, weighted.mean(range, ruca_pop), by=county_name][order(county_name)] %>% 
  .[, county_short := str_remove(county_name, ' County')] %>% 
  ggplot(., aes(forcats::fct_reorder(county_short %>% as.factor, V1), V1, fill=V1)) + 
  geom_bar(stat='identity') + 
  scale_x_discrete('') +
  scale_y_continuous('Average Uncertainty Interval (95%)') + 
  scale_fill_viridis(guide='none', option='magma') +
  coord_flip() +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 35, vjust = 0.5, hjust=1))

file.path(viz.dir, 'mean_uncertainty_county.png') %>% ggsave(height=8, width=12)

#then map the tract level deviations
range_map <- 
  cartographeR(dt=rank_stats_dt, map_varname = 'range', map_label = 'Uncertainty \n(95% Interval)',
               map_title = '',
               tag = 'gsa_results',
               scale_type='cont', 
               get_plot = T)

all_grobs <- list(range_map+
                    theme(plot.margin = unit(c(0,0,0,-1), units = "cm"),
                          legend.position = "none"), 
                  range_hist)


plot <- arrangeGrob(grobs=all_grobs, layout_matrix=lay, 
                    top=textGrob("EHD Index Uncertainty Interval", 
                                    gp = gpar(fontsize=17))
) %>% 
  grid.arrange

ggsave(plot=plot, filename=file.path(viz.dir, 'fig_2b.png'),
       width=12, height=8, units='in', dpi=900)

# Figure 2c ----------------------------------------------------------
#make a biplot showing the uncertainty
rank_stats_dt$ehd_bin <- cut(rank_stats_dt$ehd_rank, breaks = c(0,2,8,10), include.lowest = TRUE)
bi_data <- bi_class(rank_stats_dt, x = ehd_bin, y = range, style = "quantile", dim = 3)
bi_map <- 
  cartographeR(dt=bi_data, map_varname = 'bi_class', map_label = 'Rank vs Uncertainty',
               map_title = '',
               tag = 'gsa_results',
               scale_type='bivar', 
               get_plot = T)

bi_legend <- bi_legend(pal = "GrPink",
                    dim = 3,
                    xlab = "Baseline EHD Rank",
                    ylab = "Uncertainty Interval",
                    size = 8)

finalPlot <- ggdraw() +
  draw_plot(bi_map, 0, 0, 1, 1) +
  draw_plot(bi_legend, 0.01, .01, 0.2, 0.2)

file.path(viz.dir, 'figure_2c.png') %>% ggsave(finalPlot, height=8, width=12)

#plot the impact accuracy %
cartographeR(dt=brank_stats_dt, map_varname = 'Accuracy', map_label = 'Accuracy \n(Most Impacted)',
             map_title = '',
             tag = 'gsa_results',
             scale_type='cont')

# Figure 3 ----------------------------------------------------------
#plot uncertainty against the OG ranking
plot <-
ggplot(rank_stats_dt, aes(x=Nominal_rank, range, color=impacted_hierarchy %>% as.factor) ) +
  geom_jitter(width=.5, height=.5) +
  scale_color_brewer('Impact Status', palette = 'Paired') +
  scale_x_continuous('EHD Ranking') +
  scale_y_continuous('Rank Uncertainty') +
  theme_bw() 

file.path(viz.dir, 'gsa_rank_range_scatters.png') %>% ggsave(height=8, width=12)

plot <-
  ggplot(rank_stats_dt, aes(x=Nominal_rank, y=range, fill=impacted_hierarchy %>% as.factor)) +
  geom_hex(aes(alpha=log(..count..)), bins=50) +
  geom_vline(xintercept=1200, linetype='dashed', color='grey') +
  #scale_fill_brewer('Impact Status', palette = 'Paired') +
  #scale_fill_viridis() +
  scale_fill_manual('Impact\nAgreement', values=viridis::turbo(10)[c(1,7,8,10)]) +
  scale_x_continuous('Baseline EHD Ranking') +
  scale_y_continuous('Rank Uncertainty') +
  scale_alpha_continuous(guide='none', range=c(.5,1)) +
  theme_minimal() 

file.path(viz.dir, 'gsa_rank_range_hex_impacted.png') %>% ggsave(height=8, width=12)

plot <-
  ggplot(rank_stats_dt, aes(x=Nominal_rank, y=range, fill=Nominal_bin %>% as.factor)) +
  geom_hex(aes(alpha=log(..count..)), bins=50) +
  #scale_fill_brewer('Impact Status', palette = 'Paired') +
  scale_fill_viridis_d("EHD \nDeciles", option='turbo') +
  scale_x_continuous('Baseline EHD Ranking') +
  scale_y_continuous('Rank Uncertainty') +
  scale_alpha_continuous(guide='none', range=c(.5,1)) +
  theme_minimal() 

file.path(viz.dir, 'gsa_rank_range_hex.png') %>% ggsave(height=8, width=12)

plot <-
  ggplot(rank_stats_dt, aes(x=Nominal_rank, y=deviation, fill=Nominal_bin %>% as.factor)) +
  geom_hex(aes(alpha=log(..count..)), bins=40) +
  #scale_fill_brewer('Impact Status', palette = 'Paired') +
  scale_fill_viridis_d("EHD \nDeciles", option='turbo') +
  scale_x_continuous('Baseline EHD Ranking') +
  scale_y_continuous('Average Deviation') +
  scale_alpha_continuous(guide='none', range=c(.5,1)) +
  theme_minimal() 


file.path(viz.dir, 'gsa_rank_deviation_hex.png') %>% ggsave(height=8, width=12)

# Figure 4 ----------------------------------------------------------


#combine with accuracy and make the same plots
sens_stats_dt <- list(
  SA_res$Sensitivity %>% 
    .[, target := 'MARC'],
  SA_res$Accuracy %>% 
    .[, target := 'Impact \nAccuracy']
) %>% rbindlist %>% 
  .[, label := str_replace_all(parameters, '_', '\n')] %>% 
  .[label=='Norm', label := 'Normalization \nMethod']


# make stacked bar plot to show the first and total indices
ggplot(sens_stats_dt[!(parameters %like% 'Classification' | target %like% 'Impact' | sensitivity %like% 'j')], 
       aes(x=fct_reorder(label, original), y=original, color=target, group=target,
           ymax = high.ci, ymin = low.ci)) +
  geom_point(position=position_dodge(width = .9), size = 1.5,  shape=3) +
  geom_errorbar(position=position_dodge(width = .9), width = 0.2) +
  geom_hline(yintercept=0, linetype='dotted') +
  labs(
    x = NULL,
    y = NULL,
    fill = NULL) +
  scale_color_viridis_d(option='turbo') +
  facet_wrap(~sensitivity)+
  coord_flip()+
  theme_minimal() +

file.path(viz.dir, 'gsa_barplot.png') %>% ggsave(height=8, width=12)

ggplot(sens_stats_dt[!(sensitivity=='Sij')], 
       aes(x=fct_reorder(label, original), y=original, color=target, group=target,
           ymax = high.ci, ymin = low.ci)) +
  geom_point(position=position_dodge(width = .9), size = 1.5,  shape=3) +
  geom_errorbar(position=position_dodge(width = .9), width = 0.2) +
  geom_hline(yintercept=0, linetype='dotted') +
  labs(
    x = NULL,
    y = NULL,
    fill = NULL) +
  scale_color_brewer('Target Metric', palette='Dark2') + 
  facet_wrap(~sensitivity)+
  coord_flip()+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

file.path(viz.dir, 'gsa_double_barplot.png') %>% ggsave(height=8, width=12)

# make stacked bar plot to show the first and total indices
#caps at 0
SA_res$Accuracy[original<0, original:=0]
SA_res$Accuracy[low.ci<0, low.ci:=0]
ggplot(SA_res$Accuracy[!(sensitivity=='Sij')], aes(y=original, x=parameters,
                                                  ymax = high.ci, ymin = low.ci)) +
  geom_point(size = 1.5) +
  geom_errorbar(width = 0.2) +
  scale_y_continuous(limits=c(0,1.5)) +
  labs(
    x = NULL,
    y = NULL,
    fill = NULL) +
  facet_wrap(~sensitivity)+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

file.path(viz.dir, 'gsa_double_barplot.png') %>% ggsave(height=8, width=12)

#make a heatmap type chart to show the first order interactions
plot_dt <- SA_res$Sensitivity[sensitivity=='Sij'& !(parameters %like% 'Classification')] %>%
  copy %>%
  .[, c('var1', 'var2') := tstrsplit(parameters, split='.', fixed=T)] %>%
  .[, .(var1, var2, original)]

plot_dt <-
  SA_res$Sensitivity[sensitivity=='Sij'  & !(parameters %like% 'Classification')] %>% 
  copy %>% 
  .[, c('var2', 'var1') := tstrsplit(parameters, split='.', fixed=T)] %>% 
  .[, .(var1, var2, original)] %>% 
  list(.,
       plot_dt,
       SA_res$Sensitivity[sensitivity=='Si' & !(parameters %like% 'Classification'), 
                          .(var1=parameters, var2=parameters, original)]) %>% 
  rbindlist(use.names=T) %>% 
  .[as.integer(var1 %>% as.factor) < as.integer(var2 %>% as.factor), original := NA]

ggplot(plot_dt, aes(var1, var2, fill=original)) + 
  geom_tile() +
  scale_fill_viridis('Sensitivity Index', discrete=FALSE, na.value='white') +
  theme_minimal()

file.path(viz.dir, 'gsa_interactions_plot.png') %>% ggsave(height=8, width=12)

#make a heatmap type chart to show the first order interactions
plot_dt <- SA_res$Accuracy[sensitivity=='Sij']%>%
  copy %>%
  .[, c('var1', 'var2') := tstrsplit(parameters, split='.', fixed=T)] %>%
  .[, .(var1, var2, original)]

plot_dt <-
  SA_res$Accuracy[sensitivity=='Sij']%>% 
  copy %>% 
  .[, c('var2', 'var1') := tstrsplit(parameters, split='.', fixed=T)] %>% 
  .[, .(var1, var2, original)] %>% 
  list(.,
       plot_dt,
       SA_res$Accuracy[sensitivity=='Si', 
                          .(var1=parameters, var2=parameters, original)]) %>% 
  rbindlist(use.names=T) %>% 
  .[as.integer(var1 %>% as.factor) < as.integer(var2 %>% as.factor), original := NA]

ggplot(plot_dt, aes(var1, var2, fill=original)) + 
  geom_tile() +
  scale_fill_viridis('Sensitivity Index', discrete=FALSE, na.value='white') +
  theme_minimal()

file.path(viz.dir, 'gsa_acc_interactions_plot.png') %>% ggsave(height=8, width=12)

# Figure 5 ----------------------------------------------------------

#plot the MARCs by choice
#TODO needs more work
SA_res$diffs_dt[param!='Classification', .(MARC=mean(average_diff %>% abs)), by=.(param, sample)] %>% 
  .[order(MARC)] %>% 
  ggplot(., aes(forcats::fct_reorder(sample %>% as.factor, MARC), MARC, fill=MARC)) + 
  geom_bar(stat='identity') + 
  scale_x_discrete('') +
  scale_y_continuous('Average Rank Deviation') + 
  scale_fill_viridis() +
  facet_grid(~param) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

marc_dist_dt <-
SA_res$diffs_dt[param!='Classification'] %>% 
  .[, param := str_replace_all(param, '_', ' ')] %>% 
  merge(gsa_sample_map, by='sample')

marc_dist_dt %>% 
ggplot(aes(forcats::fct_reorder(label, average_diff), average_diff, fill=param)) + 
  #geom_violin(alpha=.4) + 
  #stat_summary(aes(color=param), fun.data = "mean_cl_boot", geom = "pointrange") +
  geom_boxplot(width = .2, outlier.shape = NA, coef = 1.5, notch = T) +
  scale_x_discrete('') +
  scale_y_continuous('Mean Absolute Rank Change') + 
  scale_color_manual('Parameter\nType', values=viridis::turbo(10)[c(1,2,4,7,9)]) +
  scale_fill_manual('Parameter\nType', values=viridis::turbo(10)[c(1,2,4,7,9)]) +
  coord_flip() +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

file.path(viz.dir, 'marc_distributions.png') %>% ggsave(height=8, width=12)

# Figure 7 (SI) ----------------------------------------------------------
#examine default plots
# plot_sensitivity(SA_res, ptype = "box")
# file.path(viz.dir, 'gsa_sensitivity_plot.png') %>% ggsave(height=8, width=12)
plot_uncertainty(SA_res, 'RankStats')
file.path(viz.dir, 'gsa_uncertainty_plot.png') %>% ggsave(height=8, width=12)
plot_uncertainty(SA_res, 'RankStats',
                 plot_units = dt[county_name=='King County', unique(GEOID)])
file.path(viz.dir, 'gsa_uncertainty_plot_king.png') %>% ggsave(height=8, width=12)
#***********************************************************************************************************************
 
# ---PLOT---------------------------------------------------------------------------------------------------------------
##other plots for the technical report and exploring some of the changes##
#reshape wide in order to run cor
corr_dt <- measure_ranks[, .(GEOID, item=paste0(theme, ': ', item), measure_new)] %>% 
  dcast(GEOID~item, value.var=c('measure_new')) %>% 
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
ggplot(le_reg_dt, aes(x=pollution, y=cjest_le, color=ruca_level) ) +
  geom_point(position='jitter') +
  #geom_hex(bins = 70) +
  scale_color_brewer('RUCA', palette='Set1') +
  theme_bw() 


#scatter the measure ranks to compare V1:v2
plot_dt <- merge(dt[level==3],
                 dt[level==1, .(GEOID, overall_shift=rank_shift)], 
                 by='GEOID')
ggplot(plot_dt, aes(x=rank_shift, overall_shift, color=impacted_hierarchy) ) +
  geom_point(position='jitter') +
  #geom_hex(bins = 70) +
  facet_wrap(~item_short) + 
  scale_x_continuous("Change in Indicator") +
  scale_y_continuous("Change in Overall Rank") +
  scale_color_brewer('Highly Impacted: Status', palette='Paired') +
  theme_minimal() 
#save the plot
file.path(viz.dir, 'measure_shift_scatters.png') %>% ggsave(height=8, width=12)

#scatter the measure ranks to compare V1:v2
ggplot(dt[level==3], aes(x=rank_v1, rank, color=impacted_hierarchy) ) +
  geom_point(position='jitter') +
  #geom_hex(bins = 70) +
  facet_wrap(~item_short) + 
  scale_x_continuous("Indicator Ranking, V1.1") +
  scale_y_continuous("Indicator Ranking, V2.0") +
  scale_color_brewer('Highly Impacted Status', palette='Paired') +
  theme_minimal() 
#save the plot
file.path(viz.dir, 'measure_scatters.png') %>% ggsave(height=8, width=12)

#county level plots
most_pop <- c('King', 'Pierce', 'Snohomish', 'Spokane', 'Clark', 'Thurston', 'Kitsap', 'Yakima', 'Whatcom', 'Benton',
              'Skagit', 'Cowlitz', 'Grant', 'Franklin', 'Island')
plot_dt <- merge(dt,
                 le_dt[, .(county, GEOID=geocode)],
                 by='GEOID')

#all counties-overall
ggplot(plot_dt[level==1], aes(x=rank, y=forcats::fct_reorder(county, rank, .fun=mean),
                                                     fill = stat(x))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Ranking", direction = -1, option='plasma') +
  scale_y_discrete('Counties') +
  theme_bw()
file.path(viz.dir, 'overall_ridges_counties.png') %>% ggsave(height=8, width=12)

#top 10 most populous counties
ggplot(plot_dt[level==1 & county %in% most_pop], aes(x=rank, y=forcats::fct_reorder(county, rank, .fun=mean),
                                        fill = stat(x))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Ranking", direction = -1, option='plasma') +
  scale_y_discrete('Counties (Top 10 by Population)') +
  theme_bw()

  file.path(viz.dir, 'overall_ridges_most_pop_counties.png') %>% ggsave(height=8, width=12)
  
ggplot(plot_dt[level==3 & county %in% most_pop], aes(x=rank, y=forcats::fct_reorder(county, index, .fun=mean),
                                                     fill = stat(x))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Ranking", direction = -1, option='plasma') +
  facet_wrap(~item_short) +
  geom_vline(xintercept = 8.5, linetype='dashed', color='dark red') +
  scale_y_discrete('Counties (Top 10 by Population)') +
  theme_minimal()
file.path(viz.dir, 'item_ridges_most_pop_counties.png') %>% ggsave(height=8, width=12)

#compare with other gov indicators
dt[, eji_compare:= round(rank-(eji_rpl*10), 0)]
dt[eji_compare<-5, eji_compare:=-5]
dt[eji_compare>5, eji_compare:=5]
dt[, eji_impacted_num := ((eji_rpl>.75) %>% as.numeric)]
dt[, cjest_impacted_num := (cjest_impacted %>% as.numeric)]
dt[, impacted_num := cjest_impacted]

#calculate metrics of agreement across the different 
dt[0==(is.na(eji_impacted_num)+is.na(cjest_impacted_num)), 
               tool_agreement := rgb(red=eji_impacted_num*179, 
                                     blue=cjest_impacted_num*175, 
                                     green=impacted*171,  maxColorValue = 255)]
#swap black and white, set missing to gray
dt[(eji_impacted_num+cjest_impacted_num+impacted)==0, tool_agreement := '#FFFFFF']
dt[(eji_impacted_num+cjest_impacted_num+impacted)==3, tool_agreement:='#000000']
dt[tool_agreement %>% is.na, tool_agreement:='#B3ABAF']
agreement_colors <- c('#FFFFFF',
                      '#0000AF',
                      '#00AB00',
                      '#B30000',
                      '#00ABAF',
                      '#B300AF',
                      '#B3AB00',
                      '#000000',
                      '#B3ABAF')
names(agreement_colors) <- c('None',
                             'CJEST',
                             'EHD',
                             'EJI',
                             'CJEST+EHD',
                             'CJEST+EJI',
                             'EJI+EHD',
                             'All',
                             'Missing')

#compare with the eji (after doing some rounding)
cartographeR(dt=dt, map_varname = 'eji_compare', map_label = 'EHD - EJI',
             scale_type='div_man', scale_vals=div_colors)

cartographeR(dt=dt, map_varname = 'tool_agreement', map_label = 'Cross-tool agreement',
             map_title = '',
             scale_type='identity', scale_vals=agreement_colors)

#graph LE at tract level
cartographeR(dt=dt, map_varname = 'cjest_le', map_label = 'Tract Level LE Estimates',
             map_title = '',
             scale_type='cont')

#compare ruralities
#create a manual diverging color scale to make sure that the index shifts are uniformly depicted
ruca_colors <- viridis_pal()(4)
names(ruca_colors) <- levels(dt$ruca_level)
cartographeR(dt=dt, map_varname = 'ruca_level', map_label = 'RUCA Classifications',
             map_title = '',
             tag = 'gsa_results',
             scale_type='cont_man', scale_vals=ruca_colors)
#***********************************************************************************************************************
 
# ---NORMS--------------------------------------------------------------------------------------------------------------
##transformations##
#setup some data with different types of transformation
#build coins
coin_zscore <- Treat(coin, dset = "Raw")
# coin_zscore <- Normalise(coin_zscore, dset = "Treated",
#                          global_specs=list(f_n='n_zscore'), f_n_para = list(m_sd=c(5,1)))

coin_zscore <- Normalise(coin_zscore, dset = "Treated",
          global_specs = list(f_n = "n_zscore",
                              f_n_para = list(c(5,1))))

coin_minmax <- Treat(coin, dset = "Raw")
coin_minmax <- Normalise(coin_minmax, dset = "Treated",
                         global_specs=list(f_n='n_minmax', f_n_para = list(c(0.001,10))))
coin_prank <- Normalise(coin, dset = "Raw",
                        global_specs=list(f_n='n_prank_log'))

#aggregate using base function to compare to our results
coin_minmax <- Aggregate(coin_minmax, dset = "Normalised", 
                         w='Original',
                         flatten_hierarchy = F,
                         #note that we take the arith mean until the last aggregation, where we multiply for risk score
                         f_ag = c("a_amean", "a_amean", "prod")) 

coin_zscore <- Aggregate(coin_zscore, dset = "Normalised", 
                         w='Original',
                         flatten_hierarchy = F,
                         #note that we take the arith mean until the last aggregation, where we multiply for risk score
                         f_ag = c("a_amean", "a_amean", "prod")) 



#prep a theme dt for merging to the collapsed item short for labelling and categorization
theme_labels_dt <-
  dt[level==3, .(indicator=item_short, theme)] %>% 
  .[, indicator := str_to_lower(indicator) %>% #cleanup so we can reshape wide
      str_replace_all(., ' \\(%\\)', '') %>% 
      str_replace_all(., ' ', '_') %>% 
      str_replace_all(., '\\(rsei\\)', 'rsei') %>%
      str_replace_all(., '2.5', '25')] %>% 
  unique(by='indicator')

#custom function to make the transformed dts
prepFx <- function(obj, label) {
  message('prepping', label, ' data')
  out <- obj %>% 
    as.data.table %>% 
    melt(id.vars=c('uCode'), value.name='value', variable.name='indicator') %>% 
    .[, type := label] %>% 
    setkey(uCode, indicator)
}

#merge all types together
# trans_dt <- list(
#   prepFx(coin$Data$Raw, 'raw'),
#   prepFx(coin$Data$Normalised, 'decile'),
#   prepFx(coin_prank$Data$Normalised, 'centile'),
#   prepFx(coin_zscore$Data$Normalised, 'zscore'),
#   prepFx(coin_minmax$Data$Normalised, 'minmax')
# ) %>%   
#   rbindlist %>% 
#   merge(theme_labels_dt, by='indicator')

#merge all types together
trans_dt <- prepFx(coin$Data$Raw, 'raw') %>% 
  .[, type := NULL] %>% 
  .[, minmax := n_minmax(value, c(0,10)), by=indicator] %>% 
  .[, centile := n_prank(value)*10, by=indicator] %>% 
  .[, decile := n_brank(value), by=indicator] %>% 
  .[, zscore := n_zscore(value, m_sd=c(5,2.5)), by=indicator] %>% 
  #.[, zscore := n_zscore(value), by=indicator] %>% 
  setnames('value', 'raw') %>% 
  melt(id.var=c('uCode', 'indicator'), value.name='value', variable.name='type') %>% 
  merge(theme_labels_dt, by='indicator')


#scale centile by 10 so it's easier to plot in the same space
#trans_dt[type=='centile', value := value*10]
#trans_dt[type=='zscore', value := value %>% n_minmax(l_u=c(0.001, 10))] #also scale the Zs

#make a more readable label
trans_dt[, indicator_label := str_replace_all(indicator, '_', '\n')]

#plot the raw data distributions
ggplot(trans_dt[type=='raw'], aes(value, fill=theme)) +
  geom_density(alpha=.5) +
  scale_fill_manual('Themes', values=viridis::turbo(10)[c(1,2,7,10)]) +
  scale_y_continuous('') +
  scale_x_continuous('') +
  facet_wrap(~indicator_label, scales = 'free') +
  theme_minimal()
file.path(viz.dir, 'raw_distributions.png') %>% ggsave(height=8, width=12)

#plot the raw data distributions
ggplot(trans_dt[type=='raw' & theme%like%'Exposure'], aes(value, fill=theme)) +
  geom_density(alpha=.5) +
  scale_fill_manual('Themes', values=viridis::turbo(10)[c(1,2,7,10)]) +
  scale_y_continuous('') +
  scale_x_continuous('') +
  facet_wrap(~indicator_label, scales = 'free') +
  theme_minimal()
file.path(viz.dir, 'raw_distributions_exp.png') %>% ggsave(height=8, width=12)

themDist <- function(dt, this_theme) {
  
  plot <-
    ggplot(dt[theme==this_theme & type %in% c('minmax', 'zscore', 'centile', 'decile')], 
           aes(value, fill=type)) +
    geom_density(alpha=.5) +
    scale_fill_viridis_d('Transformations', option='magma') +
    scale_y_sqrt('') +
    scale_x_continuous('') +
    facet_wrap(~indicator_label) +
    theme_minimal()
  
  plot
  
}

pdf(file.path(viz.dir, 'trans_distributions.pdf'))
lapply(unique(trans_dt$theme), themDist, dt=trans_dt)
dev.off()

ggplot(trans_dt[!(type=='raw')], aes(value, fill=type)) +
  geom_density(alpha=.5) +
  scale_fill_manual('Method', values=viridis::turbo(10)[c(1,2,7,10)]) +
  scale_y_continuous('') +
  scale_x_continuous('', limits=c(0,10)) +
  facet_wrap(~indicator_label, scales='free') +
  theme_minimal()
file.path(viz.dir, 'trans_distributions.png') %>% ggsave(height=8, width=12)

ggplot(trans_dt[!(type=='raw') & theme %like% 'Exposures'], aes(value, fill=type)) +
  geom_density(alpha=.5) +
  scale_fill_manual('Method', values=viridis::turbo(10)[c(1,2,7,10)]) +
  scale_y_continuous('') +
  scale_x_continuous('', limits=c(0,10)) +
  facet_wrap(~indicator_label, scales='free') +
  theme_minimal()
file.path(viz.dir, 'trans_distributions_exp.png') %>% ggsave(height=8, width=12)
#***********************************************************************************************************************

# ---COMPRESSION--------------------------------------------------------------------------------------------------------
#make a wide version too to calc the residuals
zscore_p80 <-  0.8416 #value of zscore for 80th pctile
trans_wide_dt <- dcast(trans_dt,
                       ...~type, 
                       value.var='value') %>% 
  #.[, zz_rank := n_zscore(zscore_rank)] %>% 
  #.[, ze_rank := n_zscore(ehd_rank)] %>% 
  .[, cent_zscore := n_zscore(centile), by=indicator] %>% 
  .[, z_med := quantile(zscore, p=.5, na.rm=T), by=indicator] %>% 
  .[, z_p25 := quantile(zscore, p=.25, na.rm=T), by=indicator] %>% 
  .[, z_p75 := quantile(zscore, p=.75, na.rm=T), by=indicator] %>% 
  .[, GEOID := uCode] %>% 
  .[, level :=1]

#fix theme names
trans_wide_dt[theme=='Environmental Effects', theme_short := 'env fx: ']
trans_wide_dt[theme=='Environmental Exposures', theme_short := 'env exp: ']
trans_wide_dt[theme=='Socioeconomic Factors', theme_short := 'socio: ']
trans_wide_dt[theme=='Sensitive Populations', theme_short := 'sensitivity: ']

#calculate the residuals by GEOID
trans_wide_dt[, minmax_resid := centile-minmax]
trans_wide_dt[, zscore_resid := centile-zscore]

agg_dt <-
merge(coin_zscore$Data$Aggregated %>% 
       as.data.table %>% 
       .[, .(uCode, zscore_rank=n_prank(ehd_rank), theme_short='EHD: ', indicator='composite')],
     coin$Data$Aggregated %>% 
       as.data.table %>% 
       .[, .(uCode, cent_rank=n_prank(ehd_rank), theme_short='EHD: ', indicator='composite')],
     by=c('uCode', 'theme_short', 'indicator')
     ) %>% 
  .[, rank_diff := cent_rank-zscore_rank] %>% 
  .[, zscore := zscore_rank %>% n_zscore] %>% 
  .[, cent_zscore := cent_rank %>% n_zscore] %>% 
  .[, classification := 'Unimpacted'] %>% 
  .[cent_rank < .8 & zscore_rank > .8, classification := 'Impacted (z-scores)'] %>% 
  .[cent_rank > .8 & zscore_rank < .8, classification := 'Impacted (EHD)'] %>% 
  .[cent_rank > .8 & zscore_rank > .8, classification := 'Impacted (both)'] 

#find cases and label
#find one each of
#geoid that is not impacted on ehd and impacted on minmax
case1_id <- agg_dt[rank_diff==agg_dt[classification%like%'z-score', min(rank_diff, na.rm=T)], uCode]
#geoid that is impacted on ehd and not impacted on minmax
case2_id <- agg_dt[rank_diff==agg_dt[classification%like%'EHD', max(rank_diff, na.rm=T)], uCode]

#bind all
agg_dt <- list(
  agg_dt,
  trans_wide_dt
) %>% rbindlist(use.names=T, fill=T) %>% 
  .[, ind := paste0(theme_short, indicator %>% str_replace_all('_', ' '))]

#label the cases
agg_dt[uCode==case1_id, case := 1]
agg_dt[uCode==case2_id, case := 2]
agg_dt[case==1, case_name := 'Okanogan']
agg_dt[case==2, case_name := 'Greenwood, Seattle']

#keep only cases
case_dt <- na.omit(agg_dt, cols='case')  %>% 
  .[, level:=1] %>% 
  .[, GEOID := uCode]

#merge dropouts to shapefile and plot
shp <- tract_sf %>% 
  merge(case_dt[, .(GEOID, case, case_name)] %>% unique, by='GEOID') 

#graph residuals at tract level
cartographeR(dt=case_dt, map_varname = 'case', map_label = 'Cases',
             map_title = '',
             scale_type='cont')

# Create a color palette for the map:
mypalette <- colorNumeric( palette="viridis", domain=shp$case, na.color="transparent")
mypalette(c(1,2, NA))

# Basic choropleth with leaflet?
m <- leaflet(shp) %>% 
  addTiles()  %>% 
  #setView( lat=10, lng=0 , zoom=2) %>%
  addPolygons( color = ~mypalette(case), stroke=FALSE )

leaflet(shp) %>% 
  addProviderTiles("CartoDB.Positron") %>% 
  addPolygons(group='geoid', color = "green")

#plot the raw data distributions
trans_wide_dt %>% 
  merge(agg_dt[, .(uCode, classification)] %>% na.omit, by='uCode') %>% 
  .[classification!='Impacted (both)'] %>% 
  .[theme_short %like% 'soc'] %>% 
ggplot(aes(raw, fill=classification)) +
  geom_density(alpha=.5) +
  scale_fill_manual('Classification', values=c('#984ea3', '#ff7f00', '#377eb8')) +
  scale_y_continuous('') +
  scale_x_continuous('') +
  facet_wrap(~indicator_label, scales = 'free') +
  theme_minimal()
file.path(viz.dir, 'comparing_distributions.png') %>% ggsave(height=8, width=12)

#collapse the residuals to GEOD
resid_dt <- trans_wide_dt %>% 
  .[, avg_zscore_resid := mean(zscore_resid %>% abs, na.rm=T), by=.(uCode)] %>% 
  .[, avg_minmax_resid := mean(minmax_resid %>% abs, na.rm=T), by=.(uCode)] %>% 
  unique(by='uCode') %>% 
  .[, level := 1] %>% 
  .[, GEOID := uCode] %>% 
  merge(rank_stats_dt[, .(GEOID, Nominal_rank, ehd_rank, impacted_hierarchy)], by='GEOID')

#graph residuals at tract level
cartographeR(dt=resid_dt, map_varname = 'minmax_resid', map_label = 'Transformation Residuals',
             map_title = '',
             scale_type='cont')

#graph difference in rank at tract level
cartographeR(dt=resid_dt, map_varname = 'rank_diff', map_label = 'Change between base/minmax',
             map_title = '',
             scale_type='cont_grad')

#graph difference in rank at tract level
cartographeR(dt=agg_dt[indicator%like%'composite'], map_varname = 'classification', 
             map_label = 'Classification',
             map_title = '',
             scale_type='class')

resid_dt$ehd_bin <- cut(resid_dt$ehd_rank, breaks = c(0,2,8,10), include.lowest = TRUE)
bi_data <- bi_class(resid_dt, x = ehd_bin, y = avg_zscore_resid, style = "quantile", dim = 3)
bi_map <- 
  cartographeR(dt=bi_data, map_varname = 'bi_class', map_label = 'Rank vs Residual',
               map_title = '',
               tag = 'gsa_results',
               scale_type='bivar', 
               get_plot = T)

bi_legend <- bi_legend(pal = "GrPink",
                       dim = 3,
                       xlab = "Baseline EHD Rank",
                       ylab = "Residual",
                       size = 8)

finalPlot <- ggdraw() +
  draw_plot(bi_map, 0, 0, 1, 1) +
  draw_plot(bi_legend, 0.01, .01, 0.2, 0.2)

file.path(viz.dir, 'resid_bimap.png') %>% ggsave(finalPlot, height=8, width=12)

#take a look at the subset that have mid EHD and high residuals
resid_ids <- bi_data[bi_class=='2-3', unique(GEOID)]
trans_wide_dt[uCode%in%resid_ids]


#graph difference in rank at tract level
  case_map <- 
  cartographeR(dt=case_dt, map_varname = 'case', map_label = 'Case studies',
               map_title = '',
               scale_type='cont_vir',
               get_plot = T)

  casePlot <- function(this_case, dt)  {
    case_dt[case %in% this_case] %>% 
      .[, .(uCode, case, x=ind %>% fct_rev,
            theme, raw, 
            value1=cent_zscore, value2=zscore, 
            p25=z_p25, med=z_med, p75=z_p75)] %>% 
      ggplot(aes(x=x, shape=theme)) +
      geom_point( aes(x=x, y=med), color='#ff7f00', size=2, shape=3 ) +
      geom_segment( aes(x=x, xend=x, y=p25, yend=p75), color='#ff7f00') +
      geom_segment( aes(x=x, xend=x, y=med, yend=value2), color='#ff7f00', linetype='dotted') +
      geom_hline(aes(yintercept=0), color='black') +
      geom_point( aes(x=x, y=value1), color='#984ea3', alpha=1, size=3) +
      geom_point( aes(x=x, y=value2), color='#ff7f00', alpha=1, size=3) +
      geom_hline(aes(yintercept=zscore_p80), color='#377eb8', size=1) +
      #geom_point( aes(x=x, y=zscore), color='#e41a1c', size=3, alpha=.6) +
      scale_color_manual('Themes', values=viridis::turbo(10)[c(1,2,7,10)]) +
      scale_shape_manual(guide='none', values=c(19,19,19,19), na.value = 17) +
      scale_x_discrete(position='top') +
      scale_y_continuous(limits=c(-5,5)) +
      coord_flip()+
      facet_wrap(~case, ncol=1, labeller = cust_label %>% as_labeller) +
      theme_minimal() +
      theme(
        legend.position = "none",
        strip.text = element_text(size = 14)
      ) +
      #ggtitle('Scaling Compression Effects') +
      xlab("") +
      ylab("")
  }
  
  casePlot(1:2, case_dt)
  file.path(viz.dir, 'scaling_compression_effect_zspace.png') %>% ggsave(height=8, width=12)


pdf(file.path(viz.dir, 'scaling_compression_effect.pdf'))
lapply(1:4, casePlot, dt=case_dt)
dev.off()

file.path(viz.dir, 'scaling_compression_effect.png') %>% ggsave(height=8, width=12)


#***********************************************************************************************************************

# ---ZSCORES------------------------------------------------------------------------------------------------------------
#move the caseplot into zscore space
casePlot2 <- function(this_case, dt)  {
  case_dt[case %in% this_case & !is.na(centile)] %>% 
    .[, .(uCode, case, x=ind, z_med, 
          theme, raw, 
          #value1_rel=abs(centile-.5), value2_rel=abs(minmax-median(minmax, na.rm=T)), 
          value1=cent_zscore, value2=zscore, 
          p25=z_p25, med=z_med, p75=z_p75)] %>% 
    ggplot(aes(x=x)) +
    geom_point( aes(x=x, y=med), color='#ff7f00', size=2, shape=3 ) +
    geom_segment( aes(x=x, xend=x, y=p25, yend=p75), color='#ff7f00') +
    geom_segment( aes(x=x, xend=x, y=med, yend=value2), color='#ff7f00', linetype='dotted') +
    #geom_point( aes(x=x, y=5), color='#984ea3', size=2, shape=3 ) +
    #geom_segment( aes(x=x, xend=x, y=5, yend=value1), color='#984ea3') +
    #geom_segment( aes(x=x, xend=x, y=value1, yend=value2, color=theme)) +
    #geom_hline(aes(yintercept=ze_rank), color='#984ea3') +
    #geom_hline(aes(yintercept=zscore_rank*10), color='#e41a1c') +
    #geom_hline(aes(yintercept=zz_rank), color='#ff7f00') +
    geom_hline(aes(yintercept=0), color='black') +
    geom_hline(aes(yintercept=zscore_p80), color='grey', linetype='dashed') +
    geom_point( aes(x=x, y=value1), color='#984ea3', alpha=.6, size=3) +
    geom_point( aes(x=x, y=value2), color='#ff7f00', alpha=.6, size=3) +
    #geom_point( aes(x=x, y=zscore), color='#e41a1c', size=3, alpha=.6) +
    scale_color_manual('Themes', values=viridis::turbo(10)[c(1,2,7,10)]) +
    scale_y_continuous(limits=c(-5,5)) +
    coord_flip()+
    #facet_wrap(~case_name, ncol=1) +
    theme_minimal() +
    theme(
      legend.position = "none",
    ) +
    ggtitle('Scaling Compression Effects') +
    xlab("") +
    ylab("Local z-score of value (Orange) vs. Local z-score of centile (Purple)")
}

casePlot2(c(2,4), case_dt)
file.path(viz.dir, 'scaling_compression_effect_zspace.png') %>% ggsave(height=8, width=12)

case_dt[case %in% c(2,4) & !is.na(centile)] %>% 
  .[, .(uCode, case, x=indicator_label, z_med, 
        theme, raw, 
        #value1_rel=abs(centile-.5), value2_rel=abs(minmax-median(minmax, na.rm=T)), 
        value1=cent_zscore, value2=zscore, 
        p25=z_p25, med=z_med, p75=z_p75,
        zscore, zscore_resid, ehd_rank, zscore_rank)] %>% 
  ggplot(aes(x=value1, y=value2, color=theme, label=x)) +
  geom_point(size=2) +
  #geom_text_repel() +
  geom_hline(aes(yintercept=ehd_rank*10), color='#984ea3') +
  #geom_hline(aes(yintercept=zscore_rank*10), color='#e41a1c') +
  geom_hline(aes(yintercept=zscore_rank*10), color='#ff7f00') +
  scale_color_manual('Themes', values=viridis::turbo(10)[c(1,2,7,10)]) +
  scale_x_continuous(limits=c(0,10)) +
  scale_y_continuous(limits=c(0,10)) +
  coord_flip()+
  facet_wrap(~case_name, nrow=1) +
  theme_minimal() +
  theme(
    legend.position = "none",
  ) +
  ggtitle('Scaling Compression Effects')

file.path(viz.dir, 'scaling_compression_effect_zspace_scatter.png') %>% ggsave(height=8, width=12)


#look at the relationship between the residuals and ehd score
ggplot(resid_dt, aes(x=Nominal_rank, y=avg_zscore_resid, fill=impacted_hierarchy %>% as.factor)) +
  geom_hex(aes(alpha=log(..count..)), bins=50) +
  #scale_fill_brewer('Impact Status', palette = 'Paired') +
  #scale_fill_viridis() +
  scale_fill_manual('Impact\nAgreement', values=viridis::turbo(10)[c(1,7,8,10)]) +
  scale_x_continuous('Baseline EHD Ranking') +
  scale_y_continuous('Rank Uncertainty') +
  scale_alpha_continuous(guide='none', range=c(.5,1)) +
  theme_minimal() 

file.path(viz.dir, 'gsa_rank_range_hex_impacted.png') %>% ggsave(height=8, width=12)

#plot residuals vs ehd score
trans_wide_dt[, ind_lab := str_replace_all(indicator, '_', ' ')]
ggplot(trans_wide_dt, aes(centile, raw, fill=theme)) +
  geom_hex(aes(alpha=log(..count..)), bins=150) +
  #geom_point(alpha=.6) +
  scale_fill_manual('Indicator', values=viridis::turbo(10)[c(1,7,8,10)]) +
  #scale_color_manual(guide='none', values=viridis::turbo(10)[c(1,2,7,10)]) +
  scale_color_viridis() +
  scale_alpha_continuous(guide='none', range=c(.5,1)) +
  scale_y_continuous('Z-score', limits=c(0,1)) +
  scale_x_continuous('Centiles') +
  #facet_wrap(~indicator_label, scales = 'free') +
  theme_minimal()

file.path(viz.dir, 'trans_scatters.png') %>% ggsave(height=8, width=12)

#***********************************************************************************************************************

# ---PM2.5--------------------------------------------------------------------------------------------------------------
##examine PM2.5 in risk space##
pm_brt_dt <- file.path(data.dir, 'cvd_ihd.csv') %>% fread
pm_dt <- dt[level==3 & item %like% 'PM 2.5 Concentration', .(GEOID, item, item_short, measure, measure_v1, cjest_le)] %>% 
  copy %>% 
  .[, exposure := measure_v1]

#merge on the CVD RR values
pm_dt <- pm_brt_dt[pm_dt, on=.(exposure), roll='nearest']

trans_dt <- list(trans_dt,
                 pm_dt[, .(indicator='pm25_concentration', value=mean, type='rr', uCode=GEOID, cjest_le)]) %>% 
  rbindlist(fill=T, use.names=T)

#***********************************************************************************************************************

# ---PCA----------------------------------------------------------------------------------------------------------------
##pca##
#also add on the PCA + individual raw data
pc_dt <- get_PCA(coin_minmax, dset='Normalised', by_groups=F, out2='preds', imputed=T) %>% 
  cbind(coin$Data$Raw, .) %>% 
  setnames(c('uCode', 'ehd_rank'), c('GEOID', 'ehd_pc')) %>% 
  .[, level := 1]

#extract the PCA model info as well for diagnostic plotting
pc_mod <- get_PCA(coin_minmax, dset='Normalised', by_groups=F, out2='mod', imputed=T)
pc_dt <- pc_mod$PCAresults$All$all_preds %>%
  .[, .(variable, id, value, combined)] %>% 
  dcast(...~variable, id.var='id', value.var='value') %>% 
  cbind(pc_dt) %>% 
  merge(dt[level==1, .(GEOID, overall=rank, ruca_pop, county_name)], by='GEOID')

pc_dt[, pca_deciles := n_brank(combined)]
pc_dt[, diff := overall-pca_deciles]
pc_dt[diff<(-5), diff := -5]
pc_dt[diff>5, diff := 5]

cartographeR(dt=pc_dt, map_varname = 'PC1', map_label = 'PC1',
             map_title = '',
             scale_type='cont')

cartographeR(dt=pc_dt, map_varname = 'PC2', map_label = 'PC2',
             map_title = '',
             scale_type='cont')

pc_mod$PCAresults$All$contrib %>% 
  merge(theme_labels_dt, by.x='name', by.y='indicator') %>% 
  ggplot(aes(x=fct_reorder(name, contrib*wt, .fun=sum), 
             y = contrib*wt, alpha=pc %>% fct_rev, fill=theme)) +
  geom_bar(stat='identity') +
  scale_fill_manual('Themes', values=viridis::turbo(10)[c(1,2,7,10)]) +
  scale_alpha_discrete('PCs', range=c(.25,1)) +
  scale_y_continuous('Weighted Contribution') +
  scale_x_discrete('') +
  coord_flip() +
  theme(legend.position = "top") +
  theme_minimal()
file.path(viz.dir, 'pca_contrib.png') %>% ggsave(height=8, width=12)

fviz_pca_var(pc_mod$PCAresults$All$PCAres,
             col.var = "contrib", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)   + ggtitle(paste0('PCA Loadings'))
file.path(viz.dir, 'pca_biplot.png') %>% ggsave(height=8, width=12)

plot <- 
fviz_pca_ind(pc_mod$PCAresults$All$PCAres,
             axes=c(1,2),
             pointshape = 19,
             label='none',
             col.var='black',
             habillage = dt[level==3 & item %like% 'Death', impacted_hierarchy], # Color by the quality of representation
             repel = TRUE     # Avoid text overlapping
) +
  ggtitle('PCA Individuals') +
  #geom_hex(data=dt[level==3 & item %like% 'Death'], aes(fill=impacted_hierarchy, alpha=log(..count..)), bins=35) +
  scale_color_manual('', values=viridis::turbo(10)[c(2,7,8,10)]) +
  theme_minimal()
file.path(viz.dir, 'pca_impacted_biplot.png') %>% ggsave(height=8, width=12)

plot$data %>% 
  na.omit() %>% 
ggplot(aes(x=x, y=y, fill=Groups)) +
   ggtitle('PCA Individuals') +
   geom_hex(aes(alpha=log(..count..)), bins=60) +
   geom_vline(xintercept=0) +
   geom_hline(yintercept = 0) +
   scale_fill_manual('', values=viridis::turbo(10)[c(2,7,8,10)]) +
   scale_alpha_continuous(guide='none', range=c(.6,1)) +
   scale_x_continuous('PC1', limits=c(-6, 6)) +
   scale_y_continuous('PC2', c(-6, 6)) +
   theme_minimal()

file.path(viz.dir, 'pca_impacted_biplot_hex.png') %>% ggsave(height=8, width=12)

#make a custom biplot
biplot_dt <-
pc_mod$PCAresults$All$PCAres$rotation %>% 
  as.data.table(keep.rownames = T) %>% 
  .[, .(indicator=rn, PC1=PC1, PC2=PC2)] %>% 
  merge(., trans_dt[, .(indicator, theme, indicator_label)] %>% unique, 
        by='indicator')
  
themedBiplot <- function(this_theme, dt) {
  
  message('labeling for ', this_theme)
  
  this_dt <- copy(dt) %>% 
    .[theme==this_theme, lab := indicator_label]
  
  plot <-
    ggplot(data=this_dt, aes(x=PC1, y=PC2, color=theme, label=lab)) +
    ggtitle('PCA Biplot') +
    geom_vline(xintercept=0) +
    geom_hline(yintercept = 0) +
    geom_point()+
    geom_segment(aes(x=0, y=0, xend=PC1, yend=PC2), linetype='dotted') +
    geom_label_repel(show.legend = F)+
    scale_color_manual('', values=viridis::turbo(10)[c(2,7,8,10)]) +
    scale_x_continuous('PC1 (28%)', limits=c(-.5, .5)) +
    scale_y_continuous('PC2 (18%)', limits=c(-.5, .5)) +
    theme_minimal()
  
  return(plot)

}

pdf(file.path(viz.dir, 'themed_biplots.pdf'), height = 8, width=12)
lapply(trans_dt$theme %>% unique, themedBiplot, dt=biplot_dt)
dev.off()

#see classification percentage
plot$data %>% as.data.table %>% .[, sum(y>=0&x>=0)/.N, by=Groups]

fviz_pca_ind(pc_mod$PCAresults$All$PCAres,
             pointshape = 19,
             label='none',
             col.ind= dt[level==3 & item %like% 'Death', cjest_le] %>% as.integer, # Color by the quality of representation
             repel = TRUE     # Avoid text overlapping
) +
  ggtitle('PCA Individuals') +
  scale_color_viridis_c('Life Expectancy', option='viridis', direction=-1) +
  scale_fill_viridis()+
  theme_minimal()

file.path(viz.dir, 'pca_impacted_le.png') %>% ggsave(height=8, width=12)

#graph aggregated PCA at tract level
cartographeR(dt=pc_dt, map_varname = 'ehd_pc', map_label = 'Combined\nPCA',
             map_title = '',
             scale_type='cont')

#first make a histogram of county level differences
diff_hist <-
  pc_dt[, weighted.mean(diff, ruca_pop), by=county_name][order(county_name)] %>% 
  .[, county_short := str_remove(county_name, ' County')] %>% 
  ggplot(., aes(forcats::fct_reorder(county_short %>% as.factor, V1), V1, fill=V1)) + 
  geom_bar(stat='identity') + 
  scale_x_discrete('') +
  scale_y_continuous('Average Difference in Deciles') + 
  scale_fill_gradient2(guide='none', na.value = "grey75", high=muted('red'), low=muted('blue')) +
  coord_flip() +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 35, vjust = 0.5, hjust=1))

file.path(viz.dir, 'pca_mean_difference_county.png') %>% ggsave(height=8, width=12)

#then map the tract level deviations
diff_map <-
  cartographeR(dt=pc_dt, map_varname = 'diff', map_label = 'Difference',
               map_title = '',
               scale_type='div_man', scale_vals = div_colors,
               get_plot=T)

all_grobs <- list(diff_map+
                    theme(plot.margin = unit(c(0,0,0,-1), units = "cm"),
                          legend.position = "none"), 
                  diff_hist)


plot <- arrangeGrob(grobs=all_grobs, layout_matrix=lay, 
                    top=textGrob("EHD Index Agreement with PCA Index", 
                                 gp = gpar(fontsize=17))
) %>% 
  grid.arrange

ggsave(plot=plot, filename=file.path(viz.dir, 'paper_2fig_5d.png'),
       width=12, height=8, units='in', dpi=900)

#***********************************************************************************************************************

# ---APP LE-------------------------------------------------------------------------------------------------------------
#generate a dataset to regress indicators on life expectancy
le_reg_dt <- dt[level==2, .(GEOID, theme, rank, cjest_le, ruca_level)] %>% 
  copy %>% 
  dcast(GEOID+cjest_le+ruca_level~theme, value.var=c('rank'), fun.aggregate = mean) %>% 
  na.omit %>% 
  setnames(c('Environmental Effects', 'Environmental Exposures', 'Sensitive Populations', 'Socioeconomic Factors'),
           c('Env_Effects', 'Env_Exposures', 'Sens_Pops', 'Soc_Factors'))

le_reg_dt[, pollution := (Env_Effects*.5+Env_Exposures)/2]
le_reg_dt[, pollution_scaled := scales::rescale(pollution)]
le_reg_dt[, pops := (Sens_Pops+Soc_Factors)/2]
le_reg_dt[, pops_scaled := scales::rescale(pops)]
le_reg_dt[, ehd_score := pollution*pops]

#also add on the IHME le data for just kingco
#ihme life expectancy data
regexp <- "[[:digit:]]+"
king_le_dt <-  file.path(data.dir, 'ihme_kingco_le_1990_2014.csv') %>%
  fread %>%
  .[year_id==2010 & sex_id==3, .(location_name, val, lower, upper)] %>%  #TODO re-evaluate year
  .[, tract_id := stringr::str_replace(location_name, 'King County Census Tract ', '')] %>%
  .[tract_id!='King County'] %>%  #drop the all county vals %>%
  .[, location_name := NULL] %>%
  setnames(., names(.), c('tract_le', 'tract_le_lower', 'tract_le_upper', 'tract_id'))

# #merge the GEOIDs onto the kingco le data
king_le_dt <- tract_sf %>%
  as.data.table %>%
  .[COUNTY=='033', .(NAME, GEOID)] %>%
  merge(., king_le_dt,
        by.x='NAME',
        by.y='tract_id')

le_reg_dt <- merge(le_reg_dt, king_le_dt, by='GEOID', all.x=T)

#graph LE at tract level
cartographeR(dt=pc_dt, map_varname = 'cjest_le', map_label = 'Life Expectancy',
             map_title = '',
             scale_type='cont')

#create plot of life expectancy distributions
# ggplot(index_dt, aes(rank %>% as.factor, le)) + 
#   geom_violin(aes(col = rank %>% as.factor, fill = rank %>% as.factor), alpha = 0.25)+
#   labs(x = "Index", y = "Life Expectancy",
#        title = "Relationship between EHD Rank and Life Expectancy",
#        subtitle = "Average Life Expectancy Between 2015-2019 (dotted line represents state average of 80.3 years)"
#   ) +
#   geom_vline(xintercept = 8.5, linetype='dashed', color='dark red') +
#   geom_hline(yintercept = index_dt[1, le_state_average], linetype='dotted', color='dark blue') +
#   scale_x_discrete() +
#   scale_color_manual('Rank', values=cont_colors, na.value = "grey75") +
#   scale_fill_manual('Rank', values=cont_colors, na.value = "grey75") +
#   geom_boxplot(color = "gray20", width = 0.15, coef = 1.5) +
#   annotate("text", x = 9.25, y = 99, label = "Highly Impacted Tracts", vjust = -0.5, color='dark red') +
#   theme_minimal() +
#   theme(legend.position = c(1, .99), legend.justification = c(1, 1),
#         plot.margin = unit(c(0, 0, 0, 0), "in"))

#models with LE
mod1 <- lm(cjest_le ~ ehd_score, data=le_reg_dt)
mod2 <- lm(cjest_le ~ ehd_pc, data=pc_dt)

stargazer(mod1, mod2, 
          type='html') %>%
  capture.output(file=file.path(out.dir, 'table_3.html'))

#make some hexplots looking at the relationship between trans and LE
themHex <- function(dt, this_theme) {
  plot <-
  ggplot(dt[theme==this_theme & type %in% c('minmax', 'zscore', 'centile', 'decile')], 
         aes(x=value, y=cjest_le)) +
    geom_hex(aes(fill=log(..count..), alpha=log(..count..)), bins=35) +
    geom_smooth(method='scam', 
                # b-spline monotonic deceasing
                # see ?shape.constrained.smooth.terms
                formula = y ~ s(x, k = 5, bs = "mpd"),
                color='black', linetype='dashed') +
    scale_fill_viridis(guide='none') +
    scale_y_continuous('Life Expectancy') +
    scale_x_continuous('Normalized Value', limits=c(0,10), breaks=c(1,5,9)) +
    scale_alpha_continuous(guide='none', range=c(.5,1)) +
    #facet_wrap(indicator~type, ncol=4, strip.position = 'right') +
    facet_grid(rows=vars(indicator_label), cols=vars(type), scales='free') +
    ggtitle('Indicator Relationship with Life Expectancy', subtitle = this_theme) +
    theme_minimal() +
    theme(strip.text.y.right = element_text(angle = 0))
  
  plot
  
}

pdf(file.path(viz.dir, 'transformations_v_le.pdf'))
lapply(unique(trans_dt$theme), themHex, dt=trans_dt)
dev.off()

##simulation study##

#bring in transformed dataset
sim_dt <- trans_wide_dt[, .(uCode, indicator, theme, indicator_label, 
                            raw, 
                            centile, decile, minmax, zscore)] 

#simulation function
runSims <- function(i, data_dt=sim_dt[theme%like%'Sens']) {

  message('rep ', i)
  
  #simulate some betas
  sim_betas <- 
    data.table(ind=unique(data_dt$indicator),
               beta=runif(uniqueN(data_dt$indicator), 0, .25))
  
  #make a response var from them
  this_sim <- data_dt[, .(uCode, indicator, minmax)] %>% 
    copy %>% 
    dcast(...~indicator, value.var = 'minmax') %>% 
    .[, sim_response := t(t(.[, sim_betas$ind, with=F])*sim_betas$beta) %>% rowSums(na.rm=T)]
  
  #merge back on
  data_dt <- merge(data_dt, 
                  this_sim[, .(uCode, sim_response)],
                  by='uCode')

  simByType <- function(type, dt) {
    
    lm(paste0('sim_response ~ indicator:', type) %>% as.formula, data=dt) %>% 
      tidy(conf.int = TRUE) %>% 
      as.data.table %>% 
      #.[, sim := sim_i] %>% 
      .[term!='(Intercept)'] %>% 
      .[, indicator := str_replace_all(term, 'indicator', '') %>% 
          str_replace_all(paste0(':', type), '')] %>% 
      .[, term := str_replace_all(indicator, '_', ' ')] %>% 
      merge(sim_betas, by.x='indicator', by.y='ind') %>% 
      .[, error:=estimate-beta] %>% 
      .[, abs_error:=abs(estimate-beta)] %>% 
      .[, coverage := (beta>conf.low & beta<conf.high)] %>% 
      .[, norm_type := type]
    
  }
  
  sim_results <- lapply(unique(trans_dt$type), simByType, dt=data_dt) %>% 
    rbindlist
  
  return(sim_results)

}

tmp <- mclapply(1:100, runSims, mc.cores=4) %>% 
  rbindlist

tmp[norm_type!='raw', .(mean_error=mean(error), mean_abs_error=mean(abs_error), coverage=sum(coverage)/.N),
    by=.(indicator, norm_type)] %>% 
  ggplot(aes(indicator, coverage, fill = norm_type)) +
  geom_bar(position = "dodge2", stat='identity') +
  coord_flip() + 
  theme_minimal()

#run models testing the relationships with 
mod1 <- lm(cjest_le ~ indicator:value, data=trans_dt[type=='raw'])
mod2 <- lm(cjest_le ~ indicator:value, data=trans_dt[type=='centile'])
mod3 <- lm(cjest_le ~ indicator:value, data=trans_dt[type=='decile'])
mod4 <- lm(cjest_le ~ indicator:value, data=trans_dt[type=='minmax'])
mod5 <- lm(cjest_le ~ indicator:value, data=trans_dt[type=='zscore'])

stargazer(mod1, mod2, mod3, mod4, mod5, 
          type='html') %>%
  capture.output(file=file.path(out.dir, 'table_4.html'))

#combine into a data.table for graphing
results_dt <- list(
  
  tidy(mod1, conf.int = TRUE) %>% as.data.table %>% .[, model := 'raw'],
  tidy(mod2, conf.int = TRUE) %>% as.data.table %>% .[, model := 'centile'],
  tidy(mod3, conf.int = TRUE) %>% as.data.table %>% .[, model := 'decile'],
  tidy(mod4, conf.int = TRUE) %>% as.data.table %>% .[, model := 'minmax'],
  tidy(mod5, conf.int = TRUE) %>% as.data.table %>% .[, model := 'zscore']
  
) %>% rbindlist %>% 
  .[, indicator := str_replace_all(term, 'indicator', '') %>% 
      str_replace_all(':value', '')] %>% 
  .[, term := str_replace_all(indicator, '_', ' ')] %>% 
  merge(theme_labels_dt, by='indicator') %>% 
  .[, theme_ind := paste0(theme, '\n', term %>% str_to_title)]

ggplot(data = results_dt[term!='(Intercept)'&model!='raw'], 
       aes(x = estimate, y = forcats::fct_reorder(theme_ind, estimate), xmin = conf.low, xmax = conf.high, 
           color = model, shape = model)) +
  geom_pointrange(position = position_dodge(width = 0.5)) +
  geom_vline(xintercept=0, linetype='dotted') +
  #coord_flip() +
  scale_y_discrete('') +
  scale_x_continuous('', limits=c(-1.2, .1)) +
  scale_color_manual('', values=viridis::turbo(10)[c(1,2,7,8,10)]) +
  scale_shape_discrete(guide='none') +
  ggtitle('Differential Relationship with Life Expectancy based on Indicator Transformations') +
  ggpubr::theme_pubclean(flip = FALSE)

file.path(viz.dir, 'trans_v_le_regression.png') %>% ggsave(height=8, width=12)

#***********************************************************************************************************************

# ---PCA SCRAP----------------------------------------------------------------------------------------------------------

##PCA modelling##
#generate dataset to run PCA 
pc_dt <-
  dt[level==3, .(GEOID, theme, item_short, rank, rank_v1,
                 measure, measure_v1, rank_shift, measure_shift, dropout_factor=dropout)] %>% 
  merge(dt[level==1, .(GEOID, overall=rank, overall_v1=rank_v1, overall_shift=rank_shift)], 
        by='GEOID')

#also pull out the theme map in order to use for plotting
theme_dt <- dt[level==3, .(item_short, theme)] %>% 
  unique(by='item_short') %>% 
  .[theme=='Sensitive Populations', theme_short := 'Pops: '] %>% 
  .[theme=='Socioeconomic Factors', theme_short := 'Socio: '] %>% 
  .[theme=='Environmental Exposures', theme_short := 'Env Exp: '] %>% 
  .[theme=='Environmental Effects', theme_short := 'Env Fx: '] %>% 
  .[, item_themes := paste0(theme_short, item_short)]

#pca model#1 = pca on current raw measure
pca_raw_mod <- pcWrappeR(pc_dt, var='measure', tag='_raw_measure_')

#pca model#2 = pca on current rank
pca_rank_mod <- pcWrappeR(pc_dt, var='rank', tag='_rank_')

#pca model#3 = pca on current rank, but:
#impute missing variables with minimum then run
#TODO try instead imputePCA command in the missMDA? 
pc_imp_dt <-
  dt[level==3, .(GEOID, theme, item_short, rank, rank_v1, ruca_level,
                 measure, measure_v1, rank_shift, measure_shift, dropout_factor=dropout)] %>% 
  merge(dt[level==1, .(GEOID, overall=rank, overall_v1=rank_v1, overall_shift=rank_shift, 
                       impacted, impacted_hierarchy)], by='GEOID') %>% 
  .[, rank := na.aggregate(rank, FUN= min) , by = item_short]

pred_imp_dt <- pc_imp_dt[, .(GEOID, item_short, rank, overall, impacted, impacted_hierarchy, ruca_level)] %>% 
  dcast(...~item_short, value.var=c('rank')) %>% 
  .[is.na(`Wastewater Discharge`), 'Wastewater Discharge' := 0] 


pca_rankimp_mod <- pcWrappeR(pc_imp_dt, pred_dt=pred_imp_dt, var='rank', tag='_rankimp_')

#also run a stratified version by urbanicity
pca_rankimp_strats <- lapply(unique(pc_imp_dt$ruca_level),
                             pcWrappeR,
                             dt=pc_imp_dt, 
                             pred_dt=pred_imp_dt, 
                             var='rank', strat_var='ruca_level',
                             tag='_rankimp_')

extractStratLoadings <- function(i, strats, lvls) {
  strats[[i]]$mod$rotation %>%
    .[,1:2] %>% 
    data.frame %>%  
    setDT(keep.rownames = T) %>% 
    setnames(., 'rn', 'var') %>% 
    .[, ruca_level := lvls[i]] %>% 
    return
}

strat_loadings <- lapply(1:length(pca_rankimp_strats), extractStratLoadings,
                         strats=pca_rankimp_strats,
                         lvls=unique(dt$ruca_level)) %>% 
  rbindlist %>% 
  melt(id.vars=c('var', 'ruca_level'), value.name='pc', variable.name='component') %>% 
  merge(theme_dt, by.x='var', by.y='item_short')

ggplot(strat_loadings, aes(fill=ruca_level, y=pc, x=item_themes)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_viridis(discrete=TRUE, name="") +
  facet_wrap(~component) +
  coord_flip() +
  ylab("PC Loading") + 
  xlab("") +
  theme_minimal()

file.path(viz.dir, 'stratified_pc_loadings.png') %>% ggsave(height=8, width=12)

#also create a facetted lolli of the centroids

ggplot( aes(x=rowname, y=V1)) +
  geom_segment( aes(x=rowname ,xend=rowname, y=V2, yend=V1), color="grey") +
  geom_point(size=5, color="#69b3a2") +
  geom_point(aes(y=V2), size=5, color="#69b3a2", alpha=0.1) +
  coord_flip() +
  theme_ipsum() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text( size=48 )
  ) +
  ylim(0,20) +
  ylab("mark") +
  xlab("")

file.path(viz.dir, 'stratified_clusters.png') %>% ggsave(height=8, width=12)

  #also create a 3d plot of the components
plot_ly(data=pca_rankimp_mod$preds, x=~PC1, y=~PC2, z=~PC3, 
        type="scatter3d", mode="markers", color=~overall, stroke=~impacted*-1)

plot_ly(data=pca_rankimp_mod$preds, x=~PC1, y=~PC2, z=~PC3, 
        type="scatter3d", mode="markers", color=~impacted_hierarchy, colors='Paired')

#regress the PCs against the overall ranking
lm_dt <-
pca_rankimp_mod$preds %>% 
  .[, .(GEOID, PC1, PC2, PC3)] %>% 
  merge( dt[level==3, .(GEOID, rank, impacted, dropout)], by='GEOID')

mod1 <- lm(rank~PC1+PC2+PC3, data=lm_dt)
mod2 <- lm(impacted~PC1+PC2+PC3, data=lm_dt)

stargazer(mod1, mod2,
          type='html') %>%
  capture.output(file=file.path(out.dir, 'table_1.html'))

#examine correlations between pcs and the variables
#reshape wide in order to run cor
corr_dt <- pca_rankimp_mod$preds[, -c('GEOID', 'overall', 'impacted', 'impacted_hierarchy', 'level'), with=F] %>% 
  na.omit

#generate correlations
corr_dt <- corr_dt %>% cor
p_mat <- corr_dt %>% cor_pmat

ggcorrplot(corr_dt, type = "upper",
           colors = c("#6D9EC1", "white", "#E46726")) +
  theme(axis.text.x=element_text(angle=45,hjust=1))  

file.path(viz.dir, paste0('pca_rankimp_corr.png')) %>% ggsave(height=8, width=12)

#create a special quadrant plot based on urb:rur wealth gradient
#find the min maxes of the first two components
xlims <- c(min(pca_rankimp_mod$mod$x[,1]) %>% floor, max(pca_rankimp_mod$mod$x[,1]) %>% ceiling)
ylims <- c(min(pca_rankimp_mod$mod$x[,2]) %>% floor, max(pca_rankimp_mod$mod$x[,2]) %>% ceiling)
crossAt <- 0
colours <- pals::brewer.seqseq2(n = 9)

fviz_pca_ind(pca_rankimp_mod$mod,
             label='none',
             col.var='black',
             habillage = pca_rankimp_mod$biplot_dt$impacted_hierarchy, # Color by the quality of representation
             repel = TRUE     # Avoid text overlapping
) +
  annotate("rect", xmin = xlims[1], xmax = crossAt, ymin = ylims[2], ymax = crossAt, fill = colours[8], alpha=.25) + 
  annotate("rect", xmin = crossAt, xmax = xlims[2], ymin = crossAt, ymax = ylims[2], fill = colours[7], alpha=.25)  +
  annotate("rect", xmin = xlims[1], xmax = crossAt, ymin = ylims[1], ymax = crossAt , fill= colours[9], alpha=.15) + 
  annotate("rect", xmin = crossAt, xmax = xlims[2], ymin = crossAt, ymax = ylims[1], fill = colours[1], alpha=.25) + 
  scale_shape_manual(values=c(16, 16, 16, 16), guide=F) +
  scale_color_brewer('Impact Status', palette = 'Paired') +
  theme_minimal()

file.path(viz.dir, paste0('pc_ind_presentation_quadrants.png')) %>% ggsave(height=8, width=12)

#plot the quadrant classifications on a map along w impacted status
quad_dt <- pca_rankimp_mod$preds %>% 
  copy %>% 
  .[, quadrant := 0] %>% 
  .[PC1<0 & PC2>0, quadrant := 1] %>% 
  .[PC1>0 & PC2>0, quadrant := 2] %>%
  .[PC1<0 & PC2<0, quadrant := 3] %>%
  .[PC1>0 & PC2<0, quadrant := 4] %>%
  .[, quadrant := factor(quadrant,
                         levels=c(1, 2, 3, 4),
                         labels=c('Urban High Deprivation',
                                  'Rural High Deprivation',
                                  'Urban Low Deprivation',
                                  'Rural Low Deprivation'))]
#build a color scheme
quad_colors <- pals::brewer.seqseq2(n = 9)[c(8,7,9,1)]

cartographeR(dt=quad_dt, map_varname = 'quadrant', map_label = 'PC Quadrant',
             map_title = paste0('Mapping the urbanicity and deprivation gradient based on PCA'),
             scale_type='cont_man',
             scale_vals=quad_colors,
             tag='pc_quadrants',
             lvl=3)

cartographeR(dt=quad_dt[impacted_hierarchy%like%'Yes'], map_varname = 'quadrant', map_label = 'PC Quadrant',
             map_title = paste0('Mapping the urbanicity and deprivation gradient based on PCA'),
             scale_type='cont_man',
             scale_vals=quad_colors,
             tag='pc_quadrants_impacted',
             lvl=3)

#plot the quadrant classifications on a map along w impacted status
#do the same for v1
quad_dt <- pca_rankimpv1_mod$preds %>% 
  copy %>% 
  .[, quadrant := 0] %>% 
  .[PC1<0 & PC2>0, quadrant := 1] %>% 
  .[PC1>0 & PC2>0, quadrant := 2] %>%
  .[PC1<0 & PC2<0, quadrant := 3] %>%
  .[PC1>0 & PC2<0, quadrant := 4] %>%
  .[, quadrant := factor(quadrant,
                         levels=c(1, 2, 3, 4),
                         labels=c('Urban High Deprivation',
                                  'Rural High Deprivation',
                                  'Urban Low Deprivation',
                                  'Rural Low Deprivation'))]
#build a color scheme
quad_colors <- pals::brewer.seqseq2(n = 9)[c(9,1,8,7)]

cartographeR(dt=quad_dt, map_varname = 'quadrant', map_label = 'PC Quadrant',
             map_title = paste0('Mapping the urbanicity and deprivation gradient based on PCA'),
             scale_type='cont_man',
             scale_vals=quad_colors,
             tag='pc_quadrants_v1',
             lvl=3)

cartographeR(dt=quad_dt[impacted_hierarchy%like%'Yes'], map_varname = 'quadrant', map_label = 'PC Quadrant',
             map_title = paste0('Mapping the urbanicity and wealth gradient based on PCA'),
             scale_type='cont_man',
             scale_vals=quad_colors,
             tag='pc_quadrants_impacted_v1',
             lvl=3)

##PLS modeling##
#pls model#1: pls on current rank against overall rank#set.seed(98118) #set seed for replication
#build a dataset, with imputation
pls_dt <-
  dt[level==3, .(GEOID, theme, item_short, rank, rank_v1,
                 measure, measure_v1, rank_shift, measure_shift, impacted,
                 le=substr(life_expectancy, 1, 4) %>% as.numeric)] %>% 
  merge(dt[level==1, .(GEOID, overall=rank, overall_v1=rank_v1, overall_shift=rank_shift)], by='GEOID') %>% 
  .[, measure := na.aggregate(measure, FUN= min) , by = item_short] %>%
  .[, rank := na.aggregate(rank, FUN= min) , by = item_short]

pls_pred_dt <- pls_dt[, .(GEOID, item_short, rank, overall)] %>% 
  dcast(...~item_short, value.var=c('rank')) %>% 
  na.omit
#run model
pls_rank_mod <- pcWrappeR(pls_dt, pred_dt=pls_pred_dt, var='rank', depvar='overall', tag='_rank_vs_ovr_')


#***********************************************************************************************************************

# ---SCRAP -------------------------------------------------------------------------------------------------------------
# hoard your scraps here 
# #***********************************************************************************************************************