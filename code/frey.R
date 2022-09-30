# ----HEADER------------------------------------------------------------------------------------------------------------
# Author: JF
# Date: 04/08/2022
# Purpose: Prepare EHD rankings from raw data, explore and analyze sensitivity
# source("/homes/jfrostad/_code/ehd_mapsense/frey.R", echo=T)
#***********************************************************************************************************************

# ----CONFIG------------------------------------------------------------------------------------------------------------
# clear memory
rm(list=ls())

#set opts
reload <- F #set true if you want to reprep all the data
options(scipen=999) #readability
#use cairo to render instead of quartz (quartz causes big slowdowns with geom_sf)
if(!identical(getOption("bitmapType"), "cairo") && isTRUE(capabilities()[["cairo"]])){
  options(bitmapType = "cairo")
}

## Set core_repo location
user            <- Sys.info()['user']
main.dir         <- ifelse(Sys.info()["sysname"] == "Linux",
                           file.path('/homes', user, ''),
                           file.path('C:/Users', user, 'Documents/ehd_mapsense/'))
my_repo <- file.path(main.dir, 'repo')
setwd(my_repo)

#load packages
#TODO only relevant to running in linux on shared cluster
package_lib    <- sprintf('%s_code/_lib/pkg_R', my_repo)
## Load libraries and  MBG project functions.
.libPaths(package_lib)

#TODO cleanup old packages
pacman::p_load(tidyverse, readxl, snakecase, janitor, data.table, naniar, visdat,
               magrittr, scales, ggplot2, ggpubr, ggridges, ggrepel, gridExtra, RColorBrewer, 
               sf, viridis, farver, reldist, ggnewscale, ggallin,
               tigris, tidycensus, ggcorrplot,
               broom.mixed, ggstance, jtools, factoextra,
               stargazer,
               caret, mlbench, randomForest,
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
viz.dir  <- file.path(main.dir, 'viz')
vizdata.dir  <- file.path(my_repo, 'code/baldR/data')
#***********************************************************************************************************************

# ---FUNCTIONS----------------------------------------------------------------------------------------------------------
#source custom functions that are relevant to this module
file.path(code.dir, '_lib', 'prep_fx.R') %>% source
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

#life expectancy data
le_dt <-  file.path(data.dir, 'le_at_birth_2015_2019.csv') %>% fread
setnames(le_dt, names(le_dt), 
         c('county', 'geocode', 'le', 'lower', 'upper'))
le_dt[, c('le', 'lower', 'upper') := lapply(.SD, as.numeric), .SDcols=c('le', 'lower', 'upper')]

#also bring in the census tracts shapefile in order to do some cartography
#can be downloaded from the census website using tigris
tract_sf <- tracts('WA', year=2010, cb=T) %>% 
  st_transform(32148) %>% 
  erase_water(area_threshold = 0.9) %>% #intersect with water overlay and remove
  mutate('GEOID'=substring(GEO_ID, 10)) #remove the excess first 9 chr and rename GEOID

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
ranks_new <- rankeR(dir=data.dir, path=data_extract_EHDv2, debug=T) 

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

#add life expectancy data
index_dt <- merge(index_dt, le_dt, by.x='GEOID', by.y='geocode')
index_dt[, le_state_average := mean(le, na.rm=T)]
index_dt[, life_expectancy := paste0(le, ' years (', lower, '-', upper, ')')]

#merge the dropouts back onto the other tables for graphing
measure_dt <- merge(measure_dt, index_dt[, .(GEOID, dropout, index=rank, 
                                             impacted, impacted_v1, life_expectancy)], by='GEOID')
theme_dt <- merge(theme_dt, index_dt[, .(GEOID, dropout, index=rank, 
                                         impacted, impacted_v1, life_expectancy)], by='GEOID')

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
dt <- list(index_dt[, -c('le', 'lower', 'upper', 'le_state_average', 'county'), with=F],
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
#update colors with allies scale
# cont_colors <- c(
#   rgb(248, 250, 251, maxColorValue = 255), # 1:     Red: 248, Green: 250, Blue: 251
#   rgb(236, 242, 246, maxColorValue = 255), # 2:     Red: 236, Green: 242, Blue: 246
#   rgb(220, 230, 239, maxColorValue = 255), # 3:     Red: 220, Green: 230, Blue: 239
#   rgb(203, 218, 233, maxColorValue = 255), # 4:     Red: 203, Green: 218, Blue: 233
#   rgb(194, 199, 223, maxColorValue = 255), # 5:     Red: 194, Green: 199, Blue: 223
#   rgb(194, 178, 213, maxColorValue = 255), # 6:     Red: 194, Green: 178, Blue: 213
#   rgb(192, 157, 203, maxColorValue = 255), # 7:     Red: 192, Green: 157, Blue: 203
#   rgb(189, 132, 186, maxColorValue = 255), # 8:     Red: 189, Green: 132, Blue: 186
#   rgb(161, 124, 177, maxColorValue = 255), # 9:     Red: 161, Green: 124, Blue: 177
#   rgb(163, 124, 162, maxColorValue = 255) # 10:   Red: 163, Green: 124, Blue: 162
# )

names(cont_colors) <- 1:10


#drop water codes with weird data from maps
drop_geocodes <- c('53057990100' #san juan water area with lots of big changes
)

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
               map_title = 'Overall Rank (V2.0)', scale_type='cont_man', scale_vals = cont_colors)
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
  
  # #loop through changes in the measures
  # lapply(unique(measure_raw$item),
  #        cartographeR,
  #        shapefile=tract_sf,
  #        dt=measure_raw, map_varname = 'measure_shift_raw', 
  #        subset_var='item', scale_type='cont')
  
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
#***********************************************************************************************************************
 
# ---PLOT--------------------------------------------------------------------------------------------------------------
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
ggplot(index_dt, aes(x=index_new, y=le, color=impacted_new %>% as.factor) ) +
  geom_point(position='jitter') +
  #geom_hex(bins = 70) +
  scale_color_brewer('Dropouts', palette='Pastel1') +
  theme_bw() 

#create plot of life expectancy distributions
ggplot(index_dt, aes(rank %>% as.factor, le)) + 
  geom_violin(aes(col = rank %>% as.factor, fill = rank %>% as.factor), alpha = 0.25)+
  labs(x = "Index", y = "Life Expectancy",
       title = "Relationship between EHD Rank and Life Expectancy",
       subtitle = "Average Life Expectancy Between 2015-2019 (dotted line represents state average of 80.3 years)"
  ) +
  geom_vline(xintercept = 8.5, linetype='dashed', color='dark red') +
  geom_hline(yintercept = index_dt[1, le_state_average], linetype='dotted', color='dark blue') +
  scale_x_discrete() +
  scale_color_manual('Rank', values=cont_colors, na.value = "grey75") +
  scale_fill_manual('Rank', values=cont_colors, na.value = "grey75") +
  geom_boxplot(color = "gray20", width = 0.15, coef = 1.5) +
  annotate("text", x = 9.25, y = 99, label = "Highly Impacted Tracts", vjust = -0.5, color='dark red') +
  theme_minimal() +
  theme(legend.position = c(1, .99), legend.justification = c(1, 1),
        plot.margin = unit(c(0, 0, 0, 0), "in"))

#save the plot
file.path(viz.dir, 'le_violin_minimal.png') %>% ggsave(height=8, width=12)

#scatter the measure ranks to compare V1:v2
ggplot(dt[level==3], aes(x=rank_v1, rank, color=dropout) ) +
  geom_point(position='jitter') +
  #geom_hex(bins = 70) +
  facet_wrap(~item_short) + 
  scale_x_continuous("Indicator Ranking, V1.1") +
  scale_y_continuous("Indicator Ranking, V2.0") +
  scale_color_brewer('Highly Impacted: \nDropouts', palette='Pastel1') +
  theme_minimal() 
#save the plot
file.path(viz.dir, 'measure_shift_scatters.png') %>% ggsave(height=8, width=12)

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
  
#***********************************************************************************************************************

# ---ANALYZE------------------------------------------------------------------------------------------------------------
##variable importance investigation##
#simulation analysis#
#first read in and calculate all the ranks using custom function
ranks_per <- lapply(unique(dt$item)[-1] %>% as.character,
                    rankeR,
                    dir=data.dir, 
                    path=data_extract_EHDv2,
                    nranks=10, 
                    clean_names=F) %>% 
    rbindlist %>% 
    .[, .(GEOID, permute_rank=index_rank_integer, permutation)]
  
ranks_per <- merge(ranks_per,
                   dt[level==1],
                   by='GEOID',
                   allow.cartesian=T)  

#also merge on the theme names for coloring

ranks_per <- merge(ranks_per[, `:=`(item = NULL, theme = NULL)],
             dt[level==3, .(item, theme)] %>% unique(by='item'),
             by.x='permutation',
             by.y = 'item')  

#identify dropout units based on impacted threshold of >8
threshold_val <- 9
ranks_per[, permute_impacted := 0]
ranks_per[permute_rank>=threshold_val, permute_impacted := 1]
ranks_per[, permute_shift := permute_rank-rank]

rank_uncertainty <- ranks_per[, range:=max(permute_rank)-min(permute_rank), by=GEOID] %>% 
  .[, rank_acc := sum(permute_rank==rank)/.N, by=GEOID] %>% 
  .[, impacted_acc := sum(impacted==permute_impacted)/.N, by=GEOID] %>% 
  unique(by='GEOID')

cartographeR(dt=rank_uncertainty, map_varname = 'range', map_label = 'Range in Ranking',
             map_title = 'Rank Range', scale_type='cont')
cartographeR(dt=rank_uncertainty, map_varname = 'rank_acc', map_label = 'Accuracy',
             map_title = 'Accuracy of Ranking', scale_type='cont')
cartographeR(dt=rank_uncertainty, map_varname = 'impacted_acc', map_label = 'Accuracy',
             map_title = 'Accuracy of Impacted Classification', scale_type='cont')

ggplot(rank_uncertainty, aes(x=rank, range, color=rank_shift %>% as.factor) ) +
  geom_point(position='jitter') +
  scale_color_viridis_d('Change in Rank') +
  theme_bw() 

file.path(viz.dir, 'rank_range_scatters.png') %>% ggsave(height=8, width=12)

ggplot(rank_uncertainty, aes(x=rank, rank_acc, color=rank_shift %>% as.factor) ) +
  geom_point(position='jitter') +
  scale_color_viridis_d('Change in Rank') +
  theme_bw() 

file.path(viz.dir, 'rank_acc_scatters.png') %>% ggsave(height=8, width=12)

var_imp <- ranks_per.[, importance_mean := mean(abs(permute_rank-rank)), by=permutation] %>% 
  .[, importance_max := max(abs(permute_rank-rank)), by=permutation] %>% 
  .[, importance_acc :=sum(impacted==permute_impacted)/.N, by=permutation] %>% 
  unique(by='permutation')

ggplot(var_imp, aes(x=forcats::fct_reorder(permutation, importance_max), 
                    y=importance_max,
                    color=theme)) +
  geom_point() + 
  geom_segment( aes(x=permutation, xend=permutation, y=0, yend=importance_max)) +
  theme_bw() +
  scale_color_brewer('Theme', palette='Paired') +
  coord_flip()
file.path(viz.dir, 'rank_importance_max_lolli.png') %>% ggsave(height=8, width=12)

ggplot(var_imp, aes(x=forcats::fct_reorder(permutation, importance_mean), 
                    y=importance_mean,
                    color=theme)) +
  geom_point() + 
  geom_segment( aes(x=permutation, xend=permutation, y=0, yend=importance_mean)) +
  theme_bw() +
  scale_color_brewer('Theme', palette='Paired') +
  coord_flip()
file.path(viz.dir, 'rank_importance_mean_lolli.png') %>% ggsave(height=8, width=12)

ggplot(var_imp, aes(x=forcats::fct_reorder(permutation, importance_acc), 
                    y=importance_acc,
                    color=theme)) +
  geom_point() + 
  geom_segment( aes(x=permutation, xend=permutation, y=0, yend=importance_acc)) +
  theme_bw() +
  scale_color_brewer('Theme', palette='Paired') +
  coord_flip()
file.path(viz.dir, 'rank_importance_acc_lolli.png') %>% ggsave(height=8, width=12)



ggplot(ranks_per, aes(x=permute_shift, y=forcats::fct_reorder(permutation, permute_shift, .fun=mean),
                                                     )) +
  geom_density_ridges(stat='binline') +
  scale_fill_viridis_c(name = "Index Ranking Change", direction = -1, option='plasma') +
  scale_y_discrete('Permuted Item') +
  theme_bw()
file.path(viz.dir, 'rank_importance_ridges.png') %>% ggsave(height=8, width=12)


#regression analysis#
#use logistic regression to find most influential vars
reg_dt <- rbind(
  dt[level==3, .(GEOID, theme, item_short,
                 measure, measure_v1, rank_shift, measure_shift, shift_type='shift',
                 impacted, impacted_v1, dropout_factor=dropout)],
  dt[level==3, .(GEOID, theme, item_short,
                 measure, measure_v1, measure_shift=measure_ratio, shift_type='ratio',
                 impacted, impacted_v1, dropout_factor=dropout)],
  fill=T) %>% 
  merge(dt[level==1, .(GEOID, overall=rank, overall_v1=rank_v1, overall_shift=rank_shift)], by='GEOID') %>% 
  .[, dropout := dropout_factor=='Dropout'] %>% 
  .[, dropin := dropout_factor=='Addition'] %>% 
  na.omit

#create a transformed version of measure shift to deal with skewness
reg_dt[, measure_shift_log := sign(measure_shift)*log(abs(measure_shift)+0.0001)]
reg_dt[, measure_shift_scale := scale(measure_shift)]
reg_dt[, item_fac := paste0(theme, ': ', item_short) %>% factor]

#
ggplot(reg_dt, aes(x=overall_shift, rank_shift, color=overall %>% as.factor) ) +
  geom_point(position='jitter') +
  #geom_hex(bins = 70) +
  facet_wrap(~item_short) + 
  facet_wrap(~item_short) + 
  scale_color_viridis_d('Overall Rank') +
  theme_bw() 

#save the plot
file.path(viz.dir, 'rank_shift_scatters.png') %>% ggsave(height=8, width=12)

#
ggplot(reg_dt, aes(x=overall_shift, measure_shift_log, color=overall %>% as.factor) ) +
  geom_point(position='jitter') +
  #geom_hex(bins = 70) +
  facet_wrap(~item_short) + 
  scale_color_viridis_d('Overall Rank') +
  theme_bw() 

#save the plot
file.path(viz.dir, 'measure_shift_scatters.png') %>% ggsave(height=8, width=12)

#generate dataset to run PCA 
pca_dt <-
  dt[level==3, .(GEOID, theme, item_short, rank, rank_v1,
                 measure, measure_v1, rank_shift, measure_shift, dropout_factor=dropout)] %>% 
  merge(dt[level==1, .(GEOID, overall=rank, overall_v1=rank_v1, overall_shift=rank_shift)], by='GEOID')

#pca model#1 = pca on current raw measure
pca_mod <- pca_dt[, .(GEOID, item_short, measure)] %>% 
  dcast(GEOID~item_short, value.var=c('measure'))  %>% 
  na.omit %>% 
  .[, -c('GEOID'), with=F] %>% 
  prcomp(scale=T)

#screeplot
fviz_eig(pca_mod)
file.path(viz.dir, 'pc_raw_eig.png') %>% ggsave(height=8, width=12)
fviz_contrib(pca_mod, choice='var', axes=1)
file.path(viz.dir, 'pc_raw_contribution_v1.png') %>% ggsave(height=8, width=12)
fviz_contrib(pca_mod, choice='var', axes=2)
file.path(viz.dir, 'pc_raw_contribution_v2.png') %>% ggsave(height=8, width=12)

fviz_pca_var(pca_mod,
             col.var = "contrib", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
file.path(viz.dir, 'pc_raw_var.png') %>% ggsave(height=8, width=12)

#pca model#2 = pca on current rank
pca_mod <- pca_dt[, .(GEOID, item_short, rank)] %>% 
  dcast(GEOID~item_short, value.var=c('rank'))  %>% 
  na.omit %>% 
  .[, -c('GEOID'), with=F] %>% 
  prcomp(scale=T)

#screeplot
fviz_eig(pca_mod)
file.path(viz.dir, 'pc_rank_eig.png') %>% ggsave(height=8, width=12)
fviz_contrib(pca_mod, choice='var', axes=1)
file.path(viz.dir, 'pc_rank_contribution_v1.png') %>% ggsave(height=8, width=12)
fviz_contrib(pca_mod, choice='var', axes=2)
file.path(viz.dir, 'pc_rank_contribution_v2.png') %>% ggsave(height=8, width=12)
fviz_contrib(pca_mod, choice='var', axes=3)
file.path(viz.dir, 'pc_rank_contribution_v3.png') %>% ggsave(height=8, width=12)

fviz_pca_var(pca_mod,
             col.var = "contrib", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
file.path(viz.dir, 'pc_rank_var.png') %>% ggsave(height=8, width=12)

#make predictions of pca #1/2s so that we can map them 
pred_dt <- pca_dt[, .(GEOID, item_short, rank)] %>% 
  dcast(GEOID~item_short, value.var=c('rank')) %>% 
  cbind(., predict(pca_mod, newdata=.)[,1:3]) %>% 
  .[,level:=3] #remind data that its most granular

#map the first 3 elements of PCA
cartographeR(dt=pred_dt, map_varname = 'PC1', map_label = 'PCA #1',
             map_title = 'Mapping the first component of PCA using the ranked measures',
             scale_type='cont_grad', 
             lvl=3)
cartographeR(dt=pred_dt, map_varname = 'PC2', map_label = 'PCA #2',
             map_title = 'Mapping the second component of PCA using the ranked measures',
             scale_type='cont_grad', 
             lvl=3)
cartographeR(dt=pred_dt, map_varname = 'PC3', map_label = 'PCA #3',
             map_title = 'Mapping the third component of PCA using the ranked measures',
             scale_type='cont_grad', 
             lvl=3)

#impute missing variables with minimum then run
#TODO try instead imputePCA command in the missMDA? 
pca_dt <-
  dt[level==3, .(GEOID, theme, item_short, rank, rank_v1,
                 measure, measure_v1, rank_shift, measure_shift, dropout_factor=dropout)] %>% 
  merge(dt[level==1, .(GEOID, overall=rank, overall_v1=rank_v1, overall_shift=rank_shift)], by='GEOID') %>% 
  .[, rank := na.aggregate(rank, FUN= min) , by = item_short]

#pca model#3 = pca on current rank with minimputation
pca_mod <- pca_dt[, .(GEOID, item_short, rank)] %>% 
  dcast(GEOID~item_short, value.var=c('rank'))  %>% 
  na.omit %>% 
  .[, -c('GEOID'), with=F] %>% 
  prcomp(scale=T)


#screeplot
fviz_eig(pca_mod)
file.path(viz.dir, 'pc_rank_min_eig.png') %>% ggsave(height=8, width=12)
fviz_contrib(pca_mod, choice='var', axes=1)
file.path(viz.dir, 'pc_rank_min_contribution_v1.png') %>% ggsave(height=8, width=12)
fviz_contrib(pca_mod, choice='var', axes=2)
file.path(viz.dir, 'pc_rank_min_contribution_v2.png') %>% ggsave(height=8, width=12)
fviz_contrib(pca_mod, choice='var', axes=3)
file.path(viz.dir, 'pc_rank_min_contribution_v3.png') %>% ggsave(height=8, width=12)

fviz_pca_var(pca_mod,
             col.var = "contrib", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
file.path(viz.dir, 'pc_rank_min_var.png') %>% ggsave(height=8, width=12)

pred_dt <- pca_dt[, .(GEOID, item_short, rank)] %>% 
  dcast(GEOID~item_short, value.var=c('rank')) %>% 
  .[is.na(`Wastewater Discharge`), 'Wastewater Discharge' := 0] %>%  #impute wastewater
  cbind(., predict(pca_mod, newdata=.)[,1:3]) %>% 
  .[,level:=3] #remind data that its most granular

#map the first 3 elements of PCA
cartographeR(dt=pred_dt, map_varname = 'PC1', map_label = 'PCA #1',
             map_title = 'Mapping the first component of PCA using the ranked measures',
             scale_type='cont_grad', 
             lvl=3)
cartographeR(dt=pred_dt, map_varname = 'PC2', map_label = 'PCA #2',
             map_title = 'Mapping the second component of PCA using the ranked measures',
             scale_type='cont_grad', 
             lvl=3)
cartographeR(dt=pred_dt, map_varname = 'PC3', map_label = 'PCA #3',
             map_title = 'Mapping the third component of PCA using the ranked measures',
             scale_type='cont_grad', 
             lvl=3)

#regress the PCs against the overall ranking
mod_dt <- pred_dt[, .(GEOID, PC1, PC2, PC3)] %>% 
  merge( dt[level==3, .(GEOID, rank, impacted, dropout)], by='GEOID') 
  
mod <- lm(rank~PC1+PC2+PC3, data=mod_dt)
mod <- lm(impacted~PC1+PC2+PC3, data=mod_dt)
mod <- lm(dropout~PC1+PC2+PC3, data=mod_dt)
#use machine learning models to # #try first with a rf model
# # prepare training scheme
# control <- trainControl(method="repeatedcv", number=5, repeats=3)
# # train the model
# model <- mod_dt[, .(GEOID, item_short, rank, overall)] %>% 
#   dcast(GEOID+overall~item_short, value.var=c('rank')) %>% 
#   .[, -c('GEOID'), with=F] %>% 
#   na.omit %>% 
#   train(overall~., data=., method="rf", preProcess="scale", trControl=control)
# # estimate variable importance
# importance <- varImp(model, scale=FALSE)
# 
# # train the model
# model <- mod_dt[, .(GEOID, item_short, rank, overall)] %>% 
#   dcast(GEOID+overall~item_short, value.var=c('rank')) %>% 
#   .[, -c('GEOID'), with=F] %>% 
#   na.omit %>% 
#   train(overall~., data=., method="pls", preProcess="scale", trControl=control)
# # estimate variable importance
# importance <- varImp(model, scale=FALSE)do some feature selection
#build a dataset
mod_dt <-
  dt[level==3, .(GEOID, theme, item_short, rank, rank_v1,
                 measure, measure_v1, rank_shift, measure_shift, dropout_factor=dropout,
                 le=substr(life_expectancy, 1, 4) %>% as.numeric)] %>% 
  merge(dt[level==1, .(GEOID, overall=rank, overall_v1=rank_v1, overall_shift=rank_shift)], by='GEOID') %>% 
  .[, measure := na.aggregate(measure, FUN= min) , by = item_short] %>%
  .[, rank := na.aggregate(rank, FUN= min) , by = item_short]

set.seed(98118)

#deepen our pls model by doing some tuning of the grid
mod_dt <- mod_dt[, .(GEOID, item_short, rank, overall)] %>% 
  dcast(GEOID+overall~item_short, value.var=c('rank')) %>% 
  .[, -c('GEOID'), with=F] %>% 
  na.omit
train_scheme <- createDataPartition(mod_dt$overall, p = .80, list = FALSE)
train_dt <- mod_dt[train_scheme,]
test_dt  <- mod_dt[-train_scheme,]
ctrl <- trainControl(
  method = "cv",
  number = 10,
)

tuneGrid <- expand.grid(
  ncomp   = seq(1, 10, by = 1)
)

#show RMSE
model <- train(
  overall ~ .,
  data = train_dt,
  method = 'pls',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid
)

#also run a simpler pls with 3 comps
y_mat <- as.matrix(mod_dt[,1])
x_mat <- as.matrix(mod_dt[,2:ncol(mod_dt)])
mod <- pls::mvr(y_mat  ~ x_mat , ncomp=10, method = "oscorespls" , scale = T)

pc_dt <- mod$Yscores %>% as.matrix %>% .[,1:3] %>% as.data.table %>% setnames(c('PLS_C1', 'PLS_C2', 'PLS_C3'))

pc_dt <- 
  dt[level==3, .(GEOID, theme, item_short, rank, rank_v1,
                 measure, measure_v1, rank_shift, measure_shift, dropout_factor=dropout,
                 le=substr(life_expectancy, 1, 4) %>% as.numeric)] %>% 
  merge(dt[level==1, .(GEOID, overall=rank, overall_v1=rank_v1, overall_shift=rank_shift)], by='GEOID') %>% 
  .[, measure := na.aggregate(measure, FUN= min) , by = item_short] %>%
  .[, rank := na.aggregate(rank, FUN= min) , by = item_short] %>% 
  dcast(GEOID+overall~item_short, value.var=c('rank')) %>% 
  na.omit %>% 
  .[,1] %>% 
  cbind(., pc_dt) %>% 
  .[, level:=3]

#map the first 3 elements of PCA
cartographeR(dt=pc_dt, map_varname = 'PLS_C1', map_label = 'PLS #1',
             map_title = 'Mapping the first component of PLS using the ranked measures',
             scale_type='cont_grad', 
             lvl=3)
cartographeR(dt=pc_dt, map_varname = 'PLS_C2', map_label = 'PLS #2',
             map_title = 'Mapping the second component of PLS using the ranked measures',
             scale_type='cont_grad', 
             lvl=3)
cartographeR(dt=pc_dt, map_varname = 'PLS_C3', map_label = 'PLS #3',
             map_title = 'Mapping the third component of PLS using the ranked measures',
             scale_type='cont_grad', 
             lvl=3)

#deepen our pls model by doing some tuning of the grid
set.seed(98118)
mod_dt <- mod_dt[, .(GEOID, item_short, rank, overall)] %>% 
  dcast(GEOID+overall~item_short, value.var=c('rank')) %>% 
  .[, -c('GEOID'), with=F]

S_plsr <- scores(mod)[, comps= 1:2, drop = FALSE]
cl_plsr <- cor(model.matrix(mod), S_plsr)
df_cor <- as.data.frame(cl_plsr)

df_depend_cor <- as.data.frame(cor(mod_dt[,1], S_plsr))

plot_loading_correlation  <-  rbind(df_cor ,
                                    df_depend_cor)
plot_loading_correlation1 <- setNames(plot_loading_correlation, c("comp1", "comp2"))


#Function to draw circle
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

dat_plsr <- circleFun(c(0,0),2,npoints = 100)

ggplot(data=plot_loading_correlation1 , aes(comp1, comp2))+
  ylab("")+xlab("")+ggtitle(" ")+
  theme_bw() +
  geom_hline(aes(yintercept = 0), size=.2, linetype = 3)+ 
  geom_vline(aes(xintercept = 0), size=.2, linetype = 3)+
  geom_text_repel(aes(label=rownames(plot_loading_correlation1), 
                      colour=ifelse(rownames(plot_loading_correlation1)!='dependent', 'red','darkblue')))+
  scale_color_manual(values=c("red","darkblue"))+
  scale_x_continuous(breaks = c(-1,-0.5,0,0.5,1))+
  scale_y_continuous(breaks = c(-1,-0.5,0,0.5,1))+
  coord_fixed(ylim=c(-1, 1),xlim=c(-1, 1))+xlab("axis 1")+ 
  ylab("axis 2")+ theme(axis.line.x = element_line(color="darkgrey"),
                        axis.line.y = element_line(color="darkgrey"))+
  geom_path(data=dat_plsr ,
            aes(x,y), colour = "darkgrey")+
  theme(legend.title=element_blank())+
  theme(axis.ticks = element_line(colour = "black"))+
  theme(axis.title = element_text(colour = "black"))+
  theme(axis.text = element_text(color="black"))+
  theme(legend.position='none')+
  theme(panel.grid.minor = element_blank()) +
  geom_segment(data=plot_loading_correlation1, aes(x=0, y=0, xend=comp1, yend=comp2), 
               arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, 
               colour=ifelse(rownames(plot_loading_correlation1)=='dependent', 'red','darkblue'))

VIP <- function(object) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}


df_vip  <- as.data.frame(VIP(model$finalModel))

df_coef <- as.data.frame(coef(model$finalModel, ncomp =1:3)) %>% 
  setnames(., c('PC1', 'PC2', 'PC3')) %>% 
  dplyr::mutate(variables = rownames(.)) %>% 
  as.data.table %>% 
  .[, variables := stringr::str_replace_all(variables, stringr::fixed("\\"), '')]  %>% 
  .[, variables := stringr::str_replace_all(variables, stringr::fixed("``"), '')]

#also merge on the theme names for coloring

df_coef <- merge(df_coef,
                   dt[level==3, .(item_short, theme)] %>% unique(by='item_short'),
                   by.x='variables',
                   by.y = 'item_short')  

p2_coef <-
  ggplot(df_coef, aes(x = forcats::fct_reorder(variables, PC1), y = PC1, group = 1, fill=theme))+  
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle=65,
                                   hjust=1,
                                   size = 8),
        axis.title.y = element_text(size = 2))+
  theme_bw()+
  scale_fill_brewer('Theme', palette='Paired') +
  coord_flip()
file.path(viz.dir, 'pls_loadings_pc1.png') %>% ggsave(height=8, width=12)

p2_coef <-
  ggplot(df_coef, aes(x = forcats::fct_reorder(variables, PC2), y = PC2, group = 1, fill=theme))+  
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle=65,
                                   hjust=1,
                                   size = 8),
        axis.title.y = element_text(size = 2))+
  theme_bw()+
  scale_fill_brewer('Theme', palette='Paired') +
  coord_flip()
file.path(viz.dir, 'pls_loadings_pc2.png') %>% ggsave(height=8, width=12)

p2_coef <-
  ggplot(df_coef, aes(x = forcats::fct_reorder(variables, PC3), y = PC3, group = 1, fill=theme))+  
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle=65,
                                   hjust=1,
                                   size = 8),
        axis.title.y = element_text(size = 2))+
  theme_bw()+
  scale_fill_brewer('Theme', palette='Paired') +
  coord_flip()
file.path(viz.dir, 'pls_loadings_pc3.png') %>% ggsave(height=8, width=12)


ggplot(var_imp, aes(x=forcats::fct_reorder(permutation, importance_acc), 
                    y=importance_acc,
                    color=theme)) +
  geom_point() + 
  geom_segment( aes(x=permutation, xend=permutation, y=0, yend=importance_acc)) +
  theme_bw() +
  scale_color_brewer('Theme', palette='Paired') +
  coord_flip()
file.path(viz.dir, 'rank_importance_acc_lolli.png') %>% ggsave(height=8, width=12)


#lets try it again using a classification model on dropouts
mod_dt <-
  dt[level==3, .(GEOID, theme, item_short, rank, rank_v1,
                 measure, measure_v1, rank_shift, measure_shift, impacted=as.factor(impacted),
                 le=substr(life_expectancy, 1, 4) %>% as.numeric)] %>% 
  merge(dt[level==1, .(GEOID, overall=rank, overall_v1=rank_v1, overall_shift=rank_shift)], by='GEOID') %>% 
  .[, measure := na.aggregate(measure, FUN= min) , by = item_short] %>%
  .[, rank := na.aggregate(rank, FUN= min) , by = item_short] %>% 
  .[, .(GEOID, item_short, rank, impacted)] %>% 
  dcast(GEOID+impacted~item_short, value.var=c('measure')) %>% 
  .[, -c('GEOID'), with=F] %>% 
  na.omit

train_dt <- mod_dt[train_scheme,]
test_dt  <- mod_dt[-train_scheme,]
ctrl <- trainControl(
  method = "cv",
  number = 10,
)

tuneGrid <- expand.grid(
  ncomp   = seq(1, 10, by = 1)
)

#show RMSE
model <- train(
  impacted ~ .,
  metric='Accuracy',
  data = train_dt,
  method = 'pls',
  preProcess = c("center", "scale"),
  trControl = ctrl,
  tuneGrid = tuneGrid
)

#run boruta using life expectancy as target
boruta_mod <- mod_dt[, .(GEOID, item_short, measure, le)] %>% 
  dcast(GEOID+le~item_short, value.var=c('measure')) %>% 
  .[, -c('GEOID'), with=F] %>% 
  na.omit %>% 
  .[, rand := rnorm(.N, mean=50, sd=10)] %>% 
  Boruta(le ~ ., data = ., doTrace = 2, maxRuns = 500)
print(boruta_mod)

png(filename=file.path(viz.dir, 'boruta_imp_le.png'),
    height=1200, width=2400, pointsize = 24)
plot(boruta_mod, las = 2, cex.axis = 0.7)
dev.off()

boruta_mod <- mod_dt[, .(GEOID, item_short, measure, overall)] %>% 
  dcast(GEOID+overall~item_short, value.var=c('measure')) %>% 
  .[, -c('GEOID'), with=F] %>% 
  na.omit %>% 
  .[, rand := rnorm(.N, mean=50, sd=10)] %>% 
  Boruta(overall ~ ., data = ., doTrace = 2, maxRuns = 500)
print(boruta_mod)

png(filename=file.path(viz.dir, 'boruta_imp.png'),
                       height=1200, width=2400, pointsize = 24)
plot(boruta_mod, las = 2, cex.axis = 0.7)
dev.off()

#run boruta using life expectancy as target
boruta_mod <- mod_dt[, .(GEOID, item_short, rank, le)] %>% 
  dcast(GEOID+le~item_short, value.var=c('rank')) %>% 
  .[, -c('GEOID'), with=F] %>% 
  na.omit %>% 
  .[, rand := rnorm(.N, mean=50, sd=10)] %>% 
  Boruta(le ~ ., data = ., doTrace = 2, maxRuns = 500)
print(boruta_mod)

png(filename=file.path(viz.dir, 'boruta_imp_rank_le.png'),
    height=1200, width=2400, pointsize = 24)
plot(boruta_mod, las = 2, cex.axis = 0.7)
dev.off()

boruta_mod <- mod_dt[, .(GEOID, item_short, rank, overall)] %>% 
  dcast(GEOID+overall~item_short, value.var=c('rank')) %>% 
  .[, -c('GEOID'), with=F] %>% 
  na.omit %>% 
  .[, rand := rnorm(.N, mean=50, sd=10)] %>% 
  Boruta(overall ~ ., data = ., doTrace = 2, maxRuns = 500)
print(boruta_mod)

png(filename=file.path(viz.dir, 'boruta_imp_rank.png'),
    height=1200, width=2400, pointsize = 24)
plot(boruta_mod, las = 2, cex.axis = 0.7)
dev.off()




# fviz_pca_biplot(pca_mod,
#                  #col.var = "contrib", # Color by the quality of representation
#                 col.ind = 'rank',
#                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#                  label='var',
#                  repel = TRUE     # Avoid text overlapping
# )

#***********************************************************************************************************************

# ---SCRAP -------------------------------------------------------------------------------------------------------------
#model cluster #1
#TODO write desc
#run model
mod1 <- glm(dropout_factor ~ item_fac:measure_shift_log, data=reg_dt[shift_type=='shift'], family='binomial')
mod2 <- glm(dropout ~ item_fac:measure_shift_log, data=reg_dt[shift_type=='shift'], family='binomial')
mod3 <- glm(dropin ~ item_fac:measure_shift_log, data=reg_dt[shift_type=='shift'], family='binomial')
# mod4 <- glm(dropout_factor ~ item:measure_shift, data=reg_dt[shift_type=='ratio'], family='binomial')
# mod5 <- glm(dropout ~ item:measure_shift, data=reg_dt[shift_type=='ratio'], family='binomial')
# mod6 <- glm(dropin ~ item:measure_shift, data=reg_dt[shift_type=='ratio'], family='binomial')
#table
stargazer(mod1, mod2, mod3, #mod4, mod5, mod6,
          type='html') %>% 
  capture.output(file=file.path(out.dir, 'table_1.html'))
# #plot coefficients
# plot_summs(mod1, mod2, mod3, mod4, mod5, mod6,
#            #plot.distributions = TRUE, inner_ci_level = .9,
#            model.names = c("drops vs shift", 
#                            "dropouts vs shift", 
#                            "dropins vs shift",
#                            'drops vs ratio',
#                            'dropouts vs ratio',
#                            'dropins vs ratio')
# )
# file.path(viz.dir, 'regression_1_coefplot.png') %>% ggsave(height=8, width=12)

plot_summs(mod1, mod2, mod3,
           model.names = c("drops vs shift",
                           "dropouts vs shift",
                           "dropins vs shift")
)
file.path(viz.dir, 'regression_1_coefplot.png') %>% ggsave(height=8, width=12)

mod1 <- glm(dropout_factor ~ item_fac:measure_shift_scale, data=reg_dt[shift_type=='shift'], family='binomial')
mod2 <- glm(dropout ~ item_fac:measure_shift_scale, data=reg_dt[shift_type=='shift'], family='binomial')
mod3 <- glm(dropin ~ item_fac:measure_shift_scale, data=reg_dt[shift_type=='shift'], family='binomial')
# mod4 <- glm(dropout_factor ~ item:measure_shift, data=reg_dt[shift_type=='ratio'], family='binomial')
# mod5 <- glm(dropout ~ item:measure_shift, data=reg_dt[shift_type=='ratio'], family='binomial')
# mod6 <- glm(dropin ~ item:measure_shift, data=reg_dt[shift_type=='ratio'], family='binomial')
#table
stargazer(mod1, mod2, mod3, #mod4, mod5, mod6,
          type='html') %>% 
  capture.output(file=file.path(out.dir, 'table_1scale.html'))
# #plot coefficients
# plot_summs(mod1, mod2, mod3, mod4, mod5, mod6,
#            #plot.distributions = TRUE, inner_ci_level = .9,
#            model.names = c("drops vs shift", 
#                            "dropouts vs shift", 
#                            "dropins vs shift",
#                            'drops vs ratio',
#                            'dropouts vs ratio',
#                            'dropins vs ratio')
# )
# file.path(viz.dir, 'regression_1_coefplot.png') %>% ggsave(height=8, width=12)

plot_summs(mod1, mod2, mod3,
           model.names = c("drops vs shift scaled",
                           "dropouts vs shift scaled",
                           "dropins vs shift scaled")
)
file.path(viz.dir, 'regression_1scale_coefplot.png') %>% ggsave(height=8, width=12)

# plot_summs(mod4, mod5, mod6, 
#            model.names = c("drops vs ratio", 
#                            "dropouts vs ratio", 
#                            "dropins vs ratio")
# )
# file.path(viz.dir, 'regression_1_coefplotb.png') %>% ggsave(height=8, width=12)

#model cluster #2
#modelling the shift in rank vs the raw shift
mod1 <- glm(dropout_factor ~ item_fac:rank_shift, data=reg_dt[shift_type=='shift'], family='binomial')
mod2 <- glm(dropout ~ item_fac:rank_shift, data=reg_dt[shift_type=='shift'], family='binomial')
mod3 <- glm(dropin ~ item_fac:rank_shift, data=reg_dt[shift_type=='shift'], family='binomial')
#table 2
stargazer(mod1, mod2, mod3,
          type='html') %>% 
  capture.output(file=file.path(out.dir, 'table_2.html'))

plot_summs(mod1, mod2, mod3,
           model.names = c("drops vs shift",
                           "dropouts vs shift",
                           "dropins vs shift")
)
file.path(viz.dir, 'regression_2_coefplot.png') %>% ggsave(height=8, width=12)

ggplot(reg_dt, aes(x=measure_shift, y=item_short, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1) +
  theme_bw()

file.path(viz.dir, 'shift_ridges.png') %>% ggsave(height=8, width=12)

ggplot(reg_dt, aes(x=measure_shift_log, y=item_short, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1) +
  theme_bw()

file.path(viz.dir, 'shift_ridges_log.png') %>% ggsave(height=8, width=12)

ggplot(reg_dt, aes(x=measure_shift_scale, y=item_short, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1) +
  theme_bw()

file.path(viz.dir, 'shift_ridges_scale.png') %>% ggsave(height=8, width=12)

densityPlot <- function(x, dt=reg_dt, var='measure_shift',
                        fill_var='theme') {
  
  message('plotting ', x)
  
  plot <- 
    ggplot(dt[item_short==x], aes_string(var, fill=fill_var)) +
    geom_density(adjust = 1/5, alpha=.5) +
    ggtitle(x) +
    scale_fill_viridis_d() +
    theme_bw()
  
}

pdf(file.path(viz.dir, 'shift_densities.pdf'), height=8, width=12)
lapply(unique(dt$item_short), densityPlot)
dev.off()

pdf(file.path(viz.dir, 'shift_densities_logged.pdf'), height=8, width=12)
lapply(unique(dt$item_short), densityPlot, var='measure_shift_log')
dev.off()

pdf(file.path(viz.dir, 'shift_densities_scale.pdf'), height=8, width=12)
lapply(unique(dt$item_short), densityPlot, var='measure_shift_scale')
dev.off()

pdf(file.path(viz.dir, 'measure_densities.pdf'), height=8, width=12)
lapply(unique(dt$item_short)[-1], densityPlot, var='value', 
       fill_var='variable',
       dt=reg_dt[, .(GEOID, theme, item_short, measure, measure_v1)] %>% 
         melt(id.vars=c('GEOID', 'theme', 'item_short'))
)
dev.off()


mod2 <- glm(overall_shift ~ item:measure_shift, data=reg_dt, family='gaussian')
stargazer(mod1, mod2, 
          type='html',
          title="Regression Results", dep.var.labels=c("Dropouts","Index Movement"), 
          keep.stat="n", ci=TRUE, ci.level=0.90, single.row=TRUE) %>% 
  capture.output(file=file.path(out.dir, 'table_2.html'))

##explore problems with wastewater data
# cartographeR(dt=ranks_old$measure_raw, map_varname = 'rank_order', scale_type='cont',
#              subset_var='item', subset_val='Wastewater Discharge')
# 
# cartographeR(dt=ranks_new$measure_raw, map_varname = 'rank_order', scale_type='cont',
#              subset_var='item', subset_val='Wastewater Discharge')
# 
# 
# #more complicated version of the violin plot (with labeled outliers)
# ggplot(index_dt, aes(index_new %>% as.factor, le)) + 
#   geom_violin(aes(col = index_new %>% as.factor, fill = index_new %>% as.factor), alpha = 0.25)+
#   labs(x = "Index", y = "Life Expectancy (Years; 2015-2019)",
#        # title = "Temperature rise in India",
#        # subtitle = "Monthly temperature fluctuations between 1961-2019",
#        # caption = "Data source: Kaggle.com"
#   ) +
#   geom_vline(xintercept = 2.5, linetype='dashed', color='dark red') +
#   scale_x_discrete() +
#   scale_color_brewer('Index', palette='PuOr', direction=-1) +
#   scale_fill_brewer('Index', palette='PuOr', direction=-1) +
#   geom_boxplot(color = "gray20", width = 0.15, coef = 1.5) +
#   geom_text_repel(aes(label = outlier_lab), na.rm = TRUE, show.legend = F) +
#   annotate("text", x = 2.1, y = 95, label = "Highly Impacted", vjust = -0.5, color='dark red') +
#   theme_minimal() +
#   theme(legend.position = c(1, 1), legend.justification = c(0, 0),
#         #legend.text=element_text(size=10),
#         plot.title = element_text(hjust=0.5), plot.margin=unit(c(0, 0, 0, 0), "in"))
# 
# #save the plot
# file.path(viz.dir, 'le_violin.png') %>% ggsave(height=8, width=12)
# 
# 
# ggplot(index_dt, aes(x = index_new %>% as.factor, y = le)) + 
#   ggdist::stat_halfeye(
#     adjust = .5, 
#     width = .6, 
#     .width = 0, 
#     justification = -.3, 
#     point_colour = NA) + 
#   geom_boxplot(
#     width = .25, 
#     outlier.shape = NA
#   ) +
#   geom_point(
#     size = 1.3,
#     alpha = .3,
#     position = position_jitter(
#       seed = 1, width = .1
#     )
#   ) + 
#   coord_cartesian(clip = "off") +
#   scale_x_discrete('Index') +
#   scale_y_continuous('Life Expectancy (Years; 2015-2019)') +
#   theme_minimal()
# 
# #save the plot
# file.path(viz.dir, 'le_distribution.png') %>% ggsave(height=8, width=12)
# 
# # Bin size control + color palette
# ggplot(measure_ranks, aes(x=item, fill=avg_shift) ) +
#   geom_bar(aes(y=measure_shift), stat = "summary", fun = "mean") +
#   #geom_hex(bins = 70) +
#   scale_fill_continuous(type = "viridis") +
#   scale_y_continuous(trans = pseudolog10_trans) +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=20,hjust=1))
# 
# # Bin size control + color palette
# ggplot(index_dt, aes(x=GEOID, index_shift, color=index_shift) ) +
#   geom_point(position='jitter') +
#   #geom_hex(bins = 70) +
#   scale_color_continuous(type = "viridis") +
#   theme_bw() 
# 
# # Bin size control + color palette
# ggplot(measure_raw[item=='Lead Risk From Housing'], aes(x=measure_old, measure_new, color=measure_shift) ) +
#   geom_point(position='jitter') +
#   #geom_hex(bins = 70) +
#   scale_color_continuous(type = "viridis") +
#   theme_bw() 
# 
# ggplot(measure_raw[item=='PM 2.5 -Diesel Emissions (Annual Tons/Km2)'], aes(x=measure_old, measure_new, color=measure_shift) ) +
#   geom_point(position='jitter') +
#   #geom_hex(bins = 70) +
#   scale_color_continuous(type = "viridis") +
#   theme_bw() 

#--------------------------------------------------------
# Purpose: re-create IBL ranks & QA observed vs expected ranks
# Author: Christopher Ahmed (christopher.ahmed@doh.wa.gov)
# Date: January 2022
#
# Note: this script reads in a data extract created from the backend
#   of WTN. The extract includes all the data IBL uses to create its ranks
#   the extract is a long data-set as opposed to wide.
#--------------------------------------------------------

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