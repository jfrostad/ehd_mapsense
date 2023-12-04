# ----HEADER------------------------------------------------------------------------------------------------------------
# Author: JF
# Date: 04/08/2022
# Purpose: Prepare EHD rankings from raw data, explore and analyze sensitivity
#***********************************************************************************************************************

# ----CONFIG------------------------------------------------------------------------------------------------------------
# clear memory
rm(list=ls())

#set opts
set.seed(98118)
options(scipen=999) #readability

#set control flow params
reload <- F #set true if you want to reprep all the data
analyze_sensitivity <- T #set true if you running sens analyses to compare against EHD map

## Set core_repo location
#TODO this is where you will change the filepaths to reflect the filestructure on your DOH PC
user            <- Sys.info()['user']
main.dir         <- ifelse(Sys.info()["sysname"] == "Linux",
                           file.path('/homes', user, '_code/ehd_mapsense/'),
                           file.path('/Users', user, 'Documents/work/ehd_mapsense/'))
my_repo <- file.path(main.dir, 'repo')
setwd(my_repo)

#TODO document which packages are used for which modules
pacman::p_load(readxl, janitor, data.table,stringr, magrittr, scales, Hmisc,
               ggplot2, viridis,
               sf,
               tigris, tidycensus, ggcorrplot,
               COINr, randtoolbox, sensobol, #sens packages
               stargazer,
               zoo)

#***********************************************************************************************************************

# ----IN/OUT------------------------------------------------------------------------------------------------------------
###Input###
#raw data
code.dir <- file.path(my_repo, 'code')
data.dir <- file.path(my_repo, 'data')
data_extract_EHDv1 <- 'ehd_data_v1_1.xlsx'
data_extract_EHDv2 <- 'ehd_data_v2.xlsx'

###Output###
out.dir <- file.path(main.dir, 'output')
viz.dir  <- file.path(main.dir, 'viz')
#***********************************************************************************************************************

# ---FUNCTIONS----------------------------------------------------------------------------------------------------------
#source custom functions that are relevant to this module
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

#merge the dropouts back onto the other tables for graphing
measure_dt <- merge(measure_dt, index_dt[, .(GEOID, dropout, index=rank, impacted_hierarchy,
                                             impacted, impacted_v1)], by='GEOID', all.x=T)
theme_dt <- merge(theme_dt, index_dt[, .(GEOID, dropout, index=rank, impacted_hierarchy,
                                         impacted, impacted_v1)], by='GEOID', all.x=T)

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
dt <- list(index_dt[, -c('county'), with=F],
           theme_dt,
           measure_dt) %>% 
  rbindlist(use.names=T, fill=T)
dt[item=='Aggregated', item_short := 'Agg'] #TODO fix earlier
dt <- dt[!(is.na(item_short))] #rows that didn't have data for v1

#save a  version of the data for the online mapping tool
out <- list(
  'dt'=dt,
  'tract_sf'=tract_sf,
  'water_sf'=water_sf,
  'road_sf'=road_sf,
  'places_sf'=places_sf
)

saveRDS(out, file=file.path(data.dir, 'prepped_data.RDS'))

} else file.path(data.dir, 'prepped_data.RDS') %>% readRDS %>% list2env(., globalenv())

#***********************************************************************************************************************

# ---VIZ/MAPPING OPTS---------------------------------------------------------------------------------------------------
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
#***********************************************************************************************************************

# ---COINr--------------------------------------------------------------------------------------------------------------
##
#this module is built on the R package COINr, which is documented below:
#https://bluefoxr.github.io/COINr/
##

#first, setup a COIN that is comparable to the baseline EHD v2.0
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

#meta
#this file contains the metadata for the coin and must follow the standard COINr paradigm
#its worth closely reading the documentation about metacoins in order to under the thinking here
#TODO i think that it would make sense to just save this file as a spreadsheet somewhere and update it whenever new
#indicators are added to the EHD. but i have kept this code inline so you can see how i made it in the first place
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

#build coin
coin_baseline <- new_coin(iData = ehd_coin_dt,
                 iMeta = ehd_coin_meta,
                 level_names = c("Indicator", "Pillar", "Sub-index", "Index"))

#normalize the data distribution by decile ranking
#note that this is a custom normalization function i wrote to help COINr emulate the EHD system
#stored in _lib/sens_fx.R
coin_baseline <- Normalise(coin_baseline, dset = "Raw",
                  global_specs=list(f_n='n_brank'))

#aggregate using base function to compare to our results
coin_baseline <- Aggregate(coin_baseline, dset = "Normalised", 
                  w='Original',
                  flatten_hierarchy = F,
                  #note that we take the arith mean until the last aggregation, where we multiply for risk score
                  f_ag = c("a_amean", "a_amean", "prod")) 

#extract results as a table 
coin_results_dt <- coin_baseline$Data$Aggregated %>% as.data.table
coin_results_dt[, ehd_decile_rank := n_brank(ehd_rank, nranks = 10)] #generate decile rankings

#make generic ehd map
cartographeR(dt=coin_results_dt, map_varname = 'ehd_decile_rank', map_label = 'EHD Index \n(Deciles)',
             geo_var='uCode',
             map_title = '',
             tag = 'coin_results', scale_type='cont_man', scale_vals = cont_colors,
             get_plot=T)
#***********************************************************************************************************************

# ---ALT COINs----------------------------------------------------------------------------------------------------------
#this section is to prepare any alternative formulations of the coins
#we will compare these to the baseline coin in our sensitivity analysis
if(analyze_sensitivity) {

#add in any new indicators for sensitivity analysis
potential_indicators_list <- c("asthma", 'wildfire_smoke', 'pesticides')

#read in the data for these new inds
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

#add the new inds
ehd_coin_meta <- ehd_coin_meta %>% 
  list(., ehd_coin_meta_new) %>% 
  rbindlist  

#now rebuild the coin with new info
#note that if you wanted to test different methods of COIN generation besides just changing indicators, you would
#do it in this section by changing the specs in the relevant step
#build coin
coin_sens <- new_coin(iData = ehd_coin_dt,
                      iMeta = ehd_coin_meta,
                      level_names = c("Indicator", "Pillar", "Sub-index", "Index"))

#normalize the data distribution by decile ranking
coin_sens <- Normalise(coin_sens, dset = "Raw",
                       global_specs=list(f_n='n_brank'))

#aggregate using base function to compare to our results
coin_sens <- Aggregate(coin_sens, dset = "Normalised", 
                       w='Original',
                       flatten_hierarchy = F,
                       #note that we take the arith mean until the last aggregation, where we multiply for risk score
                       f_ag = c("a_amean", "a_amean", "prod")) 
  
#adjustments/tests
#remove the new indicators
# remove new indicators and regenerate the coin to create a counterfactual
coin_cf <- change_ind(coin_sens, drop = potential_indicators_list, regen = TRUE)
#add/drop inds to test
coin_asthma <- change_ind(coin_cf, add = 'asthma', regen=T)
coin_smoke <- change_ind(coin_cf, add = 'wildfire_smoke', regen=T)
coin_pest <- change_ind(coin_cf, add = 'pesticides', regen=T)
coin_poc <- change_ind(coin_cf, drop = c("people_of_color"), regen = T) 

# remove two indicators and regenerate the coin
coin_list <- list(
  'coins'=list(coin_sens, coin_asthma, coin_smoke, coin_pest, coin_poc),
  'names'=list('ehd_rank', 'ehd_asthma', 'ehd_smoke', 'ehd_pest', 'ehd_white')
)

#helper function to extract data.tabled information from coins and standardize it for calculations
extractCoin <- function(i, coin_list) {
  
  #extract the values from list
  this_coin <- coin_list$coins[[i]]
  this_name <- coin_list$names[[i]]
  
  #make comparison table
  dt <- this_coin$Data$Aggregated %>% 
    as.data.table %>% 
    .[, .(uCode, ehd_rank, var=this_name)] %>% 
    #calculate deciles for comparison
    .[, ehd_rank := n_brank(ehd_rank, nranks = 10)] %>% 
    return
  
}

#run through all our comparison cases and create a results table
sens_dt <- lapply(1:length(coin_list$coins), extractCoin, coin_list=coin_list) %>% 
  rbindlist %>% 
  #cast wide for calculations
  dcast(uCode~var, value.var = 'ehd_rank') %>% 
  .[, asthma_diff := ehd_asthma - ehd_rank] %>%
  .[, smoke_diff := ehd_smoke - ehd_rank] %>%
  .[, pest_diff := ehd_pest - ehd_rank] %>%
  .[, white_diff := ehd_white - ehd_rank] %>%
  .[, level := 1] %>% 
  setnames('uCode', 'GEOID')
  
#TODO it would make sense to just turn this into an lapply loop for simplicity
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

#calculate and print relevant statistics
testChanges <- function(dt, this_col, thresh=8) {
  dt[, test_col := get(this_col)]
  message('modifying ', this_col, ' changes the ranking for:')
  dt[, sum(ehd_rank!=test_col)]  %>% print
  dt[, sum(ehd_rank!=test_col)/.N]  %>% print
  message('modifying ', this_col, ' changes the impact for:')
  dt[, sum((ehd_rank>=thresh)!=(test_col>=thresh))] %>% print
  dt[, sum((ehd_rank>=thresh)!=(test_col>=thresh))/.N]  %>% print
  dt[, test_col := NULL]
}

#TODO it would make sense to just turn this into an lapply loop for simplicity
testChanges(dt=sens_dt, this_col='ehd_white')
testChanges(dt=sens_dt, this_col='ehd_smoke')
testChanges(dt=sens_dt, this_col='ehd_pest')
testChanges(dt=sens_dt, this_col='ehd_asthma')

#output the PoC index so it can be used by the agency
file.path(out.dir, 'ehd_ranks_without_poc.csv') %>% 
  write.csv(sens_dt[, .(GEOID, ehd_rank_wo_poc=ehd_white)], file=.)

}
#***********************************************************************************************************************

# ---GSA SPECS----------------------------------------------------------------------------------------------------------
#dear future analysts:
#i kept this section in to demonstrate a variety of different ways that you can create alternative formulations of EHD
#within the COINr framework. some of them are more complicated than you will want to use and were only relevant to the
#global sensitivity analysis for my PhD but i think some of them could be relevant to EHD v3
#if you want to see more code that shows how to generate a variety of publication figures from the GSA results you can 
#look at my longer script for the gsa: ./frey.R

#now use our COIN to run a global sensitivity analysis across our parameters
#first we build out our parameters for each step

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
w_nom <- coin_baseline$Meta$Weights$Original

#build an equally weighted version
# copy original weights
w_equal <- copy(w_nom)
w_equal$Weight <- 1
coin_baseline$Meta$Weights$EqualWeights <- w_equal

#build a version weighted by the first comp of PCA
pca_results <- get_PCA(coin_baseline, dset = "Raw", Level = 1, out2 = "list", by_groups=F, weights_to='PCAComp1', impute=T)
coin_baseline$Meta$Weights$PComp1 <- pca_results$Weights
coin_baseline$Meta$Weights$PComp1$Weight <- abs(pca_results$Weights$Weight) #TODO what to do about negative loadings? arbitrary

#build a version that uses the inverted correlations of the ranked data in order to downweight highly correlated inds
#see spiegelhalter 2012
corr_results <- get_corr(coin_baseline, dset='Normalised', Levels = 1)  %>% 
  as.data.table() %>% 
  .[Var1!=Var2] %>% #drop the equivalent ones
  .[, Weight := 1/abs(sum(Correlation, na.rm = T)), by=Var1] %>% #TODO how to account for negative correlations?
  unique(by='Var1') %>% 
  .[, .(iCode=Var1, Level=1, Weight)]

coin_baseline$Meta$Weights$InvCorr <- list(
  corr_results,
  w_nom[Level>1] #use the same weights for higher levels
) %>% rbindlist

#drop indicators to proxy indicator selection
# get 100 replications
drop_cols <- names(coin_baseline$Data$Raw)[-1]
l_ind<- list(Address = "$Log$new_coin$exclude",
                Distribution = drop_cols,
                Type = "discrete")

#simulate random noise into the raw var for measurement error
# get 100 replications
noisy_dts <- get_noisy_dt(dt = coin_baseline$Data$Raw, noise_factor = .25, Nrep = 100)
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
 
  SA_res <- get_sensitivity(coin_baseline, SA_specs = SA_specs, N = 5000, SA_type = "SA", use_branks=F,
                            dset = "Aggregated", iCode = "ehd_rank", Nboot = 100, ncores = gsa_cores, sock=T)
  saveRDS(SA_res, file=file.path(scratch.dir, paste0('gsa_output_v', gsa_version, '.RDS')))

} else SA_res <- file.path(scratch.dir, paste0('gsa_output_v', gsa_version, '.RDS')) %>% readRDS
#***********************************************************************************************************************

#***********************************************************************************************************************

# ---SCRAP -------------------------------------------------------------------------------------------------------------
# hoard your scraps here 
# #***********************************************************************************************************************