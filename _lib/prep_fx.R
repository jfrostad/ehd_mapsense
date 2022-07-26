# ----HEADER------------------------------------------------------------------------------------------------------------
# Author: JF
# Date: 07/26/2022
# Purpose: Store functions used for data prep/processing/rank calculations
# source("/homes/jfrostad/_code/ehd_mapsense/_lib/prep_fx.R", echo=T)
#***********************************************************************************************************************

# ----FUNCTIONS---------------------------------------------------------------------------------------------------------
#custom function to build the EHD rankings, ported from chris' dplyr code (see scrap)
rankeR <- function(dir, path,
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
          measure_rank_val=rankcalculatedvalue,
          level=3)] #lowest level of aggregation
  
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
    .[, item := 'Aggregated'] %>%  #item no longer relevant after collapse
    .[, level := 2] %>%  #second level of aggregation
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
    .[, item := 'Aggregated'] %>%  #theme and item no longer relevant after collapse
    .[, theme := 'Aggregated'] %>%  #theme and item no longer relevant after collapse
    .[, level := 1] %>%  #highest level of aggregation
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

#***********************************************************************************************************************
#
# ---SCRAP -------------------------------------------------------------------------------------------------------------
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