# ----HEADER------------------------------------------------------------------------------------------------------------
# Author: JF
# Date: 05/31/2022
# Purpose: A shiny mapping and viztool for the mapping/sensitivity EHD proj
# source("/homes/jfrostad/_code/ehd_mapsense/baldR/app.R", echo=T)
#***********************************************************************************************************************

# ----CONFIG------------------------------------------------------------------------------------------------------------
# clear memory
rm(list=ls())

#set opts
options(scipen=999) #readability
#use cairo to render instead of quartz (quartz causes big slowdowns with geom_sf)
if(!identical(getOption("bitmapType"), "cairo") && isTRUE(capabilities()[["cairo"]])){
  options(bitmapType = "cairo")
}

## Set core_repo location
user            <- Sys.info()['user']
my_repo         <- ifelse(Sys.info()["sysname"] == "Linux",
                          file.path('/homes', user, ''),
                          file.path('C:/Users', user, 'Documents/ehd_mapsense'))

## Load libraries and  MBG project functions.
library(rsconnect)
library(ggplot2)
library(shiny)
library(shinyWidgets)
library(magrittr)
library(data.table)
library(tmap)
library(dplyr)

#capture date
today <- Sys.Date() %>% gsub("-", "_", .)

#options
run_version <- 'v1'
waLat <- 47.5 #lat for WA state
waLon <- -121.9 #long for WA state
waZoom <- 8 #starting zoom

global_data <- reactiveVal(NULL)
global_sf <- reactiveVal(NULL)
#***********************************************************************************************************************

# ----IN/OUT------------------------------------------------------------------------------------------------------------
###Input###
#TODO necessary section for shiny? Everything is pathed relative packaged within the directory of app.R
#raw data
main.dir <- file.path('C:/Users/jfrostad/Documents/ehd_mapsense/_code/baldR')
# data.dir <- file.path(main.dir, 'data')

# ###Output###
# out.dir  <- file.path(main.dir, 'output')
# viz.dir  <- file.path(main.dir, 'viz')
#***********************************************************************************************************************

# ---FUNCTIONS----------------------------------------------------------------------------------------------------------
##function lib##
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
  if(filter_geocodes %>% is.character) dt <- dt[!(geocode %in% filter_geocodes)]
  dt[, GEOID := as.character(geocode)] #setup merge var as chr to match shp
  
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
  
  return(plot)
  
}
#***********************************************************************************************************************

# ----PREP DATA------------------------------------------------------------------------------------------------------------
#read in relevant datasets
dataset <- readRDS('data/viz_data.RDS')

tract_sf <- dataset$tracts
places_sf <- dataset$places %>% 
  select(NAME) %>% 
  rename(PLACE=NAME)
water_sf <- dataset$water
road_sf <- dataset$roads
dt <- dataset$data 

#merge places onto tracts
# tract_sf <- st_join(tract_sf,  
#                     places_sf,
#                     join = st_nearest_feature, left=T)

#set values for the input slider
layer_vars <- 1:3

#merge dropouts to sf
data_sf <- tract_sf %>% 
  merge(dt[level==1], by='GEOID', allow.cartesian=T) 

#create a manual diverging color scale to make sure that the index shifts are uniformly depicted
div_colors <- RColorBrewer::brewer.pal(11, 'RdBu') %>% rev
names(div_colors) <- unique(dataset$ranks$measure_shift_capped) %>% sort
bin_colors <- c('#f1a340', '#f7f7f7', '#998ec3') %>% rev

#create a manual discrete color scale for the continuous index and make sure the color scale legend has all integers
rank_colors <- viridis::magma(n = 10)
names(rank_colors) <- 1:10

#create a manual cont color scale and have no labels
cont_colors <- viridis::magma(n = 12)

#***********************************************************************************************************************

# ---SERVER----------------------------------------------------------------------------------------------------------
##server function##
server <- function(input, output, session) {

  output$map <- renderTmap({
    
    #tmap_mode("view")

    tm_shape(shp = data_sf) +
      tm_view(set.view = c(waLon, waLat, waZoom)) +
      tm_polygons(col = "rank",
                  alpha = .9,
                  palette = rank_colors,
                  labels = names(rank_colors),
                  style='cat',
                  border.alpha = 0,
                  popup.vars=c("Measure"="item",  "Theme"="theme", 
                               #"County"="NAMELSADCO", 
                               #'Overall Rank'=index,
                               'Life Expectancy'='life_expectancy'),
                  id='name',
                  zindex=401) +
      tm_shape(shp=water_sf) +
      tm_polygons(col='gray91') + 
      tm_shape(shp=road_sf) +
      tm_lines(col='gray5')
    
  })
  
    #capture events that have to do with layer change  
    observeEvent(input$layer, {

      #capture change to dataset layer for varSelectInput for variables
      global_data(dataset$data[level==input$layer])
      
      data_sf <- tract_sf %>% 
        merge(global_data(), by='GEOID', allow.cartesian=T) 
      
      #also capture the change to the mapping layer
      global_sf(data_sf) #TODO there may be a better way using reactive({})
      
      # Make drop-down choice of measures dependent upon user input of layer
      freezeReactiveValue(input, "measure")
      updateSelectInput(session,
                        "measure",
                        choices = global_data()$item %>% unique %>% sort)
    })  
    
    #update varSelectInput options
    observeEvent(global_data(), {
      freezeReactiveValue(input, "var")
      updateVarSelectInput(session, "var", 
                           data = global_data(),
                           selected='rank')

    })

    #whenever the measure changes (last input) redraw the map
    #TODO or make it so that whenever measure or variable change?? because you could change var within measure
    toListen <- reactive({
      list(input$var,input$measure, input$alpha, input$scale)
    })
    
    observeEvent(toListen(), {

    #filter measure/theme if necessary
    data_sf <- global_sf() %>% dplyr::filter(item==input$measure)

    #setup scale
    if (input$var %like% 'rank') col_pal <- rank_colors
    else if (input$var %like% 'shift') col_pal <- div_colors
    else if (input$var %like% 'dropout') col_pal <- bin_colors
    else col_pal <- NULL
    if (input$var %like% 'rank|dropout') col_style <- 'cat'
    else col_style <- 'cont'

    #update map
    tmapProxy("map", session, {
      
      tm_remove_layer(401) +
        tm_shape(shp = data_sf) +
        tm_polygons(col = input$var %>% as.character, 
                    alpha = input$alpha %>% as.numeric,
                    palette = col_pal,
                    labels = names(col_pal),
                    style=col_style,
                    border.alpha = 0,
                    popup.vars=c("Measure"="item",  "Theme"="theme", 'Overall Rank'='rank',
                                 #"County"="NAMELSADCO", 
                                 'Life Expectancy'='life_expectancy'),
                    id='name',
                    zindex=401) 
      
    })
  })

}

#***********************************************************************************************************************
 
# --- UI ------------------------------------------------------------------------------------------------------------
##ui function ##

ui <- fluidPage(
  
  titlePanel("EHD Dataset Mapping Tool"),

  sidebarPanel(
    selectInput("layer", "Level", layer_vars,
                selected = 1),
    varSelectInput("var", "Variable", dt, 
                   selected = 'rank'),
    selectInput("measure", "Measure", 'Aggregated'),
    radioButtons("alpha",
                 label = "Transparency Level",
                 choices = c(.3, .5, .9),
                 selected = .9,
                 inline = TRUE),
    # switchInput("scale",
    #             label= 'Logarithmic Scale',
    #             value = FALSE)
  # 
  #   # sliderInput('sampleSize', 'Sample Size', min=1, max=nrow(dataset),
  #   #             value=min(1000, nrow(dataset)), step=500, round=0),
  # 
  #   selectInput('x', 'X', names(dataset)),
  #   selectInput('y', 'Y', names(dataset), names(dataset)[[2]]),
  #   selectInput('color', 'Color', c('None', names(dataset))),
  # 
  #   checkboxInput('jitter', 'Jitter'),
  #   checkboxInput('smooth', 'Smooth'),
  # 
  #   selectInput('facet_row', 'Facet Row', c(None='.', names(dataset))),
  #   selectInput('facet_col', 'Facet Column', c(None='.', names(dataset)))
  ),

  mainPanel(
    width = 12, height=8,
    tmapOutput('map', height="100vh")
  )
  
)

shinyApp(ui, server)

#***********************************************************************************************************************

# ---SCRAP DUMP----------------------------------------------------------------------------------------------------------
##SCRAP CODE
##
# ui <- fluidPage(
#   tmapOutput("map"),
#   selectInput("var", "Variable", world_vars)
# )
# 
# server <- function(input, output, session) {
#   output$map <- renderTmap({
#     tm_shape(World) +
#       tm_polygons(world_vars[1], zindex = 401)
#   })
#   
#   observe({
#     var <- input$var
#     tmapProxy("map", session, {
#       tm_remove_layer(401) +
#         tm_shape(World) +
#         tm_polygons(var, zindex = 401)
#     })
#   })
# }


# getDT <- function(query) {
#   URL = paste0('https://curtovmanagementapidev1.azurewebsites.net/Data/Run/Result?runId=', query$guid)
#   json <- fromJSON(URL)
#   dt <- as.data.table(json$results)
#   dt[, row_created := Sys.time()]
#   return(dt)
# }
# 
# observe({
#   query <- parseQueryString(session$clientData$url_search)
#   if (!is.null(query[['guid']])) {
#     updateTextInput(session, "guid", value = query[['guid']])
#   }
# })
# 
# results <- reactive({
#   query <- parseQueryString(session$clientData$url_search)
#   data.url <- paste0('https://curtovmanagementapidev1.azurewebsites.net/Data/Run/Result?runId=', query$guid)
#   
#   if (!is.null(query[['guid']])) {
#     getDT(query)
#   }
#   
# })


# dataset <- reactive({
#   diamonds[sample(nrow(diamonds), input$sampleSize),]
# })

# output$plot <- renderPlot({
# 
#   p <- ggplot(dataset(), aes_string(x=input$x, y=input$y)) + geom_point()
# 
#   if (input$color != 'None')
#     p <- p + aes_string(color=input$color)
# 
#   facets <- paste(input$facet_row, '~', input$facet_col)
#   if (facets != '. ~ .')
#     p <- p + facet_grid(facets)
# 
#   if (input$jitter)
#     p <- p + geom_jitter()
#   if (input$smooth)
#     p <- p + geom_smooth()
# 
#   print(p)
# 
# }, height=700)