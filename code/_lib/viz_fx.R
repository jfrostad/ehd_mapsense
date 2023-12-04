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
                         geo_var='GEOID', #which variable contains the data we want to map from? standard form is GEOID
                         map_varname, map_label=NA, map_title=NA, tag=NA, #use tag to label files
                         subset_var=NA, subset_val=F,
                         facet_var=F,
                         filter_geocodes=c('53057990100'), #san juan water area with lots of big changes
                         scale_type='cont',
                         scale_vals=NULL,
                         get_plot=F) {

  #filter long dataset to appropriate level
  if('level' %in% names(dt)) dt <- dt[level==lvl]

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
  
  #standardize the geography variable if not standardized already
  if(geo_var!='GEOID') dt[, GEOID := get(geo_var)]
  
  #cleanup missing data
  #remove the -1s which represent nonmissing NAs (plot as NA)
  #but note that there can be negatives for the change vars
  if(map_varname %like% 'measure') dt <- dt[map_var>=0]

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
  else if(scale_type=='class') plot <- plot + scale_fill_manual(map_label, values=c('#4daf4a', '#984ea3', '#ff7f00', '#377eb8'))
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
fviz_pca_ind2 <- function (X, axes = c(1, 2), geom = c("point", "text"), geom.ind = geom, 
          repel = FALSE, habillage = "none", palette = NULL, addEllipses = FALSE, 
          col.ind = "black", fill.ind = "white", col.ind.sup = "blue", 
          alpha.ind = 1, select.ind = list(name = NULL, cos2 = NULL, 
                                           contrib = NULL), ...) 
{
  fviz(X, element = "ind", axes = axes, geom = geom.ind, habillage = habillage, 
       palette = palette, addEllipses = addEllipses, color = col.ind, 
       fill = fill.ind, alpha = alpha.ind, col.row.sup = col.ind.sup, 
       select = select.ind, repel = repel, ...)
}

fviz2 <- function (X, element, axes = c(1, 2), geom = "auto", label = "all", 
          invisible = "none", labelsize = 4, pointsize = 1.5, pointshape = 19, 
          arrowsize = 0.5, habillage = "none", addEllipses = FALSE, 
          ellipse.level = 0.95, ellipse.type = "norm", ellipse.alpha = 0.1, 
          mean.point = TRUE, color = "black", fill = "white", alpha = 1, 
          gradient.cols = NULL, col.row.sup = "darkblue", col.col.sup = "darkred", 
          select = list(name = NULL, cos2 = NULL, contrib = NULL), 
          title = NULL, axes.linetype = "dashed", repel = FALSE, col.circle = "grey70", 
          circlesize = 0.5, ggtheme = theme_minimal(), ggp = NULL, 
          font.family = "", ...) 
{
  .check_axes(axes, .length = 2)
  facto.class <- .get_facto_class(X)
  extra_args <- list(...)
  if (!is.null(extra_args$jitter)) 
    repel <- .facto_dep("jitter", "repel", TRUE)
  lab <- .label(label)
  hide <- .hide(invisible)
  if (is.null(title)) {
    element_desc <- list(ind = "Individuals", var = "Variables", 
                         col = "Column points", row = "Row points", mca.cor = "Variables", 
                         quanti.sup = "Quantitative variables", quanti.var = "Quantitative variables", 
                         quali.var = "Qualitative variable categories", group = "Variable groups", 
                         partial.axes = "Partial axes")
    if (facto.class == "MCA") 
      element_desc$var <- "Variable categories"
    title <- paste0(element_desc[[element]], " - ", facto.class)
  }
  if (geom[1] == "auto") {
    geom <- c("point", "text")
    if (element == "var" & facto.class == "PCA") 
      geom <- c("arrow", "text")
  }
  if (facto.class %in% c("CA", "MCA")) {
    if (element %in% c("row", "ind") & missing(color)) 
      color = "blue"
    else if (element %in% c("col", "var", "mca.cor") & missing(color)) 
      color = "red"
  }
  summary.res <- c("coord", "contrib", "cos2")
  if (element == "partial.axes" | (element == "quali.var" & 
                                   facto.class == "HMFA")) 
    summary.res <- c("coord", "contrib")
  else if (element == "group" & facto.class == "HMFA") 
    summary.res <- "coord"
  df <- facto_summarize(X, element = element, axes = axes, 
                        result = summary.res)
  colnames(df)[2:3] <- c("x", "y")
  is_grouping_var_exists <- !("none" %in% habillage) | .is_grouping_var(color) | 
    .is_grouping_var(fill)
  if (!("none" %in% habillage)) {
    dd <- .add_ind_groups(X, df, habillage)
    df <- dd$ind
    color <- dd$name.quali
    if (missing(pointshape)) 
      pointshape <- dd$name.quali
  }
  if (length(color) > 1) {
    if (nrow(df) != length(color)) 
      stop("The length of color variable", "should be the same as the number of rows in the data.")
    .col.name <- "Col."
    df[[.col.name]] <- color
    if (missing(pointshape) & .is_grouping_var(color)) 
      pointshape <- .col.name
    color <- .col.name
  }
  if (length(fill) > 1) {
    if (nrow(df) != length(fill)) 
      stop("The length of fill variable", "should be the same as the number of rows in the data.")
    .col.name <- "Fill."
    df[[.col.name]] <- fill
    fill <- .col.name
  }
  if (length(pointsize) > 1) {
    if (nrow(df) != length(pointsize)) 
      pointsize <- 1.5
    df[["pointsize"]] <- pointsize
    pointsize <- "pointsize"
  }
  df.all <- df
  if (!is.null(select)) 
    df <- .select(df, select)
  if (facto.class == "PCA" & element == "var" & !is.null(extra_args$scale.)) 
    df[, c("x", "y")] <- df[, c("x", "y")] * extra_args$scale.
  if (facto.class %in% c("CA", "MCA") & !(element %in% c("mca.cor", 
                                                         "quanti.sup"))) {
    if (!is.null(extra_args$map)) 
      df <- .scale_ca(df, res.ca = X, element = element, 
                      type = extra_args$map, axes = axes)
  }
  is.pca.var <- element == "var" & facto.class == "PCA"
  point <- ("point" %in% geom) & (!hide[[element]])
  if (missing(mean.point)) 
    mean.point <- (is_grouping_var_exists & !is.pca.var) & 
    ("point" %in% geom) & (!hide[["quali"]])
  if (element == "quanti.var") 
    mean.point <- FALSE
  label <- NULL
  if (lab[[element]] & "text" %in% geom & !hide[[element]]) 
    label <- "name"
  p <- ggplot()
  if (hide[[element]]) {
    if (is.null(ggp)) 
      p <- ggplot() + geom_blank(data = df, aes_string("x", 
                                                       "y"))
    else p <- ggp
  }
  else p <- ggpubr::ggscatter(data = df, x = "x", y = "y", 
                              color = color, fill = fill, alpha = alpha, shape = pointshape, 
                              point = point, size = pointsize, mean.point = mean.point, 
                              label = label, font.label = labelsize * 3, repel = repel, 
                              ellipse = addEllipses, ellipse.type = ellipse.type, ellipse.alpha = ellipse.alpha, 
                              ellipse.level = ellipse.level, main = title, ggtheme = ggtheme, 
                              ggp = ggp, font.family = font.family, ...)
  browser()
  if (alpha %in% c("cos2", "contrib", "coord", "x", "y")) 
    p <- p + scale_alpha(limits = range(df.all[, alpha]))
  if (!is.null(gradient.cols)) 
    p <- p + ggpubr::gradient_color(gradient.cols)
  if (is.null(extra_args$legend)) 
    p <- p + theme(legend.position = "right")
  if ("arrow" %in% geom & !hide[[element]]) 
    p <- p + .arrows(data = df, color = color, alpha = alpha, 
                     size = arrowsize)
  if (facto.class == "PCA" & element == "var") {
    if (.get_scale_unit(X) & is.null(extra_args$scale.)) 
      p <- .add_corr_circle(p, color = col.circle, size = circlesize)
  }
  else if (facto.class %in% c("MCA", "MFA", "HMFA", "FAMD") & 
           element %in% c("quanti.sup", "quanti.var", "partial.axes")) {
    p <- .add_corr_circle(p, color = col.circle, size = circlesize)
  }
  if ("facet_vars" %in% colnames(df)) {
    groups <- c("facet_vars", "Groups")
    xx <- ggpubr::desc_statby(df, measure.var = "x", grps = groups)[, 
                                                                    c(groups, "mean")]
    colnames(xx)[ncol(xx)] <- "x"
    yy <- ggpubr::desc_statby(df, measure.var = "y", grps = groups)[, 
                                                                    c(groups, "mean")]
    xx$y <- yy$mean
    grp_coord <- xx
    p <- p + ggpubr::geom_exec(geom_text, data = grp_coord, 
                               x = "x", y = "y", label = "Groups", color = color)
    p <- p + facet_wrap(~facet_vars) + theme(legend.position = "none")
  }
  scale. <- ifelse(is.null(extra_args$scale.), 1, extra_args$scale.)
  esup <- .define_element_sup(X, element, geom = geom, lab = lab, 
                              hide = hide, col.row.sup = col.row.sup, col.col.sup = col.col.sup, 
                              ...)
  ca_map = extra_args$map
  if (element == "mca.cor") 
    ca_map = NULL
  if (!is.null(esup)) 
    p <- .add_supp(p, X, element = esup$name, axes = axes, 
                   select = select, geom = geom, color = esup$color, 
                   shape = esup$shape, pointsize = pointsize, labelsize = labelsize, 
                   addlabel = esup$addlabel, repel = repel, linetype = 2, 
                   scale. = scale., ca_map = ca_map, font.family = font.family)
  p <- .fviz_finish(p, X, axes, axes.linetype, ...) + labs(title = title)
  p
}


#***********************************************************************************************************************