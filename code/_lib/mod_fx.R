# ----HEADER------------------------------------------------------------------------------------------------------------
# Author: JF
# Date: 10/06/2022
# Purpose: Store functions used for modelling tasks
# source("/homes/jfrostad/_code/ehd_mapsense/_lib/mod_fx.R", echo=T)
#***********************************************************************************************************************

# ----FUNCTIONS---------------------------------------------------------------------------------------------------------
#Function to draw circle
#for biplots
circleR <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

#extract VIP scores from PLS regression
vipScoreR <- function(object) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}

#custom function to run pc type models (PCA/PLS) and generate output diagnostics and plots
pcWrappeR <- function(dt, pred_dt=NULL, var, #general opts
                      strat_var=NULL, strat_val=NULL, #var/val to stratify over 
                      #TODO build out stratification framework
                      depvar=NA, ncomps=3, #pls opts - note that i manually tested PLS to decide 3comp by elbow
                      hc_k=7, #determined visually, number of clusters for HCA
                      theme_map=theme_dt,
                      tag='') {
  
  #if there is a stratification var, we subset the data and note it in the tag
 if(strat_var %>% is.character) { 
   message('stratifying data on ', strat_var, '\n=', strat_val)
   dt <- dt[get(strat_var)==strat_val]
   tag <- paste0(tag, '_', strat_var, '_', strat_val)
 }
  
  #if there isnt a depvar, we are running PCA
  if(depvar %>% is.na) {
    message('running PCA')
    #melt data and then run the PCA
    mod_dt <- dt[, c('GEOID', 'item_short', var), with=F] %>% 
      dcast(GEOID~item_short, value.var=var)
    
    mod <- na.omit(mod_dt) %>% 
      .[, -c('GEOID'), with=F] %>% 
      prcomp(scale=T)
    
    #also run HCA 
    #TODO should scale it if we decide to use the raw variables at any point in this analysis
    message('running HCA')
    d <- dist(mod_dt%>% na.omit, method = "euclidean")
    
    # methods to assess
    m <- c( "average", "single", "complete", "ward")
    names(m) <- c( "average", "single", "complete", "ward")
    
    # function to compute coefficient
    ac <- function(x) {
      cluster::agnes(mod_dt, method = x)$ac
    }
    
    #test all functions
    purrr::map_dbl(m, ac) %>% print
    
    #here we see that ward identified the most clustering
    hc3 <- agnes(mod_dt %>% na.omit, method = "ward")
    pltree(hc3, cex = 0.6, hang = -1, main = "Dendrogram of agnes") 
    
    hc4 <- hclust(d, method = "ward.D2" )
    
    #cut into 6-8 clusters based on visual inspections
    hc_clusters <-
      hc4 %>% 
      cutree(k = hc_k)
      #.[hc4$order] #TODO troubleshoot the ordering
    
    fviz_cluster(object = list(data = mod_dt[,-1] %>% na.omit,
                               cluster = hc_clusters),
                 main = "Ward's Method Clusters")
    
    #give us the means for each var per cluster (centroids)
    hc_dt <- 
      mod_dt %>% na.omit %>% 
      .[, cluster := hc_clusters]
    hc_centroids <- hc_dt[, lapply(.SD, mean), 
        .SDcols=names(mod_dt)[-1], by=cluster] %>% 
      melt(id.vars='cluster', variable.name='variable') %>% 
      merge(theme_map,
            by.x='variable',
            by.y = 'item_short',
            all.x=T) %>% 
      .[, strat_var := strat_val]


    ggplot(data=hc_centroids, aes(x=item_themes, y=value, color=cluster, group=cluster)) +
      geom_line() + 
      scale_color_viridis(option='viridis') +
      theme_minimal() +
      theme(axis.text.x=element_text(angle=45, vjust=0.5))
    file.path(viz.dir, paste0('ggclusters', tag, '.png')) %>% ggsave(height=8, width=12)
    
    hc4 %>% 
      as.dendrogram %>% 
      ggdendrogram(rotate = TRUE)
    file.path(viz.dir, paste0('ggdend', tag, '.png')) %>% ggsave(height=8, width=12)

    #use hclust to plot more nicely
    # plot(as.phylo(hc3_dend), type = "cladogram", cex = 0.6, 
    #      label.offset = 0.5)
    
    } else {
    #if there is a depvar, we are running a pls  
    #melt data and then convert to matrices
    mod_dt <- dt[, c('GEOID', 'item_short', var, depvar), with=F] %>% 
      dcast(...~item_short, value.var=var) %>% 
      .[, -c('GEOID'), with=F] %>% 
      na.omit
    y_mat <- as.matrix(mod_dt[,1])
    x_mat <- as.matrix(mod_dt[,2:ncol(mod_dt)])
    mod <- pls::mvr(y_mat  ~ x_mat , ncomp=ncomps, method = "oscorespls" , scale = T)
    
  }
  
  #create diagnostic plots
  message('creating diganostic plots')
  if (mod %>% class == 'prcomp') {
    fviz_eig(mod) + ggtitle(paste0('PCA Screeplot ', strat_var, ' ', strat_val))
    file.path(viz.dir, paste0('pca_eig', tag, '.png')) %>% ggsave(height=8, width=12)
    fviz_contrib(mod, choice='var', axes=1)  + ggtitle(paste0('PC #1 Loading ', strat_var, ' ', strat_val))
    file.path(viz.dir, paste0('pca_contrib', tag, 'v1.png')) %>% ggsave(height=8, width=12)
    fviz_contrib(mod, choice='var', axes=2)  + ggtitle(paste0('PC #2 Loading ', strat_var, ' ', strat_val))
    file.path(viz.dir, paste0('pca_contrib', tag, 'v2.png')) %>% ggsave(height=8, width=12)
    fviz_contrib(mod, choice='var', axes=3)  + ggtitle(paste0('PC #3 Loading ', strat_var, ' ', strat_val))
    file.path(viz.dir, paste0('pca_contrib', tag, 'v3.png')) %>% ggsave(height=8, width=12)
    
    fviz_pca_var(mod,
                 col.var = "contrib", # Color by the quality of representation
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE     # Avoid text overlapping
    )   + ggtitle(paste0('PCA Loadings ', strat_var, ' ', strat_val))
    file.path(viz.dir, paste0('pc_var', tag, '.png')) %>% ggsave(height=8, width=12)

    #merge on the impacts to plot
    biplot_dt <- dt[, .(GEOID, impacted_hierarchy)] %>% 
      unique(by='GEOID') %>% 
      merge(mod_dt %>% na.omit, by='GEOID')
    
    fviz_pca_biplot(mod,
                    label='none',
                    col.var='black',
                    habillage = biplot_dt$impacted_hierarchy, # Color by the quality of representation
                    repel = TRUE     # Avoid text overlapping
    )   + 
      ggtitle(paste0('PCA Biplot ', strat_var, ' ', strat_val)) +
      scale_shape_manual(values=c(16, 16, 16, 16), guide=F) +
      scale_color_brewer('Impact Status', palette = 'Paired') +
      theme_minimal()
    
    file.path(viz.dir, paste0('pc_biplot', tag, '.png')) %>% ggsave(height=8, width=12)
    
    fviz_pca_ind(mod,
                    label='none',
                    col.var='black',
                    habillage = biplot_dt$impacted_hierarchy, # Color by the quality of representation
                    repel = TRUE     # Avoid text overlapping
    ) +
      ggtitle(paste0('PCA Individuals ', strat_var, ' ', strat_val)) +
      scale_shape_manual(values=c(16, 16, 16, 16), guide=F) +
      scale_color_brewer('Impact Status', palette = 'Paired') +
      theme_minimal()
    
    file.path(viz.dir, paste0('pc_ind', tag, '.png')) %>% ggsave(height=8, width=12)

    #make predictions of pca #1/2/3s so that we can map them 
    if(pred_dt %>% is.null) pred_dt <- mod_dt %>% copy
    if(strat_var %>% is.character) pred_dt <- pred_dt[get(strat_var)==strat_val]
    pred_dt <- pred_dt %>% 
      cbind(., predict(mod, newdata=.)[,1:9]) %>% 
      .[,level:=3] #remind data that its most granular
    
    #map the first 3 elements of PCA
    message('mapping components')
    cartographeR(dt=pred_dt, map_varname = 'PC1', map_label = 'PCA #1',
                 map_title = paste0('Mapping the first component of PCA using the ', var,
                                    strat_var, ' ', strat_val),
                 scale_type='cont_grad', 
                 tag=tag,
                 lvl=3)
    cartographeR(dt=pred_dt, map_varname = 'PC2', map_label = 'PCA #2',
                 map_title = paste0('Mapping the second component of PCA using the ', var,
                                    strat_var, ' ', strat_val),
                 scale_type='cont_grad', 
                 tag=tag,
                 lvl=3)
    cartographeR(dt=pred_dt, map_varname = 'PC3', map_label = 'PCA #3',
                 map_title = paste0('Mapping the third component of PCA using the ', var,
                                    strat_var, ' ', strat_val),
                 scale_type='cont_grad', 
                 tag=tag,
                 lvl=3)
    
    #map the clusters
    cartographeR(dt=hc_dt[, level:=3], map_varname = 'cluster', map_label = 'HCA Cluster',
                 map_title = paste0('Mapping the HCA Clusters based on', var,
                                    strat_var, ' ', strat_val),
                 scale_type='cont_vir', 
                 tag=tag,
                 lvl=3)
    
  } else if (mod %>% class == 'mvr') {
    
    #extract objects from the PLS mod in order to create a biplot
    S_plsr <- scores(mod)[, comps= 1:2, drop = FALSE]
    cl_plsr <- cor(model.matrix(mod), S_plsr)
    df_cor <- as.data.frame(cl_plsr)
    
    #extract the correlation between the yvar and the components
    #TODO how to interpet this value?
    df_depend_cor <- as.data.frame(cor(mod_dt[,1], S_plsr))
    
    #generate biplot dt
    biplot_dt  <-  rbind(df_cor ,
                         df_depend_cor) %>% 
      as.data.table(keep.rownames=T) %>% 
      setNames(c('variable', "PC1", "PC2")) %>% 
      merge(theme_map,
            by.x='variable',
            by.y = 'item_short',
            all.x=T) %>% 
      .[variable=='overall', theme:= 'The Overall Ranking']

    #create custom biplot
    ggplot(data=biplot_dt , aes(x=PC1, y=PC2, color=theme))+
      geom_hline(aes(yintercept = 0), size=.2, linetype = 3)+ 
      geom_vline(aes(xintercept = 0), size=.2, linetype = 3)+
      geom_text_repel(aes(label=variable))+
      geom_path(data=circleR(c(0,0),2,npoints = 100),
                aes(x,y), colour = "darkgrey")+
      geom_segment(aes(x=0, y=0, xend=PC1, yend=PC2), 
                   arrow=arrow(length=unit(0.2,"cm")), alpha=0.75) +
      scale_color_brewer('Theme', palette='Paired') +
      scale_x_continuous(breaks = c(-1,-0.5,0,0.5,1))+
      scale_y_continuous(breaks = c(-1,-0.5,0,0.5,1))+
      coord_fixed(ylim=c(-1, 1),xlim=c(-1, 1))+
      xlab("PC1")+ 
      ylab("PC2")+ 
      ggtitle("PLS Biplot")+
      theme_bw() +
      theme(axis.line.x = element_line(color="darkgrey"),
                            axis.line.y = element_line(color="darkgrey"))+
      theme(axis.ticks = element_line(colour = "black"))+
      theme(axis.title = element_text(colour = "black"))+
      theme(axis.text = element_text(color="black"))+
      theme(panel.grid.minor = element_blank())
      
    file.path(viz.dir, 'pls_biplot.png') %>% ggsave(height=8, width=12)
    
    #extract the VIP scores
    #TODO how to interpret this value
    vip_dt  <- mod %>% vipScoreR %>% as.data.frame
    
    #extract the coefficients (PLS loadings)
    coef_dt <- as.data.frame(coef(mod, ncomp =1:3)) %>% 
      setnames(., c('PC1', 'PC2', 'PC3')) %>% 
      dplyr::mutate(variable = rownames(.)) %>% 
      as.data.table %>% 
      .[, variable := stringr::str_replace_all(variable, stringr::fixed("\\"), '')]  %>% 
      .[, variable := stringr::str_replace_all(variable, stringr::fixed("``"), '')] %>% 
      #merge on themes for plotting colors
      merge(theme_map,
            by.x='variable',
            by.y = 'item_short')  %>% 
      #melt long for plotting
      melt(coef_dt, 
           id.vars=c('variable', 'theme'),
           variable.name = 'component')

    #create plots of the PLS loadings
    #PC1
    ggplot(coef_dt_long, aes(x = forcats::fct_reorder(variable, value), y = value, group = 1, fill=theme))+  
      geom_bar_pattern(aes(pattern=component),
                       stat = "identity", position='stack',
                       color = "black", 
                       pattern_fill = "black",
                       pattern_angle = 45,
                       pattern_density = 0.1,
                       pattern_spacing = 0.025,
                       pattern_key_scale_factor = 0.6)+
      theme(axis.text.x = element_text(angle=65,
                                       hjust=1,
                                       size = 8),
            axis.title.y = element_text(size = 2))+
      theme_bw()+
      scale_fill_brewer('Theme', palette='Paired') +
      scale_pattern_manual(values = c(PC1 = "none", PC2 = "stripe", PC3='wave')) +
      coord_flip()
    file.path(viz.dir, 'pls_loadings.png') %>% ggsave(height=8, width=12)

    #make predictions of pca #1/2/3s so that we can map them 
    if(pred_dt %>% is.null) pred_dt <- mod_dt %>% copy
    
    #extract the component values for each geoid from the mod obj
    pc_dt <- mod$Yscores %>% 
      as.matrix %>% 
      .[,1:ncomps] %>% 
      as.data.table %>% 
      setnames(c('PLS_C1', 'PLS_C2', 'PLS_C3'))
    
    #bind those values to the GEOID dt 
    #TODO how to verify there is no resorting happening here
    #TODO could test the covariate/yvar values to verify
    pc_dt <- pred_dt[, .(GEOID, overall)] %>% 
      cbind(., pc_dt) %>% 
      .[, level:=3]

    #map the first 3 elements of PLS
    cartographeR(dt=pc_dt, map_varname = 'PLS_C1', map_label = 'PLS #1',
                 map_title = 'Mapping the first component of PLS using the ranked measures',
                 scale_type='cont_grad', 
                 tag=tag,
                 lvl=3)
    cartographeR(dt=pc_dt, map_varname = 'PLS_C2', map_label = 'PLS #2',
                 map_title = 'Mapping the second component of PLS using the ranked measures',
                 scale_type='cont_grad', 
                 tag=tag,
                 lvl=3)
    cartographeR(dt=pc_dt, map_varname = 'PLS_C3', map_label = 'PLS #3',
                 map_title = 'Mapping the third component of PLS using the ranked measures',
                 scale_type='cont_grad', 
                 tag=tag,
                 lvl=3)
    
  }
  
  list('mod'=mod,
       'preds'=pred_dt,
       'biplot_dt'=biplot_dt,
       'hc_centroids'=hc_centroids) %>% 
    return
  
}
