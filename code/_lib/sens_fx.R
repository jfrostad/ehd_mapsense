# ----HEADER------------------------------------------------------------------------------------------------------------
# Author: JF
# Date: 10/06/2022
# Purpose: Store functions used for modelling tasks
# source("/homes/jfrostad/_code/ehd_mapsense/_lib/mod_fx.R", echo=T)
#***********************************************************************************************************************

# ----FUNCTIONS---------------------------------------------------------------------------------------------------------

#pulled from COINr github
#' Sensitivity and uncertainty analysis of a coin
#'
#' This function performs global sensitivity and uncertainty analysis of a coin. You must specify which
#' parameters of the coin to vary, and the alternatives/distributions for those parameters.
#'
#' COINr implements a flexible variance-based global sensitivity analysis approach, which allows almost any assumption
#' to be varied, as long as the distribution of alternative values can be described. Variance-based "sensitivity indices"
#' are estimated using a Monte Carlo design (running the composite indicator many times with a particular combination of
#' input values). This follows the methodology described in \doi{10.1111/j.1467-985X.2005.00350.x}.
#'
#' To understand how this function works, please see `vignette("sensitivity")`. Here, we briefly recap the main input
#' arguments.
#'
#' First, you can select whether to run an uncertainty analysis `SA_type = "UA"` or sensitivity analysis `SA_type = "SA"`.
#' The number of replications (regenerations of the coin) is specified by `N`. Keep in mind that the *total* number of
#' replications is `N` for an uncertainty analysis but is `N*(d + 2)` for a sensitivity analysis due to the experimental
#' design used.
#'
#' To run either types of analysis, you must specify *which* parts of the coin to vary and *what the distributions/alternatives are*
#' This is done using `SA_specs`, a structured list. See `vignette("sensitivity")` for details and examples.
#'
#' You also need to specify the target of the sensitivity analysis. This should be an indicator or aggregate that can be
#' found in one of the data sets of the coin, and is specified using the `dset` and `iCode` arguments.
#'
#' Finally, if `SA_type = "SA"`, it is advisable to set `Nboot` to e.g. 100 or more, which is the number of bootstrap samples
#' to take when estimating confidence intervals on sensitivity indices. This does *not* perform extra regenerations of the
#' coin, so setting this to a higher number shouldn't have much impact on computational time.
#'
#' This function replaces the now-defunct `sensitivity()` from COINr < v1.0.
#'
#' @param coin A coin
#' @param SA_specs Specifications of the input uncertainties
#' @param N The number of regenerations
#' @param SA_type The type of analysis to run. `"UA"` runs an uncertainty analysis. `"SA"` runs a sensitivity
#' analysis (which anyway includes an uncertainty analysis).
#' @param dset The data set to extract the target variable from (passed to [get_data()]).
#' @param iCode The variable within `dset` to use as the target variable (passed to [get_data()]).
#' @param quietly Set to `TRUE` to suppress progress messages.
#' @param Nboot Number of bootstrap samples to take when estimating confidence intervals on sensitivity
#' indices.
#' @param check_addresses Logical: if `FALSE` skips the check of the validity of the parameter addresses. Default `TRUE`,
#' but useful to set to `FALSE` if running this e.g. in a Rmd document (because may require user input).
#'
#' @importFrom stats runif
#'
#' @return Sensitivity analysis results as a list, containing:
#' * `.$Scores` a data frame with a row for each unit, and columns are the scores for each replication.
#' * `.$Ranks` as `.$Scores` but for unit ranks
#' * `.$RankStats` summary statistics for ranks of each unit
#' * `.$Para` a list containing parameter values for each run
#' * `.$Nominal` the nominal scores and ranks of each unit (i.e. from the original COIN)
#' * `.$Sensitivity` (only if `SA_type = "SA"`) sensitivity indices for each parameter. Also confidence intervals if `Nboot`
#' was specified.
#' * Some information on the time elapsed, average time, and the parameters perturbed.
#' * Depending on the setting of `store_results`, may also contain a list of Methods or a list of COINs for each replication.
#'
#' @export
#'
#' @examples
#' # for examples, see `vignette("sensitivity")`
#' # (this is because package examples are run automatically and this function can
#' # take a few minutes to run at realistic settings)
#'
get_sensitivity <- function(coin, SA_specs, N, SA_type = "UA", dset, iCode, Nboot = NULL, quietly = FALSE,
                            use_branks=F, #give option to use the binned ranks like EHD
                            ncores=1, sock=F,
                            check_addresses = TRUE){
  
  t0 <- proc.time()
  # CHECKS ------------------------------------------------------------------
  
  # check_coin_input(coin)
  # stopifnot(is.list(SA_specs),
  #           is.numeric(N),
  #           length(N) == 1,
  #           N > 2,
  #           SA_type %in% c("SA", "UA"))
  # 
  # # check format of SA_specs
  # check_specs <- sapply(SA_specs, function(li){
  #   !is.null(li$Address) & !is.null(li$Distribution) & !is.null(li$Type)
  # })
  # 
  # if(any(!check_specs)){
  #   stop("One or more entries in SA_specs is missing either the $Name, $Address or $Distribution entries.")
  # }
  
  
  # PREP --------------------------------------------------------------------
  
  # number of uncertain input paras
  d <- length(SA_specs)
  
  message('beginning sobol sampling process')
  
  # get sample
  if(SA_type == "UA"){
    
    # a random (uniform) sample
    XX <- matrix(stats::runif(d*N), nrow = N, ncol = d)
    
  } else {
    
    if(d==1){
      stop("Only one uncertain input defined. It is not meaningful to run a sensitivity analysis
      with only one input variable. Consider changing SA_type to \"UA\".")
    }
    
    # use standard MC estimators of sensitivity indices
    #XX <- SA_sample(N, d)
    XX <- sobol_matrices(N=N, params=names(SA_specs), order='second')
    
  }
  # covert to df
  XX <- as.data.frame(XX)
  # total number of regens
  NT <- nrow(XX)
  
  # convert sample to parameters (data frame with list cols?)
  XX_p <- mapply(function(x, spec){
    sample_2_para(x, distribution = spec$Distribution, dist_type = spec$Type)
  }, XX, SA_specs, SIMPLIFY = FALSE)
  # name list according to parameters
  names(XX_p) <- names(SA_specs)
  
  # also get addresses
  addresses <- sapply(SA_specs, `[[`, "Address")
  
  # check addresses for validity
  if(check_addresses){
    a_check <- lapply(addresses, check_address, coin)
  }
  
  # RUN COINS ---------------------------------------------------------------
  
  # at this point the parameters are stored in a list where each entry of the list is a parameter,
  # and the entry contains N instances of each parameter
  
  # first get nominal results
  SA_scores <- get_data(coin, dset = dset, iCodes = iCode) %>% 
    as.data.table %>% 
    setnames(iCode, 'Nominal') %>% 
    setkey(uCode) %>% 
    .[, Nominal_rank := frank(Nominal, 
                              na.last='keep', 
                              ties.method = 'min')] %>% 
    .[, bin_size := floor(sum(!is.na(Nominal))/10)] %>% 
    #bin the ranks and number them by rounding up
    .[, Nominal_brank := (Nominal_rank / bin_size) %>% ceiling] %>% 
    .[Nominal_brank>10, Nominal_brank := 10] %>% 
    .[, Nominal_impacted := 0] %>% 
    .[Nominal_brank>=8, Nominal_impacted := 1]

  # make a df of NAs in case a coin regen fails
  #TODO but also i dont believe in failure
  # v_fail <- SA_scores
  # v_fail[names(v_fail) == iCode] <- NA
  

  # looping over each replication in the SA
  regenLoop <- function(irep) {
    
    # list of parameters for current rep
    l_para_rep <- lapply(XX_p, `[[`, irep)
    
    if (!quietly){
      message(paste0("Rep ",irep," of ",NT," ... ", round(irep*100/NT,1), "% complete" ))
    }
    
    # regenerate coin using parameter list
    coin_rep <- regen_edit(l_para_rep, addresses, coin)

    # extract variable of interest
    if(is.coin(coin_rep)){
      v_out <- get_data(coin_rep, dset = dset, iCodes = iCode) %>% 
        as.data.table %>% 
        .[, rep := irep] %>% 
        setkey(uCode)
    } else {
      # df with just NAs
      v_out <- v_fail
    }
    
  return(v_out)
    
  }

  #loop over all runs in parallel and then merge the results
  message('beginning simulation')
  if(ncores>1 & sock!=F) v_out <- pbmclapply(1:NT, regenLoop, mc.cores=ncores) %>% rbindlist
  else if (ncores>1 & sock==T) {
    #prep cluster
    cl <- makeCluster(ncores)
    clusterExport(cl, c('XX_p', 'NT', 'l_para_rep', 'addresses', 'coin', 'dset'))
    #run
    v_out <- parLapply(cl, 1:5000, regenLoop)
    #cleanup
    stopCluster(cl)
    
  } else  v_out <- lapply(1:NT, regenLoop) %>% rbindlist

  # merge onto nominal results
  SA_scores <- merge(SA_scores, v_out, by = "uCode", all = TRUE) %>% 
    setnames(iCode, 'score')
  
  
  # POST --------------------------------------------------------------------
  
  message('calculating uncertainty stats')
  
  # get simranks
  SA_scores <-
  SA_scores %>% 
    .[, rank := frank(score, 
                      na.last='keep', 
                      ties.method = 'min'), by=rep]
    
  #also create bins using your chosen classification method
  if('Classification' %in% names(XX_p)) {
    SA_scores <- 
    XX_p$Classification %>% 
      as.data.table() %>% 
      .[, rep:=.I] %>% 
      setnames('.', 'classifier') %>% 
      merge(SA_scores, by='rep')

    SA_scores <-
      SA_scores %>% 
      .[, bin_size := floor(sum(!is.na(score))/10),by=rep] %>% 
      #bin the ranks and number them by rounding up
      .[classifier=='deciles', brank := (rank / bin_size) %>% ceiling,by=rep] %>% 
      .[classifier=='deciles' & brank>10, brank := 10] %>% 
      .[classifier=='equal_int', int_brank := cut(score, breaks = 10, labels = 1:10), by=rep] %>% 
      .[classifier=='equal_int', brank := int_brank %>% as.integer] %>% 
      .[classifier=='zscore', brank := (score - mean(score)) / sd(score), by=rep] %>% 
      .[classifier!='zscore', pct_80 := quantile(brank, .8), by=rep] %>% 
      .[classifier=='zscore', pct_80 := 0.842] %>% #80th percentile z score
      .[, impacted := 0] %>% 
      .[, impacted := (brank>=pct_80)] #test against the 80th pct cutoff 

  } else {
    #use standard method if not permuting the classifier
    SA_scores <-
      SA_scores %>% 
      .[, bin_size := floor(sum(!is.na(score))/10),by=rep] %>% 
      #bin the ranks and number them by rounding up
      .[, brank := (rank / bin_size) %>% ceiling,by=rep] %>% 
      .[brank>10, brank := 10] 
      .[, impacted := 0] %>% 
      .[brank>8, impacted := 1]
  }
  
  #also make a binned version to correspond with EHD
  # rank_cols <- names(SA_ranks)[-1]
  # SA_branks <- copy(SA_ranks) %>% 
  #   as.data.table %>% 
  #   .[, (rank_cols) := lapply(.SD, n_brank), .SDcols=rank_cols]
  
  # # get ranks, but just the ones from the SA/UA. If SA, only keep first 2N cols
  # # which correspond to random sampling.
  # SA_ranks_ <- SA_ranks[names(SA_ranks) %nin% c("uCode", "Nominal")]
  # SA_branks_ <- SA_branks[, c(names(SA_branks) %nin% c("uCode", "Nominal")), with=F]
  # if(SA_type == "SA"){
  #   SA_ranks_ <- SA_ranks_[, 1:(2*N)]
  #   SA_branks_ <- SA_branks_[, 1:(2*N)]
  # }

  #rank stats
  RankStats <- 
    #corresponds to random sampling (2N), .(Nominal=Nominal_rank,
    SA_scores[rep %in% c(1:(2*N)), 
              .(Nominal_score=Nominal,
                Nominal_rank=Nominal_rank,
                Nominal_bin=Nominal_brank,
                Q5=quantile(rank, probs = 0.05, na.rm = TRUE),
                Q25=quantile(rank, probs = 0.25, na.rm = TRUE),
                Mean=mean(rank, na.rm=TRUE),
                Median=median(rank, na.rm=TRUE),
                Q75=quantile(rank, probs = 0.75, na.rm = TRUE),
                Q95=quantile(rank, probs = 0.95, na.rm = TRUE)), by=uCode] %>% 
    unique
  
  #brank stats
  BRankStats <- 
    #corresponds to random sampling (2N), .(Nominal=Nominal_rank,
    SA_scores[rep %in% c(1:(2*N)), 
              .(Nominal_rank=Nominal_rank,
                Nominal_bin=Nominal_brank,
                Q5=quantile(brank, probs = 0.05, na.rm = TRUE),
                Q25=quantile(brank, probs = 0.25, na.rm = TRUE),
                Mean=mean(brank, na.rm=TRUE),
                Median=median(brank, na.rm=TRUE),
                Q75=quantile(brank, probs = 0.75, na.rm = TRUE),
                Q95=quantile(brank, probs = 0.95, na.rm = TRUE),
                Accuracy=sum(Nominal_impacted==impacted)/(.N)), by=uCode] %>% 
    unique(by='uCode')

  # Build list to output
  
  # First remove the measurement error DTs from memory because they take up too much space and are needless simulation
  #TODO if necessary we could try to just preserve the index to recreate 
  XX_p$Measurement_Error <- NULL
  
  SA_out <- list(
    # Scores = SA_scores,
    # Ranks = SA_ranks,
    # Branks = SA_branks,
    RankStats = RankStats,
    BRankStats = BRankStats,
    Para = XX_p
  )
  
  # get sensitivity indices if SA
  if(SA_type == "SA"){
    
    message('beginning sensitivity analysis')

    # An easy target is the mean absolute rank change
    # decide whether you want to use the raw rank or a binned ranking like EHD
    if(use_branks) y_AvDiffs <- SA_scores[, abs(Nominal_brank - brank) %>% mean, by=rep][, V1]
    else y_AvDiffs <- SA_scores[, abs(Nominal_rank - rank) %>% mean, by=rep][, V1]
    
    class_accuracy <- SA_scores[, sum(Nominal_impacted==impacted)/(.N), by=rep][, V1]

    #calculate the average differences
    # message('calculating differences in rank')
    # y_AvDiffs <- 
    # y_AvDiffs <- apply(target_obj[names(target_obj) %nin% c("uCode", "Nominal")], 2,
    #                    FUN = function(x) mean(abs(x-target_obj$Nominal), na.rm = TRUE) )
    # SA_out$diff <- y_AvDiffs

    message('calculating sensitivity stats')
    SA_out$Sensitivity <- 
      sobol_indices(Y=y_AvDiffs, params=names(SA_specs), order='second', N=N, boot=T, R=Nboot,
                    parallel = 'multicore', ncpus = ncores) %T>%
      print() %>% 
      .$results
    
    SA_out$Accuracy <- 
      sobol_indices(Y=class_accuracy, params=names(SA_specs), order='second', N=N, boot=T, R=Nboot,
                    parallel = 'multicore', ncpus = ncores) %T>%
      print() %>% 
      .$results
    
    #also add to a data.table for more functionality splitting out
    #scrape out all the parameter combinations for a more detailed SA
    extractParams <- function(x,i,diffs, the_list) {
      
      samp <- 
        the_list[[x]][[i]][[1]]

      if(samp %>% is.data.table) samp <- paste0('sim_', i)
      else if (samp %>% is.character) samp <- paste(samp, collapse=', ' )
      else samp <- unlist(samp)

      #TODO, also setup to calculate Si / STi at this level
      data.table(sim=i,
                 param=x, 
                 sample=ifelse(length(samp)==1, 
                               samp, 
                               paste0('sim_',i)),
                 average_diff=diffs[i]) %>% 
        return
      
      
    }

    sims <- expand.grid(x=names(XX_p),
                        y=1:nrow(XX))

    para_dt <-
      mapply(extractParams, 
             x=sims$x,
             i=sims$y,
             MoreArgs = list(diffs=y_AvDiffs, the_list=XX_p),
             SIMPLIFY=F) %>% 
      rbindlist
    
    SA_out$diffs_dt <- para_dt
    
  }

  # SA_out$Nominal <- data.frame(uCode = SA_scores$uCode,
  #                              Score = SA_scores$Nominal,
  #                              Rank = SA_ranks$Nominal,
  #                              BRank = SA_ranks$Nominal %>% n_brank)
  
  # timing
  tf <- proc.time()
  tdiff <- tf-t0
  telapse <- as.numeric(tdiff[3])
  taverage <- telapse/NT
  
  if(!quietly){
    message(paste0("Time elapsed = ", round(telapse,2), "s, average ", round(taverage,2), "s/rep."))
  }
  
  SA_out
}


# Regenerate an edited coin
#
# This is similar to [edit_coin()] but works with a list of parameters to change, rather than one, and also outputs
# a regenerated coin.
#
# @param l_para A list of parameter values to change. Should be of the format `list(para_name = new_value)`, where
# `new_value`
# @param addresses A list or character vector of addresses. `names(addresses)` must correspond to `names(l_para)`.
# @param coin A coin, to be edited.
#
# @return A regenerated coin
# @export
#
# @examples
# #
regen_edit <- function(l_para, addresses, coin){
  
  d <- length(l_para)
  p_names <- names(l_para)
  stopifnot(length(addresses) == d,
            is.coin(coin),
            setequal(names(addresses), p_names))
  
  # copy
  coin_i <- coin

  # modify parameters
  for(ii in 1:d){
    coin_i <- edit_coin(coin_i, address = addresses[names(addresses) == p_names[ii]],
                        new_value = l_para[[ii]])
  }

  # regenerate the results
  tryCatch(
    expr = Regen(coin_i, quietly = TRUE),
    error = function(e){
      message("Regen failed. Probably a conflict between methods.")
      message(l_para)
      message(addresses)
      stop()
      return(NA)
    }
  )
  
}


# Convert a numeric sample to parameter values
#
# Converts a numeric sample `x`, which should have values between 0 and 1, to a corresponding vector or list of
# parameter values, based on `distribution`.
#
# The `distribution` argument specifies how to map `x` to parameter values and can be used in two different ways,
# depending on `dist_type`. If `dist_type = "discrete"`, then `distribution` should be a vector or list of alternative
# parameter values (or objects). Each entry of `x` is mapped to an entry from `distribution` by treating `distribution`
# as a discrete uniform distribution over its entries.
#
# If `dist_type = "continuous"`, `distribution` is assumed to be a continuous uniform distribution, such that
# `distribution` is a 2-length numeric vector with the first value being the lower bound, and the second value the
# upper bound. For example, if `distribution = c(5, 10)`, then `x` will be mapped onto a continuous uniform distribution
# bounded by 5 and 10.
#
# @param distribution The distribution to sample using `x` - see details.
# @param dist_type Either `"discrete"` or `"continuous"` - see details.
# @param checks Logical: if `TRUE` runs some checks on inputs, set to `FALSE` to increase speed slightly.
# @param x A numeric vector with values between 0 and 1
#
# @return A vector or list of parameter values.
#
# @examples
# #
sample_2_para <- function(x, distribution, dist_type = "discrete", checks = TRUE){
  
  if(checks){
    stopifnot(is.numeric(x),
              any(x >= 0),
              any(x <= 1),
              is.character(dist_type),
              length(dist_type) == 1,
              dist_type %in% c("discrete", "continuous"))
  }
  
  # specs can be a set of discrete alternatives, or else a uniform distribution
  if(dist_type == "discrete"){
    
    # the number of discrete alternatives
    n_alt <- length(distribution)
    # convert x to indexes of the discrete parameters
    i_para <- cut(x, n_alt, 1:n_alt)
    # now get the output vector/list
    l_out <- distribution[i_para]
    
  } else {
    
    # here we assume a uniform distribution
    if(checks){
      stopifnot(is.numeric(distribution),
                length(distribution) == 2,
                distribution[2] > distribution[1])
    }
    
    # we simply scale x up to the interval covered by the distribution
    l_out <- x*(distribution[2] - distribution[1]) + distribution[1]
    
  }
  
  # output
  l_out
  
  
}


# Edit objects inside a coin
#
# Changes the object found at `address` to `new_value`.
#
# @param coin A coin
# @param address A string specifying the location in the coin of the object to edit. This should begin with `"$"`, omitting the coin itself
# in the address. E.g. if you target `coin$x$y$z` enter `"$x$y$z"`.
# @param new_value The new value to assign at `address`.
# @param checks Logical: if `TRUE`, runs some basic checks, otherwise omitted if `FALSE`. Setting `FALSE` may speed
# things up a bit in sensitivity analysis, for example.
#
# @return An updated coin
# @export
#
# @examples
# #
edit_coin <- function(coin, address, new_value, checks = FALSE){
  
  # checks
  if(checks){
    check_coin_input(coin)
    stopifnot(is.character(address),
              length(address) == 1,
              substr(address,1,1) == "$")
  }
  
  # this is the call to evaluate, as a string
  expr_str <- paste0("coin", address, " <- new_value")
  # evaluate the call
  eval(str2lang(expr_str))
  
  # output
  coin
  
}

# Check address in coin
check_address <- function(address, coin){
  
  # checks
  stopifnot(is.character(address),
            length(address) == 1)
  
  # check address begins with $
  if(substr(address,1,1) != "$"){
    stop("Address must begin with '$'! Your address: ", address, call. = FALSE)
  }
  
  # this is the call to evaluate, as a string
  expr_str <- paste0("coin", address)
  # evaluate the call
  address_value <- eval(str2lang(expr_str))
  
  if(is.null(address_value)){
    xx <- readline(paste0("Address ", address, " is not currently present in the coin or else is NULL. Continue anyway (y/n)?  "))
    
    if(xx %nin% c("y", "n")){
      stop("You didn't input y or n. I'm taking that as a no.", call. = FALSE)
    }
    
    if(xx == "n"){
      stop("Exiting sensitivity analysis. Please check the address: ", address, call. = FALSE)
    }
  }
  
}


#' Estimate sensitivity indices
#'
#' Post process a sample to obtain sensitivity indices. This function takes a univariate output
#' which is generated as a result of running a Monte Carlo sample from [SA_sample()] through a system.
#' Then it estimates sensitivity indices using this sample.
#'
#' This function is built to be used inside [get_sensitivity()].
#'
#' @param yy A vector of model output values, as a result of a \eqn{N(d+2)} Monte Carlo design.
#' @param N The number of sample points per dimension.
#' @param d The dimensionality of the sample
#' @param Nboot Number of bootstrap draws for estimates of confidence intervals on sensitivity indices.
#' If this is not specified, bootstrapping is not applied.
#'
#' @importFrom stats var
#'
#' @examples
#' # This is a generic example rather than applied to a COIN (for reasons of speed)
#'
#' # A simple test function
#' testfunc <- function(x){
#' x[1] + 2*x[2] + 3*x[3]
#' }
#'
#' # First, generate a sample
#' X <- SA_sample(500, 3)
#'
#' # Run sample through test function to get corresponding output for each row
#' y <- apply(X, 1, testfunc)
#'
#' # Estimate sensitivity indices using sample
#' SAinds <- SA_estimate(y, N = 500, d = 3, Nboot = 1000)
#' SAinds$SensInd
#' # Notice that total order indices have narrower confidence intervals than first order.
#'
#' @seealso
#' * [get_sensitivity()] Perform global sensitivity or uncertainty analysis on a COIN
#' * [SA_sample()] Input design for estimating sensitivity indices
#'
#' @return A list with the output variance, plus a data frame of first order and total order sensitivity indices for
#' each variable, as well as bootstrapped confidence intervals if `!is.null(Nboot)`.
#'
#' @export
SA_estimate <- function(yy, N, d, Nboot = NULL){
  
  # checks
  stopifnot(is.numeric(yy))
  if(length(yy) != N*(d+2)){
    stop("The length of 'yy' does not correspond to the values of 'N' and 'd'. The vector 'yy' should be of length N(d+2).")
  }
  
  # put into matrix format: just the ABis
  yyABi <- matrix(yy[(2*N +1):length(yy)], nrow = N)
  # get yA and yB
  yA <- yy[1:N]
  yB <- yy[(N+1) : (2*N)]
  
  # calculate variance
  varY <- stats::var(c(yA,yB))
  
  # calculate Si
  Si <- apply(yyABi, 2, function(x){
    mean(yB*(x - yA))/varY
  })
  
  # calculate ST
  STi <- apply(yyABi, 2, function(x){
    sum((yA - x)^2)/(2*N*varY)
  })
  
  # make a df
  SensInd <- data.frame(Variable = paste0("V", 1:d),
                        Si = Si,
                        STi = STi)

  ## BOOTSTRAP ## -----
  
  if (!is.null(Nboot)){
    
    # Get the "elements" to sample from
    STdiffs <- apply(yyABi, 2, function(x){
      yA - x
    })
    Sidiffs <- apply(yyABi, 2, function(x){
      yB*(x - yA)
    })
    
    # prep matrices for bootstrap samples
    Si_boot <- matrix(NA, d, Nboot)
    STi_boot <- Si_boot
    
    # do the bootstrapping bit
    for (iboot in 1:Nboot){
      
      # calculate Si
      Si_boot[,iboot] <- apply(Sidiffs, 2, function(x){
        mean(sample(x, replace = TRUE))/varY
      })
      
      # calculate ST
      STi_boot[,iboot] <- apply(STdiffs, 2, function(x){
        sum(sample(x, replace = TRUE)^2)/(2*N*varY)
      })
      
    }
    
    # get quantiles of sensitivity indices to add to
    SensInd$Si_q5 <- apply(Si_boot, MARGIN = 1,
                           function(xx) stats::quantile(xx, probs = 0.05, na.rm = TRUE))
    SensInd$Si_q95 <- apply(Si_boot, MARGIN = 1,
                            function(xx) stats::quantile(xx, probs = 0.95, na.rm = TRUE))
    SensInd$STi_q5 <- apply(STi_boot, MARGIN = 1,
                            function(xx) stats::quantile(xx, probs = 0.05, na.rm = TRUE))
    SensInd$STi_q95 <- apply(STi_boot, MARGIN = 1,
                             function(xx) stats::quantile(xx, probs = 0.95, na.rm = TRUE))
    
  }
  
  # return outputs
  list(Variance = varY,
       SensInd = SensInd)
  
}


#' Generate sample for sensitivity analysis
#'
#' Generates an input sample for a Monte Carlo estimation of global sensitivity indices. Used in
#' the [get_sensitivity()] function. The total sample size will be \eqn{N(d+2)}.
#'
#' This function generates a Monte Carlo sample as described e.g. in the [Global Sensitivity Analysis: The Primer book](https://onlinelibrary.wiley.com/doi/book/10.1002/9780470725184).
#'
#' @param N The number of sample points per dimension.
#' @param d The dimensionality of the sample
#'
#' @importFrom stats runif
#'
#' @examples
#' # sensitivity analysis sample for 3 dimensions with 100 points per dimension
#' X <- SA_sample(100, 3)
#'
#' @return A matrix with \eqn{N(d+2)} rows and `d` columns.
#'
#' @seealso
#' * [get_sensitivity()] Perform global sensitivity or uncertainty analysis on a COIN.
#' * [SA_estimate()] Estimate sensitivity indices from system output, as a result of input design from SA_sample().
#'
#' @export
SA_sample <- function(N, d){
  
  # a random (uniform) sample
  Xbase <- matrix(stats::runif(d*N*2), nrow = N, ncol = d*2)
  # get first half
  XA <- Xbase[, 1:d]
  # get second half
  XB <- Xbase[, (d+1):(2*d)]
  # make big matrix (copy matrix d times on the bottom)
  XX <- matrix(rep(t(XA), d), ncol = ncol(XA), byrow = TRUE )
  
  # now substitute in columns from B into A
  for (ii in 1:d){
    XX[(1 + (ii-1)*N):(ii*N), ii] <- XB[, ii]
  }
  
  # add original matrices on the beginning
  XX <-  rbind(XA, XB, XX)
  
  XX
}


#' Plot ranks from an uncertainty/sensitivity analysis
#'
#' Plots the ranks resulting from an uncertainty and sensitivity analysis, in particular plots
#' the median, and 5th/95th percentiles of ranks.
#'
#' To use this function you first need to run [get_sensitivity()]. Then enter the resulting list as the
#' `SAresults` argument here.
#'
#' See `vignette("sensitivity")`.
#'
#' This function replaces the now-defunct `plotSARanks()` from COINr < v1.0.
#'
#' @param SAresults A list of sensitivity/uncertainty analysis results from [get_sensitivity()].
#' @param plot_units A character vector of units to plot. Defaults to all units. You can also set
#' to `"top10"` to only plot top 10 units, and `"bottom10"` for bottom ten.
#' @param order_by If set to `"nominal"`, orders the rank plot by nominal ranks
#' (i.e. the original ranks prior to the sensitivity analysis). Otherwise if `"median"`, orders by
#' median ranks.
#' @param dot_colour Colour of dots representing median ranks.
#' @param line_colour Colour of lines connecting 5th and 95th percentiles.
#'
#' @importFrom ggplot2 geom_line geom_point scale_shape_manual scale_size_manual labs guides
#' @importFrom ggplot2 scale_color_manual theme_classic theme scale_y_discrete scale_x_reverse
#' @importFrom ggplot2 coord_flip element_text
#'
#' @examples
#' # for examples, see `vignette("sensitivity")`
#' # (this is because package examples are run automatically and sensitivity analysis
#' # can take a few minutes to run at realistic settings)
#'
#' @seealso
#' * [get_sensitivity()] Perform global sensitivity or uncertainty analysis on a coin
#' * [plot_sensitivity()] Plot sensitivity indices following a sensitivity analysis.
#'
#' @return A plot of rank confidence intervals, generated by 'ggplot2'.
#'
#' @export
plot_uncertainty <- function(SAresults, 
                             plot_var='RankStats',
                             plot_units = NULL, order_by = "Nominal_Rank"){
  
  rnks <- SAresults[[plot_var]] %>% 
    as.data.table

  if(!is.null(plot_units)){
    if(length(plot_units) == 1){
      
      if (plot_units == "top10"){
        unit_include <- SAresults$Nominal$uCode[SAresults$Nominal$Rank <= 10]
      } else if (plot_units == "bottom10"){
        unit_include <- SAresults$Nominal$uCode[
          SAresults$Nominal$Rank >= (max(SAresults$Nominal$Rank, na.rm = TRUE)-10)]
      } else {
        stop("plot_units not recognised: should be either a character vector of unit codes or else
      \"top10\" or \"bottom10\" ")
      }
    } else {
      # vector, so this should be a vector of unit codes
      unit_include <- plot_units
      if(any(unit_include %nin% rnks$uCode)){
        stop("One or more units in 'plot_units' not found in SA results.")
      }
    }
    rnks <- rnks[rnks$uCode %in% unit_include,]
    
  }

  #calculate difference 
  rnks[, diff := Median - Nominal_rank]

  # generate plot
  #TODO find out why .data won't work with function args as strings..its supposed to maybe version issues
    ggplot(rnks, aes(x = fct_reorder(.data[['uCode']], .data[['Nominal_rank']]), y = .data[['Median']])) +
      geom_errorbar(aes(ymin=.data[['Q5']], ymax=.data[['Q95']]), width=.05,
                    position=position_dodge(0.05)) +
      geom_point(aes(color = diff)) +
      theme(legend.position="top") +
      scale_color_gradient2('Deviation \n from EHD Rank', na.value = "grey75") +
      scale_x_discrete('Census Tract') +
      scale_y_continuous('Median Ranking') +
      theme_minimal() +
      theme(axis.text.x = element_blank()) 

}


#' Plot sensitivity indices
#'
#' Plots sensitivity indices as bar or pie charts.
#'
#' To use this function you first need to run [get_sensitivity()]. Then enter the resulting list as the
#' `SAresults` argument here.
#'
#' See `vignette("sensitivity")`.
#'
#' This function replaces the now-defunct `plotSA()` from COINr < v1.0.
#'
#' @param SAresults A list of sensitivity/uncertainty analysis results from [plot_sensitivity()].
#' @param ptype Type of plot to generate - either `"bar"`, `"pie"` or `"box"`.
#'
#' @importFrom ggplot2 ggplot aes geom_bar labs theme_minimal geom_errorbar coord_polar theme_void
#' @importFrom ggplot2 facet_wrap
#' @importFrom rlang .data
#'
#' @examples
#' # for examples, see `vignette("sensitivity")`
#' # (this is because package examples are run automatically and sensitivity analysis
#' # can take a few minutes to run at realistic settings)
#'
#' @return A plot of sensitivity indices generated by ggplot2.
#'
#' @seealso
#' * [get_sensitivity()] Perform global sensitivity or uncertainty analysis on a COIN
#' * [plot_uncertainty()] Plot confidence intervals on ranks following a sensitivity analysis
#'
#' @export
plot_sensitivity <- function(SAresults, ptype = "bar"){
  
  stopifnot(is.list(SAresults))
  # prep data first
  Sdf <- SAresults$Sensitivity
  if(is.null(Sdf)){
    stop("Sensitivity indices not found. Did you run get_sensitivity with SA_type = 'SA'?")
  }
  numcols <- Sdf[names(Sdf) != "Variable"]
  
  # set any negative values to zero. By definition they can't be negative.
  numcols[numcols < 0] <- 0
  Sdf[names(Sdf) != "Variable"] <- numcols
  
  if(ptype == "bar"){
    
    # the full bar is STi. It is divided into Si and the remainder, so we need STi - Si
    Sdf$Interactions = Sdf$STi - Sdf$Si
    Sdf$Interactions[Sdf$Interactions < 0] <- 0
    
    # rename col to improve plot
    colnames(Sdf)[colnames(Sdf) == "Si"] <- "MainEffect"
    # now pivot to get in format for ggplot
    bardf <- lengthen(Sdf, cols = c("MainEffect", "Interactions"))
    # bardf <- tidyr::pivot_longer(Sdf,
    #                              cols = c("MainEffect", "Interactions"))
    
    # make stacked bar plot
    plt <- ggplot2::ggplot(bardf, ggplot2::aes(fill=.data$name, y=.data$Value, x=.data$Variable)) +
      ggplot2::geom_bar(position="stack", stat="identity") +
      ggplot2::labs(
        x = NULL,
        y = NULL,
        fill = NULL) +
      ggplot2::theme_minimal()
    
    
  } else if (ptype == "pie"){
    
    # we are plotting first order sensitivity indices. So, also estimate interactions.
    Sis <- Sdf[c("Variable", "Si")]
    intsum <- max(c(1 - sum(Sis$Si, na.rm = TRUE), 0))
    Sis <- rbind(Sis, data.frame(Variable = "Interactions", Si = intsum))
    
    # Basic piechart
    plt <- ggplot2::ggplot(Sis, ggplot2::aes(x = "", y = .data$Si, fill = .data$Variable)) +
      ggplot2::geom_bar(stat="identity", width=1, color="white") +
      ggplot2::coord_polar("y", start=0) +
      ggplot2::theme_void() # remove background, grid, numeric labels
    
  } else if (ptype == "box"){
    
    if(any(c("Si_q5", "Si_q95", "STi_q5", "STi_q95") %nin% names(Sdf))){
      stop("Quantiles not found for sensitivity indices (required for box plot). Did you forget to set Nboot when running get_sensitivity()?")
    }
    Sdf <- lengthen(Sdf, cols = c("Si", "STi"))
    #Sdf1 <- tidyr::pivot_longer(Sdf, cols = c("Si", "STi"))
    Sdf$q5 <- ifelse(Sdf$name == "STi", Sdf$STi_q5, Sdf$Si_q5)
    Sdf$q95 <- ifelse(Sdf$name == "STi", Sdf$STi_q95, Sdf$Si_q95)
    Sdf$q5[Sdf$q5 > 1] <- 1
    Sdf$q95[Sdf$q95 > 1] <- 1
    Sdf$Value[Sdf$Value > 1] <- 1
    
    plt <- ggplot2::ggplot(Sdf, ggplot2::aes(x = .data$Variable, y = .data$Value, ymax = .data$q95, ymin = .data$q5)) +
      ggplot2::geom_point(size = 1.5) +
      ggplot2::geom_errorbar(width = 0.2) +
      ggplot2::theme_bw() +
      facet_wrap(~name) +
      ggplot2::labs(
        x = NULL,
        y = NULL)
    
  }
  
  plt +
    ggplot2::theme(text=element_text(family="sans"))
  
}


#' Noisy replications of weights
#'
#' Given a data frame of weights, this function returns multiple replicates of the weights, with added
#' noise. This is intended for use in uncertainty and sensitivity analysis.
#'
#' Weights are expected to be in a data frame format with columns `Level`, `iCode` and `Weight`, as
#' used in `iMeta`. Note that no `NA`s are allowed anywhere in the data frame.
#'
#' Noise is added using the `noise_specs` argument, which is specified by a data frame with columns
#' `Level` and `NoiseFactor`. The aggregation level refers to number of the aggregation level to target
#' while the `NoiseFactor` refers to the size of the perturbation. If e.g. a row is `Level = 1` and
#' `NoiseFactor = 0.2`, this will allow the weights in aggregation level 1 to deviate by +/- 20% of their
#' nominal values (the values in `w`).
#'
#' This function replaces the now-defunct `noisyWeights()` from COINr < v1.0.
#'
#' @param w A data frame of weights, in the format found in `.$Meta$Weights`.
#' @param noise_specs a data frame with columns:
#'  * `Level`: The aggregation level to apply noise to
#'  * `NoiseFactor`: The size of the perturbation: setting e.g. 0.2 perturbs by +/- 20% of nominal values.
#' @param Nrep The number of weight replications to generate.
#'
#' @examples
#' # build example coin
#' coin <- build_example_coin(up_to = "new_coin", quietly = TRUE)
#'
#' # get nominal weights
#' w_nom <- coin$Meta$Weights$Original
#'
#' # build data frame specifying the levels to apply the noise at
#' # here we vary at levels 2 and 3
#' noise_specs = data.frame(Level = c(2,3),
#'                          NoiseFactor = c(0.25, 0.25))
#'
#' # get 100 replications
#' noisy_wts <- get_noisy_weights(w = w_nom, noise_specs = noise_specs, Nrep = 100)
#'
#' # examine one of the noisy weight sets, last few rows
#' tail(noisy_wts[[1]])
#'
#' @return A list of `Nrep` sets of weights (data frames).
#'
#' @seealso
#' * [get_sensitivity()] Perform global sensitivity or uncertainty analysis on a COIN
#'
#' @export
get_noisy_weights <- function(w, noise_specs, Nrep){
  
  # CHECKS
  
  stopifnot(is.data.frame(noise_specs))
  
  if(any(is.na(w))){
    stop("NAs found in w: NAs are not allowed.")
  }
  
  if(any(c("iCode", "Level", "Weight") %nin% names(w))){
    stop("One or more required columns (iCode, Level, Weight) not found in w.")
  }
  
  if (length(unique(noise_specs$Level)) < nrow(noise_specs)){
    stop("Looks like you have duplicate Level values in the noise_specs df?")
  }
  
  # make list for weights
  wlist <- vector(mode = "list", length = Nrep)
  
  for (irep in 1:Nrep){
    
    # make fresh copy of weights
    wrep <- w
    
    for (ii in noise_specs$Level){
      # weights for this level
      wts <- wrep$Weight[w$Level == ii]
      # vector of noise: random number in [0,1] times 2, -1. This interprets NoiseFactor as
      # a +/-% deviation.
      wnoise <- (stats::runif(length(wts))*2 - 1)*noise_specs$NoiseFactor[noise_specs$Level == ii]*wts
      # add noise to weights and store
      wts <- wts + wnoise
      wrep$Weight[w$Level == ii] <- wts
    }
    wlist[[irep]] <- wrep
  }
  
  return(wlist)
}

# Data frame or matrix to long form
#
# This is a substitute function for tidyr's 'pivot_longer' to avoid dependencies, and behaves in more or
# less the same way.
#
# If `cols` is not specified, assumes a square correlation matrix to convert to long form. If `cols` is
# specified, this behaves like pivot_longer's "cols" argument.
#
# @param X A data frame or square correlation matrix
# @param cols Columns to pivot into longer format.
#
# @importFrom utils stack
#
# @return A long format data frame
lengthen <- function(X, cols = NULL){
  
  # make df
  X <- as.data.frame(X)
  
  if(!is.null(cols)){
    
    stopifnot(all(cols %in% names(X)))
    X_ <- X[cols]
    X <- X[names(X) %nin% cols]
    X$V_to_pivot <- rownames(X)
    
  } else {
    X_ <- X
  }
  
  # stack and add names
  X1 <- cbind(utils::stack(X_), rownames(X_))
  names(X1) <- c("Value", "V2", "V1")
  X1$V2 <- as.character(X1$V2)
  X1 <- rev(X1)
  
  if(!is.null(cols)){
    X1 <- merge(X, X1, by.x = "V_to_pivot", by.y = "V1", all = TRUE)
    X1 <- X1[names(X1) != "V_to_pivot"]
    names(X1)[names(X1) == "V2"] <- "name"
  }
  
  X1
  
}

#TODO i don't know why this function by default inverts the rank order
#made my own version to keep the ranks ascending like EHD
rank_df <-
function (df, use_group = NULL) 
{
  dfo <- df
  if (is.null(use_group)) {
    df <- data.frame(lapply(df, function(y) if (is.numeric(y)) 
      rank(y, na.last = "keep", ties.method = "min")
      else y))
  }
  else {
    stopifnot(use_group %in% colnames(df))
    grps <- unique(unlist(df[[use_group]]))
    dfold <- df
    for (grp in grps) {
      grprows <- df[[use_group]] == grp
      grprows[is.na(grprows)] <- FALSE
      df[grprows, ] <- data.frame(lapply(dfold[grprows, 
      ], function(y) if (is.numeric(y)) 
        rank(y, na.last = "keep", ties.method = "min")
      else y))
    }
    if (any(is.na(df[[use_group]]))) {
      df[is.na(df[[use_group]]), ] <- data.frame(lapply(df[is.na(df[[use_group]]), 
      ], function(y) if (is.numeric(y)) 
        NA
      else y))
    }
  }
  rownames(df) <- NULL
  names(df) <- names(dfo)
  df
}

#' Aggregate indicators
#'
#' Aggregates a named data set specified by `dset` using aggregation function `f_ag`, weights `w`, and optional
#' function parameters `f_ag_para`. Note that COINr has a number of aggregation functions built in,
#' all of which are of the form `a_*()`, e.g. [a_amean()], [a_gmean()] and friends.
#'
#' Aggregation is performed row-wise using the function `f_ag`, such that for each row `x_row`, the output is
#' `f_ag(x_row, f_ag_para)`, and for the whole data frame, it outputs a numeric vector. The data frame `x` must
#' only contain numeric columns.
#'
#' The function `f_ag` must be supplied as a string, e.g. `"a_amean"`, and it must take as a minimum an input
#' `x` which is either a numeric vector (if `by_df = FALSE`), or a data frame (if `by_df = TRUE`). In the former
#' case `f_ag` should return a single numeric value (i.e. the result of aggregating `x`), or in the latter case
#' a numeric vector (the result of aggregating the whole data frame in one go).
#'
#' `f_ag` can optionally have other parameters, e.g. weights, specified as a list in `f_ag_para`.
#'
#' Note that COINr has a number of aggregation functions built in,
#' all of which are of the form `a_*()`, e.g. [a_amean()], [a_gmean()] and friends. To see a list browse COINr functions alphabetically or
#' type `a_` in the R Studio console and press the tab key (after loading COINr).
#'
#' Optionally, a data availability threshold can be assigned below which the aggregated value will return
#' `NA` (see `dat_thresh` argument). If `by_df = TRUE`, this will however be ignored because aggregation is not
#' done on individual rows. Note that more complex constraints could be built into `f_ag` if needed.
#'
#' @param x A coin class object.
#' @param dset The name of the data set to apply the function to, which should be accessible in `.$Data`.
#' @param f_ag The name of an aggregation function, a string. This can either be a single string naming
#' a function to use for all aggregation levels, or else a character vector of function names of length `n-1`, where `n` is
#' the number of levels in the index structure. In this latter case, a different aggregation function may be used for each level
#' in the index: the first in the vector will be used to aggregate from Level 1 to Level 2, the second from Level 2 to Level 3, and
#' so on.
#' @param w An optional data frame of weights. If `f_ag` does not require accept weights, set to `"none"`.
#' @param f_ag_para Optional parameters to pass to `f_ag`, other than `x` and `w`. As with `f_ag`, this can specified to have different
#' parameters for each aggregation level by specifying as a nested list of length `n-1`.
#' @param dat_thresh An optional data availability threshold, specified as a number between 0 and 1. If a row
#' within an aggregation group has data availability lower than this threshold, the aggregated value for that row will be
#' `NA`. Data availability, for a row `x_row` is defined as `sum(!is.na(x_row))/length(x_row)`, i.e. the
#' fraction of non-`NA` values.
#' @param by_df Controls whether to send a numeric vector to `f_ag` (if `FALSE`, default) or a data frame (if `TRUE`) - see
#' details.
#' @param out2 Either `"coin"` (default) to return updated coin or `"df"` to output the aggregated data set.
#' @param write_to If specified, writes the aggregated data to `.$Data[[write_to]]`. Default `write_to = "Aggregated"`.
#' @param ... arguments passed to or from other methods.
#'
#' @examples
#' # build example up to normalised data set
#' coin <- build_example_coin(up_to = "Normalise")
#'
#' # aggregate normalised data set
#' coin <- Aggregate(coin, dset = "Normalised")
#'
#' @return An updated coin with aggregated data set added at `.$Data[[write_to]]` if `out2 = "coin"`,
#' else if `out2 = "df"` outputs the aggregated data set as a data frame.
#'
#' @export
Aggregate.coin <- function(x, dset, f_ag = NULL, w = NULL, f_ag_para = NULL, dat_thresh = NULL,
                           flatten_hierarchy = F,
                           by_df = FALSE, out2 = "coin", write_to = NULL, ...){
  
  # Write to Log ------------------------------------------------------------
  
  coin <- write_log(x, dont_write = "x")
  
  # CHECK AND SET f_ag ------------------------------------------------------
  
  if(flatten_hierarchy) nlev <- 2
  else nlev <- max(coin$Meta$Ind$Level, na.rm = TRUE)
  
  # default and check
  if(is.null(f_ag)){
    f_ag <- "a_amean"
    f_ag_para <- NULL
  } else {
    if(!is.character(f_ag)){
      stop("f_ag must be specified as a character string or vector (function name(s) in inverted commas).")
    }
  }
  stopifnot(length(f_ag) > 0)
  
  # if same for all levels, repeat
  if(length(f_ag) == 1){
    f_ags <- rep(f_ag, nlev - 1)
  } else {
    if(length(f_ag) != (nlev - 1) & !flatten_hierarchy){
      stop("f_ag must have either length 1 (same function for all levels) or length equal to number of levels - in your case: ", nlev)
      
      
    } else if (flatten_hierarchy) { 
      f_ags <- f_ag[[1]]; message('flattening f_ag hierarchy to first option provided')
    } else f_ags <- f_ag
  }
  
  # CHECK AND SET w ---------------------------------------------------------
  
  # if weights is supplied we have to see what kind of thing it is
  # NULL indicates that we should use metadata weights
  if(!is.null(w)){
    
    if(is.data.frame(w)){
      
      stopifnot(exists("iCode", w),
                exists("Weight", w))
      w1 <- w
      
    } else if(is.character(w)){
      
      if(length(w) != 1){
        stop("w must be either a string indicating a name of a weight set, or a data frame of weights, or 'none', or NULL (to use weights from metadata).")
      }
      
      if(w != "none"){
        
        # we look for a named weight set
        w1 <- coin$Meta$Weights[[w]]
        if(is.null(w1)){
          stop("Weight set with name '", w, "' not found in .$Meta$Weights.")
        }
        stopifnot(is.data.frame(w1),
                  exists("iCode", w1),
                  exists("Weight", w1))
        
      } else {
        # convert w1 to NULL
        w1 <- NULL
      }
    } else {
      stop("w must be either a string indicating a name of a weight set, or a data frame of weights, or 'none', or NULL (to use weights from metadata).")
    }
    
  } else{
    # if w was NULL, get from metadata
    w1 <- coin$Meta$Ind[c("iCode", "Weight")]
  }
  
  # from this point, w1 is either a data frame of weights, or NULL (don't pass weights to f_ag)
  
  # CHECK AND SET f_ag_para -------------------------------------------------
  
  # if f_ag_para is NULL, repeat for all levs
  if(!is.null(f_ag_para)){
    if(!is.list(f_ag_para)){
      stop("f_ag_para must be specified as a list or list of lists")
    }
    stopifnot(length(f_ag_para) > 0)
    
    # if same for all levels, repeat
    if(length(f_ag_para) == 1){
      f_ag_paras <- rep(f_ag_para, nlev - 1)
    } else {
      if(length(f_ag_para) != (nlev - 1)){
        stop("f_ag_para must have either length 1 (same parameters for all levels) or length equal to number of levels - in your case: ", nlev)
      }
      f_ag_paras <- f_ag_para
    }
  } else {
    f_ag_paras <- rep(list(NULL), 4)
  }
  
  # Other Prep --------------------------------------------------------------------
  
  if(is.null(dat_thresh)){
    dat_threshs <- rep(list(NULL), 4)
  } else {
    if(!is.numeric(dat_thresh)){
      stop("dat_thresh must be a numeric value or vector of length (number of levels - 1) - in your case: ", nlev)
    }
    if(any((dat_thresh < 0) | (dat_thresh > 1))){
      stop("dat_thresh must only contain numeric values between 0 and 1.")
    }
    if(length(dat_thresh) == 1){
      dat_threshs <- rep(dat_thresh, nlev - 1)
    } else {
      if(length(dat_thresh) != (nlev - 1)){
        stop("dat_thresh must have either length 1 (same for all levels) or length equal to number of levels - in your case: ", nlev)
      }
      dat_threshs <- dat_thresh
    }
  }
  
  # Aggregate ---------------------------------------------------------------
  # Here we apply the aggregation by level
  
  # get data (also performing checks)
  indat <- get_dset(coin, dset)
  
  #if flattening hierarchy, only keep first and last level
  if(flatten_hierarchy) { 
    message('flattening hierarchy in the meta obj')
    
    #TODO maybe can use a central fx for this like edit_coin or regen
    #edit the meta object itself
    imeta <- coin$Meta$Ind[!is.na(coin$Meta$Ind$Level), ] %>% 
      copy %>% 
      .[Level %in% range(Level)] %>% 
      .[Level==max(Level), Level := 2] %>% 
      .[Level==min(Level), Parent := .[Level==max(Level), iCode]]

  } else imeta <- coin$Meta$Ind[!is.na(coin$Meta$Ind$Level), ] 
  
  # get metadata
  

  # Function that aggregates from Level = lev to the next level up
  # calls the function specified by f_ag.
  aggregate_level <- function(lev){
    
    # filter metadata to level
    imeta_l <- imeta[imeta$Level == (lev-1), ]

    if(is.null(w1)){
      aggs <- tapply(imeta_l$iCode, imeta_l$Parent, function(codes){
        # call func
        do.call("Aggregate",
                list(x = indat_ag[codes],
                     f_ag = f_ags[lev-1],
                     f_ag_para = f_ag_paras[[lev-1]],
                     dat_thresh = dat_threshs[lev-1],
                     by_df = by_df))
      })
    } else {
      aggs <- tapply(imeta_l$iCode, imeta_l$Parent, function(codes){
        # get weights
        wts <- w1$Weight[match(codes, w1$iCode)]
        
        # call func
        do.call("Aggregate",
                list(x = indat_ag[codes],
                     f_ag = f_ags[lev-1],
                     f_ag_para = c(list(w = wts), f_ag_paras[[lev-1]]),
                     dat_thresh = dat_threshs[[lev-1]],
                     by_df = by_df))
      })
    }
    
    if(is.data.frame(aggs)) { 
      return(aggs)
    } else if (is.list(aggs)){
      # aggs comes out as a list of vectors, have to make to df
      as.data.frame(do.call(cbind, aggs))
    } else if (is.numeric(aggs)){
      # in this case there is only one row, and it comes out as an array which needs to be converted
      as.data.frame(t(aggs))
    }  
    
  }

  # run the above function for each level
  indat_ag <- indat #TODO not sure why we need to define this before fx, seems weird namespace
  # indat_ag <- sapply(2:nlev, aggregate_level) %>% 
  #   .[!unlist(lapply(., is.null))] %>% #remove any NULLs 
  #   do.call(cbind, .) %>% #bind the results of the list
  #   cbind(indat_ag, .) #bind the results to our input data
  
  if('pca' %in% f_ag){
    
    #only use PCA to collapse non hierarchically
    agg <- get_PCA(coin, dset='Normalised', by_groups=F, out2='preds', imputed=T)
    indat_ag <- cbind(indat_ag, agg)

  } else {
  
    #TODO there is probably a better way to do this but each level needs the previous level data to agg further
    #lapply cant modify outside of scope so for loop works better
    for(lev in 2:nlev) {
      
      agg <- aggregate_level(lev)
      if(agg %>% is.data.frame) indat_ag <- cbind(indat_ag, agg)

    }
  }
  
  # Output ------------------------------------------------------------------
  
  # output list
  if(out2 == "df"){
    indat_ag
  } else {
    if(is.null(write_to)){
      write_to <- "Aggregated"
    }
    write_dset(coin, indat_ag, dset = write_to)
  }
  
  
}


# WRITING TO COINS

# Write function arguments to log
#
# This used inside `build_*` functions. It takes the coin object as an input, then writes the arguments of the
# current function into the `.$Log` list of the coin object. This is then used as a record of the operations used
# to build the coin, and can be edited.
#
# @param coin A coin class object
# @param dont_write Any variables not to write to the coin
# @param write2log If `FALSE`, just passes the coin back without writing anything.
#
# @examples
# #
#
# @return Updated coin object with function arguments written to `.$Log`
write_log <- function(coin, dont_write = NULL, write2log = TRUE){
  
  if(!write2log){
    return(coin)
  }
  
  # get calling function name and its arguments
  func_args <- as.list(sys.frame(-1))
  func_name <- deparse(as.list(sys.call(-1))[[1]])
  
  # check whether call is of type COINr::func_name, or just func_name
  if(grepl("::", func_name)){
    # split at ::, this returns a list (to deal with character vectors) - we only have one string so take first
    xx <- strsplit(func_name, "::")[[1]]
    if(xx[1] != "COINr"){
      stop("Attempt to write log from non-COINr function!")
    } else {
      # take function name excluding COINr:: bit
      func_name <- xx[2]
    }
  }
  
  # tweak list first (exclude args we don't want)
  dont_write <- c(dont_write, "coin", "*tmp*")
  func_args <- func_args[!(names(func_args) %in% dont_write)]
  
  # check that we are getting function arguments and nothing else
  if(!all(names(func_args) %in% names(formals(func_name)))){
    stop(paste0("Mismatch between function arguments of ", func_name, " and attempt to write to .$Log."))
  }
  
  # remove method .coin or similar from func_name
  func_name2 <- unlist(strsplit(func_name, "\\.")[[1]])[1]
  
  # make sure this is a builder function calling
  builders <- c("Aggregate", "Denominate", "Impute", "new_coin", "Normalise", "qNormalise", "qTreat", "Screen", "Treat")
  if(func_name2 %nin% builders){
    stop("The calling function ", func_name2, " is not one of the functions allowed to write to log. Authorised functions are ", builders)
  }
  
  # write to coin
  coin$Log[[func_name2]] <- func_args
  coin
  
}

# Direct function outputs
#
# Shortcut to be used at end of functions, either attach to coin or output as list.
#
# @param coin A coin class object
# @param l The list to direct
# @param out2 Whether list or attach to coin
# @param lev1 Address at first lev of coin
# @param lev2 Address at second lev of coin
# @param lev3 Address at third lev of coin
#
# @examples
# #
#
# @return Either a list or an updated coin
write2coin <- function(coin, l, out2, lev1, lev2, lev3 = NULL){
  
  check_coin_input(coin)
  
  if(out2 %nin% c("list", "coin")){
    stop("out2 not recognised, should be either 'coin' or 'list'")
  }
  
  if(out2 == "list"){
    l
  } else if (out2 == "coin"){
    if(is.null(lev3)){
      coin[[lev1]][[lev2]] <- l
    } else {
      coin[[lev1]][[lev2]][[lev3]] <- l
    }
    coin
  }
  
}


# Write a named data set to coin
#
# Writes a data set to the coin, and performs some checks in the process.
#
# @param coin A coin class object
# @param x The data to write
# @param dset A character string for naming the data, e.g. `Raw`.
# @param quietly If `TRUE`, suppresses messages.
# @param ignore_class If `TRUE` ignores the class of the input (used for [new_coin()]).
#
# @examples
# #
#
# @return Updated coin
write_dset <- function(coin, x, dset, quietly = FALSE, ignore_class = FALSE){
  
  # checks
  if(!ignore_class){
    stopifnot(is.coin(coin))
  }
  stopifnot(is.character(dset),
            length(dset)==1,
            is.data.frame(x))
  
  # further checks
  if(is.null(x$uCode)){
    stop("Required col uCode not found in data set to write to coin.")
  }
  icodes <- names(x)[names(x) != "uCode"]
  not_numeric <- !(sapply(x[icodes], is.numeric))
  if(any(not_numeric)){
    stop("Non-numeric cols detected in data set to be written to coin (excluding uCode).")
  }
  if(any(icodes %nin% coin$Meta$Ind$iCode[coin$Meta$Ind$Type %in% c("Indicator", "Aggregate")])){
    stop("Names of columns do not correspond to entries in .$Meta$Ind$iCode in data to write to coin.")
  }
  
  # flag if dset exists
  dset_exists <- !is.null(coin$Data[[dset]])
  
  # write to coin
  coin$Data[[dset]] <- x
  
  if(!quietly){
    message("Written data set to .$Data$", dset)
    if(dset_exists){
      message("(overwritten existing data set)")
    }
  }
  
  coin
  
}

#' Regenerate a coin
#'
#' Regenerates the `.$Data` entries in a coin by rerunning the construction functions according to the specifications in `.$Log`.
#' This effectively regenerates the results. Different variations of coins can be quickly achieved by editing the
#' saved arguments in `.$Log` and regenerating.
#'
#' The `from` argument allows partial regeneration, starting from a
#' specified function. This can be helpful to speed up regeneration in some cases. However, keep in mind that
#' if you change a `.$Log` argument from a function that is run before the point that you choose to start running
#' from, it will not affect the results.
#'
#' Note that while sets of weights will be passed to the regenerated COIN, anything in `.$Analysis` will be removed
#' and will have to be recalculated.
#'
#' See also `vignette("adjustments")` for more info on regeneration.
#'
#' @param x A coin class object
#' @param from Optional: a construction function name. If specified, regeneration begins from this function, rather
#' than re-running all functions.
#' @param quietly If `TRUE` (default), messages are suppressed during building.
#' @param ... arguments passed to or from other methods.
#'
#' @examples
#' # build full example coin
#' coin <- build_example_coin(quietly = TRUE)
#'
#' # copy coin
#' coin2 <- coin
#'
#' # change to prank function (percentile ranks)
#' # we don't need to specify any additional parameters (f_n_para) here
#' coin2$Log$Normalise$global_specs <- list(f_n = "n_prank")
#'
#' # regenerate
#' coin2 <- Regen(coin2)
#'
#' # compare index, sort by absolute rank difference
#' compare_coins(coin, coin2, dset = "Aggregated", iCode = "Index",
#'               sort_by = "Abs.diff", decreasing = TRUE)
#'
#' @return Updated coin object with regenerated results (data sets).
#'
#' @export
Regen.coin <- function(x, from = NULL, quietly = TRUE, ...){
  
  coin <- x
  
  stopifnot(is.coin(coin))
  
  
  # GATHER PARAMS -----------------------------------------------------------
  
  # the full list of function arguments, for each build_ function
  f_logs <- coin$Log
  f_names <- names(f_logs)
  
  # check if can regenerate
  stopifnot(!is.null(f_logs$can_regen))
  if(!f_logs$can_regen){
    stop("Cannot regenerate coin. This may be because it has been normalised with global = TRUE, or
         it has been converted from an older COIN class.")
  }
  # remove can_regen from here
  f_names <- setdiff(f_names, "can_regen")
  f_logs <- f_logs[f_names]
  
  # here we exclude any function names that are before "from", if it is specified
  if(!is.null(from)){
    if(from %nin% f_names){
      stop("Function name specified by 'from' is not found in the coin log.")
    }
    i_name <- which(f_names == from) - 1
    if(i_name > 0){
      f_names <- f_names[-(1:i_name)]
    }
  }

  # RERUN FUNCS -------------------------------------------------------------
  
  # looping over build_ functions
  for (func in f_names){
    
    # the arguments of the same func, stored in Log
    f_log <- f_logs[[func]]
    
    # the declared arguments of the function
    # NOTE this doesn't work here since construction funcs are now methods, so args are all (x, ...)
    #f_args <- names(formals(func))
    # check if what is in Log agrees with function arguments
    # if(!all(names(f_log) %in% f_args)){
    #   stop(paste0("Mismatch between function arguments of ", func, " and .$Log entry. Cannot regenerate."))
    # }
    
    # run function at arguments
    if(func == "new_coin"){
      
      if(quietly){
        coin <- suppressMessages( do.call(func, args = f_log) )
      } else {
        coin <- do.call(func, args = f_log)
      }
      
      # we also need to pass old weights to new coin
      wlist_old <- x$Meta$Weights[names(x$Meta$Weights) != "Original"]
      
      # the only thing to check is whether the iCodes are the same. If not, means that something has happened
      same_codes <- sapply(wlist_old, function(w){
        setequal(coin$Meta$Weights$Original$iCode, w$iCode)
      })

      if(any(!same_codes)){
        message(".$Meta$Weights iCodes do not match new coin, filtering them down")
        wlist_old <- lapply(wlist_old, function(x) x[iCode %in% unique(coin$Meta$Weights$Original$iCode)])
        coin$Meta$Weights <- c(coin$Meta$Weights, wlist_old)
      } else {
        coin$Meta$Weights <- c(coin$Meta$Weights, wlist_old)
      }
      
    } else {
      # add coin obj to arg list (not logged for obvious inception reasons)
      if(quietly){
        coin <- suppressMessages( do.call(func, args = c(list(x = coin), f_log) ) )
      } else {
        coin <- do.call(func, args = c(list(x = coin), f_log) )
      }
      
    }
    
  }
  
  
  # WEIGHTS -----------------------------------------------------------------
  
  
  
  coin
}

get_results_custom <-
function (coin, dset, tab_type = "Summ", also_get = NULL, use = "scores", 
          order_by = NULL, nround = 2, use_group = NULL, out2 = "df") 
{
  stopifnot(tab_type %in% c("Summ", "Aggs", "Full"), use %in% 
              c("scores", "ranks", "groupranks"), is.numeric(nround), 
            out2 %in% c("df", "coin"))
  also_get <- union(use_group, also_get)
  iData <- get_data(coin, dset = dset, also_get = also_get, 
                    use_group = use_group)
  mcols <- extract_iData(coin, iData, GET = "mCodes")
  iMeta <- coin$Meta$Ind
  iMeta_ia <- iMeta[iMeta$Type %in% c("Indicator", "Aggregate"), 
  ]
  iMeta_ia <- iMeta_ia[order(-iMeta_ia$Level, iMeta_ia$Parent), 
  ]
  #removing this for now until i find a better way to handle flattening the hierarchy in aggregation
  #TODO
  # if (any(iMeta_ia$iCode %nin% names(iData))) {
  #   stop("The data set extracted by 'dset' does not seem to be an aggregated data set (indicator or aggregate codes are missing).")
  # }
  iData <- iData[c(mcols, iMeta_ia$iCode)]
  if (is.null(order_by)) {
    sortcode <- iMeta_ia$iCode[iMeta_ia$Level == coin$Meta$maxlev]
  }
  else {
    if (order_by %nin% names(iData)) {
      stop("'order_by' is not found in the selected data set.")
    }
    sortcode <- order_by
  }
  iData$Rank <- rank(-1 * iData[[sortcode]], na.last = "keep", 
                     ties.method = "min")
  if (tab_type %in% c("Summ", "Summary")) {
    tabout <- iData[c(mcols, sortcode, "Rank")]
  }
  else if (tab_type %in% c("Aggs", "Aggregates")) {
    tabout <- iData[c(mcols, "Rank", iMeta_ia$iCode[iMeta_ia$Type == 
                                                      "Aggregate"])]
  }
  else if (tab_type %in% c("Full", "FullWithDenoms")) {
    othercodes <- coin$Meta$Lineage[[1]]
    stopifnot(any(othercodes %in% names(iData)))
    tabout <- iData[c(mcols, "Rank", iMeta_ia$iCode[iMeta_ia$Type == 
                                                      "Aggregate"], othercodes)]
  }
  tabout <- tabout[order(-tabout[[sortcode]]), ]
  tabout <- round_df(tabout, nround)
  if (use == "ranks") {
    tabout <- tabout[colnames(tabout) != "Rank"]
    tabout <- rank_df(tabout)
  }
  else if (use == "groupranks") {
    if (is.null(use_group)) {
      stop("If groupranks is specified, you need to also specify use_group.")
    }
    tabout <- tabout[colnames(tabout) != "Rank"]
    tabout <- rank_df(tabout, use_group = use_group)
    tabout <- tabout[order(tabout[[use_group]]), ]
  }
  if (out2 == "df") {
    return(tabout)
  }
  else if (out2 == "coin") {
    if (use == "scores") {
      coin$Results[[paste0(tab_type, "Score")]] <- tabout
    }
    else if (use == "ranks") {
      coin$Results[[paste0(tab_type, "Rank")]] <- tabout
    }
    else if (use == "groupranks") {
      coin$Results[[paste0(tab_type, "GrpRnk", use_group)]] <- tabout
    }
    return(coin)
  }
  else {
    stop("out2 not recognised!")
  }
}

get_PCA <- function(coin, dset = "Raw", iCodes = NULL, Level = NULL, by_groups = TRUE, imputed=F,
                    nowarnings = FALSE, weights_to = NULL, out2 = "list"){
  
  if(is.null(Level)){
    Level <- 1
  }
  
  # There is a catch here because we might want to do PCA weights across one level, but that level
  # may have multiple groups. This means we have to call PCA separately for each group.
  
  # first we define a function which returns weights for a given set of indicator data
  # this function implicitly calls other variables from the environment inside getPCA() so we don't need
  # to explicitly pass everything to it.
  PCAwts <- function(icodes1){
    
    # get ind data
    iData_ <- get_data(coin, dset = dset, iCodes = icodes1, Level = Level, also_get = "none")
    
    #impute data
    if(imputed) iData_ <- Impute(iData_, f_i='i_median', impute_by='column')
    
    # check for missing vals
    nNA <- sum(is.na(iData_))
    
    # remove any rows with missing data
    if (nNA > 0){
      dat4PCA <- stats::na.omit(iData_)
      if(!nowarnings){
        warning(paste0(nNA, " missing values found. Removing ", nrow(iData_)-nrow(dat4PCA), " rows with missing values in order to perform
PCA. You can also try imputing data first to avoid this."))
      }
      
    } else {
      dat4PCA <- iData_
    }
    
    # perform PCA
    PCAres <- stats::prcomp(dat4PCA, center = TRUE, scale = TRUE)
    
    # just for writing results - if Level not specified then we are working at ind level
    if(is.null(Level)){Level<-1}
    
    # weight from first PC should be the max variance weights
    wts <- as.numeric(PCAres$rotation[,1])

    #also make a weighted prediction from any PCs that have eigenvalue>1 (kaiser criterion)
    pca_kaiser <- PCAres$sd %>% 
      .^2 %>% 
      .[which(.>1)] %>% 
      as.data.table %>% 
      prop.table %>% 
      as.data.table %>% 
      .[, variable := paste0('PC', (.I))] %>% 
      setnames('.', 'wt') 
    
    preds <- 
    predict(PCAres, newdata=dat4PCA) %>% 
      .[,1:nrow(pca_kaiser)] %>%  #note here is where we could keep additional PCAs to permute the aggs
      as.data.table %>% 
      .[, id := .I] %>% 
      melt(id.vars='id') %>% 
      merge(pca_kaiser, by='variable') %>% 
      .[, combined := sum(value*wt), by=id] %>% #from cutter 2003 - SoVI analysis
      unique(by='id') %>% 
      .[, .(combined)] %>% 
      setnames('combined', coin$Meta$Lineage$Index %>% unique)

    #check direction of the PCAs association with inputs
    #sometimes PCA can be arbitrarily inverted
    sign <- cbind(dat4PCA, preds) %>% 
      cor(method='spearman') %>% 
      as.data.table %>% 
      .[ehd_rank!=1, sum(ehd_rank)] %>% 
      sign
    
    #invert if it has the wrong direction
    preds[, (coin$Meta$Lineage$Index %>% unique) := get(coin$Meta$Lineage$Index %>% unique)*sign]
    
    list(wts = wts, PCAres = PCAres, iCodes = names(iData_), preds=preds)
  }
  
  # We need to know the codes of the inds/aggs to get weights from
  iData_full <- get_data(coin, dset = dset, iCodes = iCodes, Level = Level, also_get = "none")
  IndCodes <- names(iData_full)
  
  if(by_groups){
    # OK, first thing is to find what groups we have
    # Get index structure
    lin <- coin$Meta$Lineage
    # Get cols of interest: the present one plus the parents
    lin <- lin[c(Level, Level + 1)]
    # Get parents of these codes
    parents <- unlist(unique(lin[(lin[[1]] %in% IndCodes) ,2]))
  } else {
    parents = "All"
  }
  
  
  # Right, now we need to cycle through these groups and do PCA on each group.
  # List for general PCA results
  PCAlist <- vector(mode = "list", length = length(parents))
  # copy of weights to modify
  wlist <- coin$Meta$Weights$Original
  
  for (ii in 1: length(parents)){
    if(by_groups){
      # get PCA results for group
      outPCA <- PCAwts(parents[ii])
    } else {
      # get PCA results for group
      outPCA <- PCAwts(NULL)
    }
    # attach weights to list
    # wts should be in the same order as out$iCodes. We have to make sure they match exactly here as
    # sometimes things get reordered. This is done with match() rather than %in% for this reason.
    wlist$Weight[match(outPCA$iCodes, wlist$iCode)] <- outPCA$wts
    # add general results to list
    PCAlist[[ii]] <- outPCA
  }
  # rename list
  names(PCAlist) <- parents
  
  # write results
  if(out2 == "coin"){
    
    if(!is.null(weights_to)){
      #w_name <- paste0("PCA_",dset,"L",Level)
      # write weights
      coin$Meta$Weights[[weights_to]] <- wlist
      message("Weights written to .$Meta$Weights$", weights_to)
    }
    
    # write other info
    coin$Analysis[[dset]][[paste0("$PCA$L",Level)]] <- PCAlist
    
    coin
    
  } else if (out2=='preds') {

    return(outPCA$preds)

    
  } else {
    
    list("Weights" = wlist,
         "PCAresults" = PCAlist)
  }
  
}

i_EM <- function(x){
  # impute
  amOut <- Amelia::amelia(x, m = 1, p2s = 0, boot.type = "none")
  # return imputed data
  amOut$imputations[[1]]
}
