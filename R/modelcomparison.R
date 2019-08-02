# Description:
# Given a set of model names (M1-M3, M1*-M3*),
# run each and return a list of posteriors
runModelComparison <- function(genotype,
                               traits,
                               trait.type,
                               strains,
                               models,
                               mcmc.chains = 2,
                               mcmc.steps = 1500,
                               mcmc.warmup = 500,
                               mcmc.cores = 1,
                               hdi.level = 0.95) {


  # in case the user provides duplicate models
  models <- unique(models)


  # check inputs
  checkInputModelComparison(genotype = genotype,
                            traits = traits,
                            trait.type = trait.type,
                            strains = strains,
                            models = models,
                            mcmc.chains = mcmc.chains,
                            mcmc.steps = mcmc.steps,
                            mcmc.warmup = mcmc.warmup,
                            mcmc.cores = mcmc.cores,
                            hdi.level = hdi.level)


  # convert input data to genetic-trait data (for stan)
  gt.data <- getStanData(genotype = genotype,
                         traits = traits,
                         trait.type = trait.type,
                         strains = strains)


  ps <- vector(mode = "list", length = length(models))
  names(ps) <- models


  cat("1) Bayesian inference ... \n")
  # run each model
  for(i in 1:length(models)) {

    # TODO: uncomment
    # stan.model <- getStanModel(model.name = models[i],
    #                            comparison.mode = TRUE)
    stan.model <- getStanModelDebug(model.name = models[i])


    # collect fits
    ps[[i]] <- runBayesianInference(gt.data = gt.data,
                                    mcmc.chains = mcmc.chains,
                                    mcmc.steps = mcmc.steps,
                                    mcmc.warmup = mcmc.warmup,
                                    cores = cores,
                                    stan.model = stan.model,
                                    dot.param = dot.param,
                                    comparison = TRUE)
  }


  # compute LOOIC
  cat("2) Computing LOOIC ... \n")
  ic <- NA
  # ic <- getIC(ps = ps)



  # compute PPC
  cat("3) Posterior prediction ... \n")
  ppc <- NA
  # ppc <- getPpc(ps = ps,
  #               gt.data = gt.data,
  #               models = models,
  #               hdi.level = hdi.level)
  # ppc <- getPpcMc(ps = ps,
  #                 gt.data = gt.data,
  #                 models = models,
  #                 hdi.level = hdi.level)


  return (list(ps = ps,
               ic = ic,
               ppc = ppc,
               gt.data = gt.data))
}



# loo information criterion
getIC <- function(ps) {

  # loo and waic
  loo.list <- vector(mode = "list", length = length(ps))
  waic.list <- vector(mode = "list", length = length(ps))
  names(loo.list) <- names(ps)
  names(waic.list) <- names(ps)

  for(i in 1:length(ps)) {
    # no need to declare loo at all
    loo.list[[i]] <- loo::loo(x = loo::extract_log_lik(stanfit = ps[[i]]))
    waic.list[[i]] <- loo::waic(x = loo::extract_log_lik(stanfit = ps[[i]]))
  }

  return(list(loo = loo.list,
              waic = waic.list))
}




# Function:
# Posterior prediction multicore
getPpcMc <- function(ps,
                     gt.data,
                     models,
                     hdi.level,
                     cores = 4,
                     debug = F) {





  # runPpcMono <- function(ext, gt.data,
  #                        model, hdi.level) {
  #
  #
  #   cat("Model:", model, "\n", sep = '')
  #
  #   # get posterior  #TODO: fix 100->500 or 1,000
  #   ext <- data.frame(ext)
  #   ext <- ext[sample(x = 1:nrow(ext),
  #                     size = min(c(250, nrow(ext))),
  #                     replace = TRUE), ]
  #
  #   ppc.out <- c()
  #   ppc.out <- rbind(ppc.out, getPpcVertical(ext = ext,
  #                                            gt.data = gt.data,
  #                                            model = model,
  #                                            hdi.level = hdi.level))
  #
  #   ppc.out <- rbind(ppc.out, getPpcHorizontal(ext = ext,
  #                                              gt.data = gt.data,
  #                                              model = model,
  #                                              hdi.level = hdi.level))
  #
  #   return (ppc.out)
  # }



  runPpcMono <- function(x,
                         gt.data,
                         hdi.level) {

    model <- x@model_name
    cat("Model:", model, "\n", sep = '')

    # get posterior  #TODO: fix 100->500 or 1,000
    x <- rstan::extract(object = x)
    x <- data.frame(x)
    x <- x[sample(x = 1:nrow(x),
                      size = min(c(250, nrow(x))),
                      replace = TRUE), ]

    ppc.out <- getPpcVertical(ext = x,
                              gt.data = gt.data,
                              model = model,
                              hdi.level = hdi.level)
    ppc.out <- rbind(ppc.out,
                     getPpcHorizontal(ext = x,
                                      gt.data = gt.data,
                                      model = model,
                                      hdi.level = hdi.level))

    return (ppc.out)
  }



  ppc.list <- lapply(X = ps,
                     FUN = runPpcMono,
                     gt.data = gt.data,
                     hdi.level = hdi.level)

  ppc.list <- parallel::mclapply(X = ps,
                                 FUN = runPpcMono,
                                 gt.data = gt.data,
                                 hdi.level = hdi.level,
                                 mc.cores = cores)


  # multicore classification
  # cl <- parallel::makeCluster(cores)
  # doParallel::registerDoParallel(cl)
  # if(debug) {
  #   ppc.list <- (foreach(i = 1:length(ps),
  #                        .export = c("getHdi",
  #                                    "getPpcVertical",
  #                                    "getPpcHorizontal")) %dopar%
  #                  runPpcMono(ext = ps[[i]],
  #                             gt.data = gt.data,
  #                             model = models[i],
  #                             hdi.level = hdi.level))
  # }
  # else {
  #   ppc.list <- (foreach(i = 1:length(ps),
  #                        .export = c("getHdi",
  #                                    "getPpcVertical",
  #                                    "getPpcHorizontal")) %dopar%
  #                  runPpcMono(ext = rstan::extract(object = ps[[i]],
  #                                                  include = FALSE,
  #                                                  pars = c("z", "log_lik",
  #                                                           "log_lik2")),
  #                             gt.data = gt.data,
  #                             model = models[i],
  #                             hdi.level = hdi.level))
  # }
  # top cluster
  # parallel::stopCluster(cl = cl)
  # doParallel::stopImplicitCluster()

  ppc.list <- do.call(rbind, ppc.list)
  return (ppc.list)
}



# Function:
# Posterior prediction
getPpc <- function(ps,
                   gt.data,
                   models,
                   hdi.level,
                   debug = F) {

  ppc.out <- c()
  for(i in 1:length(models)) {
    cat("Model:", models[i], "\n", sep = '')


    if(debug) {
      ext <- data.frame(ps[[i]])
    }
    else {
      # get posterior
      ext <- data.frame(rstan::extract(object = ps[[i]],
                                       pars = c("z", "log_lik", "log_lik2"),
                                       include = FALSE))
      ext <- ext[sample(x = 1:nrow(ext),
                        size = min(c(500, nrow(ext))),
                        replace = TRUE), ]
    }

    # ppc.out <- rbind(ppc.out, getPpcVertical(ext = ext,
    #                                          gt.data = gt.data,
    #                                          model = models[i],
    #                                          hdi.level = hdi.level))

    ppc.out <- rbind(ppc.out, getPpcHorizontal(ext = ext,
                                               gt.data = gt.data,
                                               model = models[i],
                                               hdi.level = hdi.level))

    cat("\n")
  }
  return (ppc.out)
}




