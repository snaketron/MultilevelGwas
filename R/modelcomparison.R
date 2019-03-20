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
                               cores = 1,
                               hdi.level = 0.95,
                               ...) {


  # check optional (dot) inputs
  dot.param <- checkDotParModelComparison(...)


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
                            cores = cores,
                            hdi.level = hdi.level)


  # convert input data to genetic-trait data (for stan)
  gt.data <- getStanData(genotype = genotype,
                         traits = traits,
                         trait.type = trait.type,
                         strains = strains)


  ps <- vector(mode = "list", length = length(models))
  names(ps) <- models

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
                                    dot.param = dot.param)
  }


  # compute LOO_ic
  loo.ic <- getLooIc(ps = ps)


  # compute PPC
  ppc <- getPpc(ps = ps, gt.data = gt.data, models = models,
                hdi.level = hdi.level, cores = cores)


  return (list(ps = ps,
               loo.ic = loo.ic,
               ppc = ppc))
}



# TODO:
# loo information criterion
getLooIc <- function(ps) {

  loo.list <- vector(mode = "list", length = length(ps))
  names(loo.list) <- names(ps)

  for(i in 1:length(ps)) {
    # no need to declare loo at all
    # loo.list[[i]] <- loo::loo(loo::extract_log_lik(stanfit = ps[[i]]))
    loo.list[[i]] <- rstan::loo(ps[[i]], pars = "log_lik")
  }

  return(loo.list)
}


# Function:
# Posterior prediction
getPpc <- function(ps, gt.data, models,
                   hdi.level, cores) {




  ppc.list <- vector(mode = "list", length = length(ps))
  names(ppc.list) <- names(ps)


  for(i in 1:length(models)) {
    # get posterior
    ext <- data.frame(rstan::extract(object = ps[[i]]))
    ext <- ext[, regexpr(pattern = "z\\.|log_lik",
                         text = colnames(ext)) == -1]
    ext <- ext[sample(x = 1:nrow(ext),
                      size = max(c(500, nrow(ext))),
                      replace = TRUE), ]



    # register multicore ppc
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)


    if(gt.data$Ntq > 0 & gt.data$Ntd > 0) {
      ppc.list[[i]] <- (foreach(s = 1:gt.data$Ns,
                                .export = c("getHdi", "getPpcQD")) %dopar%
                          getPpcQD(ext = ext,
                                   gt.data = gt.data,
                                   model = models[i],
                                   hdi.level = hdi.level,
                                   s = s))
    }
    else if(gt.data$Ntq > 0 & gt.data$Ntd == 0) {
      ppc.list[[i]] <- (foreach(s = 1:gt.data$Ns,
                                .export = c("getHdi", "getPpcQ")) %dopar%
                          getPpcQ(ext = ext,
                                  gt.data = gt.data,
                                  model = models[i],
                                  hdi.level = hdi.level,
                                  s = s))
    }
    else if(gt.data$Ntq == 0 & gt.data$Ntd > 0) {
      ppc.list[[i]] <- (foreach(s = 1:gt.data$Ns,
                                .export = c("getHdi", "getPpcD")) %dopar%
                          getPpcD(ext = ext,
                                  gt.data = gt.data,
                                  model = models[i],
                                  hdi.level = hdi.level,
                                  s = s))
    }

    # release cluster
    parallel::stopCluster(cl = cl)
    doParallel::stopImplicitCluster()
  }




  return (ppc.list)
}





