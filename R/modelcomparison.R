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
  ic <- getIC(ps = ps)



  # compute PPC
  cat("3) Posterior prediction ... \n")
  ppc <- getPpc(ps = ps,
                gt.data = gt.data,
                models = models,
                hdi.level = hdi.level)


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

    # loo.list[[i]] <- rstan::loo(ps[[i]], pars = "log_lik")
  }

  return(list(loo = loo.list,
              waic = waic.list))
}




# Function:
# Posterior prediction
getPpc <- function(ps,
                   gt.data,
                   models,
                   hdi.level) {


  ppc.list <- vector(mode = "list", length = length(ps))
  names(ppc.list) <- names(ps)


  for(i in 1:length(models)) {
    # get posterior
    ext <- data.frame(rstan::extract(object = ps[[i]]))
    ext <- ext[, regexpr(pattern = "z\\.|log_lik",
                         text = colnames(ext)) == -1]
    ext <- ext[sample(x = 1:nrow(ext),
                      size = min(c(500, nrow(ext))),
                      replace = TRUE), ]


    ppc.out <- c()
    for(s in 1:gt.data$Ns) {

      if(gt.data$Ntq > 0 & gt.data$Ntd > 0) {
        ppc.out <- rbind(ppc.out, getPpcQD(ext = ext,
                                           gt.data = gt.data,
                                           model = models[i],
                                           hdi.level = hdi.level,
                                           s = s))
      }
      else if(gt.data$Ntq > 0 & gt.data$Ntd == 0) {
        ppc.out <- rbind(ppc.out, getPpcQ(ext = ext,
                                          gt.data = gt.data,
                                          model = models[i],
                                          hdi.level = hdi.level,
                                          s = s))
      }
      else if(gt.data$Ntq == 0 & gt.data$Ntd > 0) {
        ppc.out <- rbind(ppc.out, getPpcD(ext = ext,
                                          gt.data = gt.data,
                                          model = models[i],
                                          hdi.level = hdi.level,
                                          s = s))
      }

      if(models[i] != "M0" & models[i] != "M0c") {
        ppc.out <- rbind(ppc.out, getPpcSnpBeta(p = ext, gt.data = gt.data,
                                                s = s, hdi.level = hdi.level))

      }

      # progress
      if(s %% 5 == 0) {
        cat(models[i], ":", s, "/", gt.data$Ns, "\n", sep = '')
      }

    }
    ppc.list[[i]] <- ppc.out
  }

  # release cluster
  # parallel::stopCluster(cl = cl)
  # doParallel::stopImplicitCluster()

  return (ppc.list)
}





