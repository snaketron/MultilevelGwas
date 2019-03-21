


# Function:
# Main analysis
runMLgwas <- function(genotype,
                      traits,
                      trait.type,
                      strains,
                      model,
                      mcmc.chains = 2,
                      mcmc.steps = 2500,
                      mcmc.warmup = 500,
                      cores = 1,
                      hdi.level = 0.95,
                      stat.learn.method = "rf",
                      cv.steps = 1000,
                      ...) {


  # check optional (dot) inputs
  dot.param <- checkDotParRun(...)


  # check inputs
  checkInputMain(genotype = genotype,
                 traits = traits,
                 trait.type = trait.type,
                 strains = strains,
                 model = model,
                 mcmc.chains = mcmc.chains,
                 mcmc.steps = mcmc.steps,
                 mcmc.warmup = mcmc.warmup,
                 cores = cores,
                 hdi.level = hdi.level,
                 stat.learn.method = stat.learn.method,
                 cv.steps = cv.steps)


  # convert input data to genetic-trait data (for stan)
  gt.data <- getStanData(genotype = genotype,
                         traits = traits,
                         trait.type = trait.type,
                         strains = strains)


  # TODO: uncomment in final version
  # get the appropriate stan model
  # stan.model <- getStanModel(model.name = model,
  #                            comparison.mode = FALSE)
  stan.model <- getStanModelDebug(model.name = model)


  cat("======== Bayesian Inference ======== \n")
  p <- runBayesianInference(gt.data = gt.data,
                            mcmc.chains = mcmc.chains,
                            mcmc.steps = mcmc.steps,
                            mcmc.warmup = mcmc.warmup,
                            cores = cores,
                            stan.model = stan.model,
                            dot.param = dot.param)



  # cat("======== Statistical Learning ======== \n")
  # s <- runStatLearn(gt.data = gt.data,
  #                   method = stat.learn.method,
  #                   cv.steps = cv.steps,
  #                   hdi.level = hdi.level,
  #                   cores = cores,
  #                   dot.param = dot.param)


  # cat("======== Collecting Results ======== \n")
  # o <- getScores(p = p, s = s$results,
  #                hdi.level = hdi.level,
  #                gt.data = gt.data)
  #
  #
  # # format scores
  # o <- do.call(rbind, o)
  # o <- o[, c("site", "ref", "alt", "refN", "altN", "p", "mean",
  #            "se_mean", "sd", "X2.5.", "X97.5.", "n_eff", "Rhat",
  #            "ca", "ca.L", "ca.H", "k", "k.L", "k.H")]
  # colnames(o) <- c("site", "ref", "alt", "refN", "altN", "phenotype.id",
  #                  "beta.mean", "beta.se", "beta.sd", "beta.hdi.low",
  #                  "beta.hdi.high", "Neff", "Rhat",
  #                  "ca.mean", "ca.hdi.low", "ca.hdi.high",
  #                  "kappa.mean", "kappa.hdi.low", "kappa.hdi.high")
  #
  #
  #
  # cat("======== Pareto Optimization ======== \n")
  # o <- getParetoScores(scores = o)
  #
  #
  #
  # # ppc
  # cat("======== Posterior Prediction ======== \n")
  # ppc <- getPpc(posterior = p$posterior,
  #               gt.data = gt.data,
  #               hdi.level = hdi.level)

  return(p)

  return (list(scores = o,
               ppc = ppc,
               complete.posterior = p$posterior))
}

