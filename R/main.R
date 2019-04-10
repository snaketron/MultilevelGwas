


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
  stan.model <- getStanModelDebug(model.name = model,
                                  comparison = FALSE)



  cat("======== Bayesian Inference ======== \n")
  p <- runBayesianInference(gt.data = gt.data,
                            mcmc.chains = mcmc.chains,
                            mcmc.steps = mcmc.steps,
                            mcmc.warmup = mcmc.warmup,
                            cores = cores,
                            stan.model = stan.model,
                            dot.param = dot.param,
                            comparison = FALSE)



  cat("======== Statistical Learning ======== \n")
  s <- runStatLearn(gt.data = gt.data,
                    method = stat.learn.method,
                    cv.steps = cv.steps,
                    hdi.level = hdi.level,
                    cores = cores,
                    dot.param = dot.param)



  cat("======== Pareto Optimization ======== \n")
  r <- getParetoRanks(p = p, s = s,
                      model = model,
                      hdi.level = hdi.level)



  cat("======== Collecting Results ======== \n")
  scores <- getScores(p = p, s = s, r = r,
                      model = model,
                      gt.data = gt.data,
                      hdi.level = hdi.level)



  cat("======== Posterior Predictive Checks ======== \n")
  ppc <- getPpc(ps = list(p = p),
                gt.data = gt.data,
                models = model,
                hdi.level = hdi.level)



  return (list(scores = scores,
               ppc = ppc,
               gt.data = gt.data))
}

