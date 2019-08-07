


# Description:
# Main analysis
runGwas <- function(genotype,
                    traits,
                    trait.type,
                    strains,
                    model,
                    mcmc.warmup = 500,
                    mcmc.steps = 1500,
                    mcmc.chains = 4,
                    mcmc.cores = 4,
                    hdi.level = 0.95,
                    adapt.delta = 0.95,
                    max.treedepth = 10,
                    write.samples = TRUE) {


  # check inputs
  checkInput(genotype = genotype,
             traits = traits,
             trait.type = trait.type,
             strains = strains,
             model = model,
             mcmc.chains = mcmc.chains,
             mcmc.steps = mcmc.steps,
             mcmc.warmup = mcmc.warmup,
             mcmc.cores = mcmc.cores,
             hdi.level = hdi.level,
             adapt.delta = adapt.delta,
             max.treedepth = max.treedepth,
             write.samples = write.samples)


  # convert input data to genetic-trait data (for stan)
  gt.data <- getStanData(genotype = genotype,
                         traits = traits,
                         trait.type = trait.type,
                         strains = strains)


  # compile stan model
  cat("Compiling model ... \n")
  stan.model <- getStanModel(model.name = model)


  cat("Running inference ... \n")
  glm <- runInference(gt.data = gt.data,
                      mcmc.chains = mcmc.chains,
                      mcmc.steps = mcmc.steps,
                      mcmc.warmup = mcmc.warmup,
                      mcmc.cores = mcmc.cores,
                      stan.model = stan.model,
                      adapt.delta = adapt.delta,
                      max.treedepth = max.treedepth,
                      comparison = FALSE,
                      write.samples = write.samples)


  cat("Compute scores ... \n")
  scores <- getScores(glm = glm$glm,
                      model = model,
                      hdi.level = hdi.level,
                      gt.data = gt.data)


  return (list(glm = glm, scores = scores))
}




# Description:
# Comparison analysis
runComparison <- function(genotype,
                          traits,
                          trait.type,
                          strains,
                          models,
                          mcmc.warmup = 500,
                          mcmc.steps = 1500,
                          mcmc.chains = 4,
                          mcmc.cores = 4,
                          hdi.level = 0.95,
                          adapt.delta = 0.95,
                          max.treedepth = 10,
                          write.samples = TRUE) {

  # check inputs
  checkInput(genotype = genotype,
             traits = traits,
             trait.type = trait.type,
             strains = strains,
             model = "M0",
             mcmc.chains = mcmc.chains,
             mcmc.steps = mcmc.steps,
             mcmc.warmup = mcmc.warmup,
             mcmc.cores = mcmc.cores,
             hdi.level = hdi.level,
             adapt.delta = adapt.delta,
             max.treedepth = max.treedepth,
             write.samples = write.samples)


  # check models input
  models <- unique(models)
  if(length(models) <= 1) {
    stop("at least two models must be specified.")
  }
  if(all(models %in% c("M0", "M0c", "M1", "M1c")) == FALSE) {
    stop("allowed models are 'M0', 'M0c', 'M1' and 'M1c'.")
  }


  # convert input data to genetic-trait data (for stan)
  gt.data <- getStanData(genotype = genotype,
                         traits = traits,
                         trait.type = trait.type,
                         strains = strains)

  # results
  out <- vector(mode = "list", length = length(models))
  names(out) <- models

  # scores
  scores <- vector(mode = "list", length = length(models))
  names(scores) <- models

  # ppc
  ppc <- vector(mode = "list", length = length(models))
  names(scores) <- models



  # loop through models
  for(model in models) {

    # compile stan model
    cat(paste("Compiling model (", model, ")", "..., \n", sep = ''))
    stan.model <- getStanModel(model.name = model)

    # inference
    cat(paste("Running inference (", model, ")", "..., \n", sep = ''))
    out[[model]] <- runInference(gt.data = gt.data,
                                 mcmc.chains = mcmc.chains,
                                 mcmc.steps = mcmc.steps,
                                 mcmc.warmup = mcmc.warmup,
                                 mcmc.cores = mcmc.cores,
                                 stan.model = stan.model,
                                 adapt.delta = adapt.delta,
                                 max.treedepth = max.treedepth,
                                 comparison = TRUE)

    cat(paste("Computing scores (", model, ")", "..., \n", sep = ''))
    scores[[model]] <- getScores(glm = out[[model]]$glm,
                                 model = model,
                                 hdi.level = hdi.level,
                                 gt.data = gt.data)

  }

  return (list(out = out, scores = scores))
}


