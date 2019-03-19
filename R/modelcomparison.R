# Description:
# Given a set of model names (M1-M3, M1*-M3*),
# run each and return a list of posteriors
compareModels <- function(models,
                          genotype.matrix,
                          phenotype.matrix) {

  # Function:
  # Run stan for selected model, data, and mcmc param
  # TODO: default = dots
  runStanModel <- function(model, data.list, ...) {
    p <- rstan::sampling(object = model,
                          data = data.list,
                          iter = 1500,
                          warmup = 500,
                          chains = 4,
                          cores = 4,
                          control = list(adapt_delta = 0.99,
                                         max_treedepth = 10))
    return (p)
  }


  # TODO:
  # here convert genotype.matrix & phenotype matrix into data.list


  # TODO:
  # how about covariates
  # Cov matrix NxC, C can be 0

  p <- vector(mode = "vector", length = length(models))

  # run each model
  for(i in 1:length(models)) {
    model <- getStanModel(model.name = models[i],
                          comparison.mode = TRUE)

    # collect fits
    p[[i]] <- runStanModel(model = model, data.list = data.list)

  }

  names(p) <- models
  return (p)
}



# TODO:
# loo operation here
