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

  # Function:
  # Given model.name, find appropriate stan model
  getStanModel <- function(model.name) {
    # M0
    if(model.name == "M0") {
      model <- stanmodels$M0_loglik
    }
    if(model.name == "M0c") {
      model <- stanmodels$M0c_loglik
    }

    # M1
    if(model.name == "M1") {
      model <- stanmodels$M1_loglik
    }
    if(model.name == "M1c") {
      model <- stanmodels$M1c_loglik
    }

    # M2
    if(model.name == "M2") {
      model <- stanmodels$M2_loglik
    }
    if(model.name == "M2c") {
      model <- stanmodels$M2c_loglik
    }

    return (model)
  }


  # TODO:
  # here convert genotype.matrix & phenotype matrix into data.list


  # TODO:
  # how about covariates
  # Cov matrix NxC, C can be 0

  p <- vector(mode = "vector", length = length(models))

  # run each model
  for(i in 1:length(models)) {
    model <- getStanModel(model.name = models[i])

    # collect fits
    p[[i]] <- runStanModel(model = model, data.list = data.list)

  }

  names(p) <- models
  return (p)
}


