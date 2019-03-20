

# Description:
# Bayesian inference
runBayesianInference <- function(gt.data,
                                 mcmc.chains,
                                 mcmc.steps,
                                 mcmc.warmup,
                                 cores,
                                 stan.model,
                                 dot.param) {

  # extra conversion
  Yd <- c()
  if(gt.data$Ntd > 0) {
    for(i in 1:gt.data$Ntd) {
      Yd <- cbind(Yd, as.numeric(gt.data$Yd[, i]))
    }
    gt.data$Yd <- Yd
  }


  data.list <- list(N = gt.data$N,
                    Ns = gt.data$Ns,
                    Nk = gt.data$Nk,
                    Ntq = gt.data$Ntq,
                    Ntd = gt.data$Ntd,
                    Yq = gt.data$Yq,
                    Yd = gt.data$Yd,
                    X = gt.data$X,
                    M = gt.data$Ms,
                    K = gt.data$K)


  # get initial parameter values
  control <- list(adapt_delta = dot.param$adapt_delta,
                  max_treedepth = dot.param$max_treedepth)

  # create names for sample files to use if sample.file is specified
  sample.file <- NULL
  if(dot.param$with.sample.file == TRUE) {
    sample.file <- paste(stan.model$model_name, "posterior", sep = '.')
  }

  # run
  posterior <- rstan::sampling(object = stan.model,
                               data = data.list,
                               iter = mcmc.steps,
                               warmup = mcmc.warmup,
                               chains = mcmc.chains,
                               cores = cores,
                               control = control,
                               verbose = dot.param$verbose,
                               refresh = dot.param$refresh,
                               sample_file = sample.file)

  # collect results from csv files if sample.file is specified
  if(dot.param$with.sample.file == TRUE) {
    posterior.files <- paste(sample.file, 1:mcmc.chains, sep = '_')
    if(all(file.exists(posterior.files)) == TRUE) {
      posterior <- rstan::read_stan_csv(csvfiles = posterior.files)
    }
    else {
      posterior.files <- paste(tempdir(), posterior.files, sep = '/')
      if(all(file.exists(posterior.files)) == TRUE) {
        posterior <- rstan::read_stan_csv(csvfiles = posterior.files)
      }
      else {
        stop("Cannot find posterior files")
      }
    }
  }

  # return
  return (posterior)
}
