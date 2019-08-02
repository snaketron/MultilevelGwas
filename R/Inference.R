

# Description:
# Bayesian inference
runInference <- function(gt.data,
                         mcmc.chains,
                         mcmc.steps,
                         mcmc.warmup,
                         mcmc.cores,
                         stan.model,
                         adapt.delta,
                         max.treedepth,
                         comparison) {


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



  # get appropriate parameters to monitor
  pars <- getStanModelPars(model.name = stan.model@model_name,
                           comparison = comparison)

  # run
  p <- rstan::sampling(object = stan.model,
                       data = data.list,
                       iter = mcmc.steps,
                       warmup = mcmc.warmup,
                       chains = mcmc.chains,
                       cores = mcmc.cores,
                       control = list(adapt_delta = adapt.delta,
                                      max_treedepth = max.treedepth),
                       pars = pars,
                       include = TRUE,
                       refresh = 250)
  # sample_file = sample.file,
  return (p)
}
