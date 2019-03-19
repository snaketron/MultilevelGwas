

# Description:
# Bayesian inference
runBayesianInference <- function(genphen.data,
                                 mcmc.chains,
                                 mcmc.steps,
                                 mcmc.warmup,
                                 cores,
                                 model.stan,
                                 ...) {

  # extra conversion
  Yd <- c()
  if(genphen.data$Ntd > 0) {
    for(i in 1:genphen.data$Ntd) {
      Yd <- cbind(Yd, as.numeric(genphen.data$Yd[, i]))
    }
    genphen.data$Yd <- Yd
  }


  data.list <- list(N = genphen.data$N,
                    Ntq = genphen.data$Ntq,
                    Ntd = genphen.data$Ntd,
                    Ns = genphen.data$Ns,
                    Nsk = genphen.data$Nsk,
                    Yq = genphen.data$Yq,
                    Yd = genphen.data$Yd,
                    X = genphen.data$X,
                    M = genphen.data$Ms)


  # get initial parameter values
  control <- list(adapt_delta = list(...)[["adapt_delta"]],
                  max_treedepth = list(...)[["max_treedepth"]])
  refresh <- list(...)[["refresh"]]
  verbose <- list(...)[["verbose"]]
  posterior <- rstan::sampling(object = model.stan,
                               data = data.list,
                               iter = mcmc.steps,
                               warmup = mcmc.warmup,
                               chains = mcmc.chains,
                               cores = cores,
                               control = control,
                               verbose = verbose,
                               refresh = refresh)

  # return
  return (list(posterior = posterior))
}
