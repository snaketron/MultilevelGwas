

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
                    K = gt.data$K,
                    Xk = gt.data$Xk)


  # get appropriate parameters to monitor
  pars <- getStanModelPars(model.name = stan.model@model_name)



  # sample file
  sample.file <- paste(round(x = runif(n=1,min=0,max=10^6), digits=0),
                       "sampling", stan.model@model_name, sep = '_')


  # run
  glm <- rstan::sampling(object = stan.model,
                         data = data.list,
                         iter = mcmc.steps,
                         warmup = mcmc.warmup,
                         chains = mcmc.chains,
                         cores = mcmc.cores,
                         control = list(adapt_delta = adapt.delta,
                                        max_treedepth = max.treedepth),
                         pars = pars,
                         include = TRUE,
                         sample_file = sample.file,
                         refresh = 100,
                         algorithm = "NUTS")

  files <- paste(paste(sample.file, "_", 1:mcmc.chains, sep = ''),".csv",sep='')
  return (list(glm = glm, sample.files = files))
}
