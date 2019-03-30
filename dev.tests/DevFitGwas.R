source("R/main.R")
source("R/util.R")
source("R/bayesian.R")
source("R/statlearn.R")
source("R/modelcomparison.R")
source("R/ppc.R")

require(rstan)
require(loo)
require(parallel)
require(foreach)
require(doParallel)
require(Biostrings)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)



fit.m0 <- rstan::stan_model(file = "src/stan_files/fitM0_loglik.stan")
fit.m0c <- rstan::stan_model(file = "src/stan_files/fitM0c_loglik.stan")
fit.m1 <- rstan::stan_model(file = "src/stan_files/fitM1_loglik.stan")
fit.m1c <- rstan::stan_model(file = "src/stan_files/fitM1c_loglik.stan")
fit.mX <- rstan::stan_model(file = "src/stan_files/fitMX_loglik.stan")


# data.list <- get(load(file = "~/MiceGwas/output/posterior/vsv/data.list.RData"))
# data.list <- get(load(file = "~/MiceGwas/output/posterior/lcmv/data.list.RData"))
data.list <- get(load(file = "~/MiceGwas/output/posterior/tnf/data.list.RData"))
data.list$Ntq <- data.list$Ntc
data.list$Yq <- data.list$Yc
data.list$X <- NULL
data.list$Ntd <- 0
data.list$Yd <- matrix(data = 0, nrow = data.list$N, ncol = 0)




fitted.m0 <- rstan::sampling(object = fit.m0, data = data.list,
                             iter = 2500, warmup = 1500, chains = 4,
                             cores = 4, control = list(adapt_delta = 0.99,
                                                       max_treedepth = 15))

fitted.m0c <- rstan::sampling(object = fit.m0c, data = data.list,
                              iter = 2500, warmup = 1500, chains = 4,
                              cores = 4, control = list(adapt_delta = 0.99,
                                                        max_treedepth = 15))

fitted.m1 <- rstan::sampling(object = fit.m1, data = data.list,
                             iter = 2500, warmup = 1500, chains = 4,
                             cores = 4, control = list(adapt_delta = 0.99,
                                                       max_treedepth = 15))

fitted.m1c <- rstan::sampling(object = fit.m1c, data = data.list,
                              iter = 2500, warmup = 1500, chains = 4,
                              cores = 4, control = list(adapt_delta = 0.99,
                                                        max_treedepth = 15))

fitted.X <- rstan::sampling(object = fit.mX, data = data.list,
                              iter = 2500, warmup = 1500, chains = 4,
                              cores = 4, control = list(adapt_delta = 0.99,
                                                        max_treedepth = 10))



fit.comparison <- list(fitted.m0 = fitted.m0, fitted.m0c = fitted.m0c,
                       fitted.m1 = fitted.m1, fitted.m1c = fitted.m1c,
                       fitted.X = fitted.X)
save(fit.comparison, file = "dev.tests/fit.comparison.lcmv.RData")





loo::loo(loo::extract_log_lik(stanfit = fit.comparison$fitted.m0))
loo::loo(loo::extract_log_lik(stanfit = fit.comparison$fitted.m0c))
loo::loo(loo::extract_log_lik(stanfit = fit.comparison$fitted.m1))
loo::loo(loo::extract_log_lik(stanfit = fit.comparison$fitted.m1c))
loo::loo(loo::extract_log_lik(stanfit = fit.comparison$fitted.X))



loo::compare(loo::loo(loo::extract_log_lik(stanfit = fit.comparison$fitted.m0)),
             loo::loo(loo::extract_log_lik(stanfit = fit.comparison$fitted.m1)))

loo::compare(loo::loo(loo::extract_log_lik(stanfit = fit.comparison$fitted.m0c)),
             loo::loo(loo::extract_log_lik(stanfit = fit.comparison$fitted.m1c)))

loo::compare(loo::loo(loo::extract_log_lik(stanfit = fit.comparison$fitted.m0)),
             loo::loo(loo::extract_log_lik(stanfit = fit.comparison$fitted.m0c)))
loo::compare(loo::loo(loo::extract_log_lik(stanfit = fit.comparison$fitted.m1)),
             loo::loo(loo::extract_log_lik(stanfit = fit.comparison$fitted.m1c)))

# MU comparison
mu <- data.frame(summary(fitted.m0, par = "mu", prob = c((1-0.95)/2, 1-(1-0.95)/2))$summary)
mu$par <- rownames(mu)
mu$t <- NA
mu$t[regexpr(pattern = "\\[1", text = mu$par) != -1] <- 1
mu$t[regexpr(pattern = "\\[2", text = mu$par) != -1] <- 2
mu$t[regexpr(pattern = "\\[3", text = mu$par) != -1] <- 3
mu$t[regexpr(pattern = "\\[4", text = mu$par) != -1] <- 4

muc <- data.frame(summary(fitted.m0c, par = "mu", prob = c((1-0.95)/2, 1-(1-0.95)/2))$summary)
muc$par <- rownames(muc)
muc$t <- NA
muc$t[regexpr(pattern = "\\[1", text = muc$par) != -1] <- 1
muc$t[regexpr(pattern = "\\[2", text = muc$par) != -1] <- 2
muc$t[regexpr(pattern = "\\[3", text = muc$par) != -1] <- 3
muc$t[regexpr(pattern = "\\[4", text = muc$par) != -1] <- 4


plot(mu$mean[mu$t == 2], mu$mean[mu$t == 4])
points(muc$mean[muc$t == 2], muc$mean[muc$t == 4], col = "blue", type = "p")
abline(a = 0, b = 0.649660)


# RHO
summary(fitted.m0c, par = "rho", prob = c((1-0.95)/2, 1-(1-0.95)/2))$summary


summary(fitted.m0, par = "grand_mu", prob = c((1-0.95)/2, 1-(1-0.95)/2))$summary
summary(fitted.m0c, par = "grand_mu", prob = c((1-0.95)/2, 1-(1-0.95)/2))$summary



summary(fitted.m0, par = "sigma", prob = c((1-0.95)/2, 1-(1-0.95)/2))$summary
summary(fitted.m0c, par = "sigma", prob = c((1-0.95)/2, 1-(1-0.95)/2))$summary


summary(fitted.m1, par = "sigma", prob = c((1-0.95)/2, 1-(1-0.95)/2))$summary
summary(fitted.m1c, par = "sigma", prob = c((1-0.95)/2, 1-(1-0.95)/2))$summary



sigmas <- data.frame(summary(fitted.X, par = "sigma", prob = c((1-0.95)/2, 1-(1-0.95)/2))$summary)
sigmas$par <- rownames(sigmas)
sigmas$t <- NA
sigmas$t[regexpr(pattern = "\\[1", text = sigmas$par) != -1] <- 1
sigmas$t[regexpr(pattern = "\\[2", text = sigmas$par) != -1] <- 2
plot(sigmas$mean[sigmas$t == 1], sigmas$mean[sigmas$t == 2], col = "blue", type = "p")
abline(a = 0, b = 1)


