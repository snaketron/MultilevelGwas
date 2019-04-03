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







# homozygous strains
N <- 5
K <- 10

genotype <- cbind(rep(x = c("A", "C"), each = N*K/2),
                  replicate(n = 24, rep(x = sample(x = c("A", "C"),
                                                   size = K, replace = T), each = N)))
unique(apply(X = genotype, MARGIN = 2, FUN = base::paste, collapse = ''))


rho.12 <- rnorm(n = K, mean = 1, sd = 0.1)
rho.13 <- rnorm(n = K, mean = 0.8, sd = 0.1)
rho.14 <- rnorm(n = K, mean = 1, sd = 0.1)
rho.12 <- sapply(X = rho.12, FUN = getMax <- function(x) {min(x, 1)})
rho.13 <- sapply(X = rho.13, FUN = getMax <- function(x) {min(x, 1)})
rho.14 <- sapply(X = rho.14, FUN = getMax <- function(x) {min(x, 1)})



mu1.k <- c(rnorm(n = K/2, mean = 2.5, sd = 1),
           rnorm(n = K/2, mean = 0, sd = 1))
mu2.k <- mu1.k * rho.12
mu3.k <- mu1.k * rho.13
mu4.k <- mu1.k * rho.14

library("Hmisc")
library("corrplot")
corrplot(rcorr(cbind(mu1.k, mu2.k, mu3.k, mu4.k))$r,
         type = "upper", order = "hclust",
         tl.col = "black", tl.srt = 45)

sigma.k <- abs(rnorm(n = K, mean = 0.5, sd = 1))

traits <- matrix(data = 0, nrow = N*K, ncol = 4)
for(i in 1:K) {
  traits[(N*(i-1)+1):(N*i), 1] <- rnorm(n = N, mean = mu1.k[i], sd = sigma.k[i])
  traits[(N*(i-1)+1):(N*i), 2] <- rnorm(n = N, mean = mu2.k[i], sd = sigma.k[i])
  traits[(N*(i-1)+1):(N*i), 3] <- rnorm(n = N, mean = mu3.k[i], sd = sigma.k[i])
  traits[(N*(i-1)+1):(N*i), 4] <- rnorm(n = N, mean = mu4.k[i], sd = sigma.k[i])
}
rm(i)
strains <- rep(x = 1:K, each = N)




gt.data <- getStanData(genotype = genotype,
                       traits = traits,
                       trait.type = c("Q", "Q", "Q", "Q"),
                       strains = strains)





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


om0 <- rstan::stan_model(file = "src/stan_files/M2_loglik.stan")
om0c <- rstan::stan_model(file = "src/stan_files/M2c_loglik.stan")
m0 <- rstan::stan_model(file = "src/stan_files/noise/tM2_loglik.stan")
m0c <- rstan::stan_model(file = "src/stan_files/noise/tM2c_loglik.stan")



# run
op0 <- rstan::sampling(object = om0,
                      data = data.list,
                      iter = 1500,
                      warmup = 500,
                      chains = 4,
                      cores = 4,
                      control = list(adapt_delta = 0.97,
                                     max_treedepth = 10))

# run
op0c <- rstan::sampling(object = om0c,
                       data = data.list,
                       iter = 1500,
                       warmup = 500,
                       chains = 4,
                       cores = 4,
                       control = list(adapt_delta = 0.97,
                                      max_treedepth = 10))

# run
p0 <- rstan::sampling(object = m0,
                      data = data.list,
                      iter = 1500,
                      warmup = 500,
                      chains = 4,
                      cores = 4,
                      control = list(adapt_delta = 0.97,
                                     max_treedepth = 10))

# run
p0c <- rstan::sampling(object = m0c,
                      data = data.list,
                      iter = 1500,
                      warmup = 500,
                      chains = 4,
                      cores = 4,
                      control = list(adapt_delta = 0.95,
                                     max_treedepth = 10))




loo::compare(loo::loo(loo::extract_log_lik(stanfit = op0)),
             loo::loo(loo::extract_log_lik(stanfit = op0c)),
             loo::loo(loo::extract_log_lik(stanfit = p0)),
             loo::loo(loo::extract_log_lik(stanfit = p0c)))

x <- data.frame(summary(p0, pars = "mu_beta")$summary)
x$par <- rownames(x)
x[x$par == "mu_beta[1,1]", ]
y <- regx$par

y <- data.frame(summary(op0, pars = "mu_beta")$summary)
y$par <- rownames(y)
y[y$par == "mu_beta[1,1]", ]

plot(x = x$mean, y  = y$mean)
abline(0, 1)
