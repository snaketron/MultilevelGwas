source("R/main.R")
source("R/util.R")
source("R/bayesian.R")
source("R/statlearn.R")
require(rstan)
require(loo)
require(parallel)
require(foreach)
require(doParallel)
require(Biostrings)



set.seed(seed = 12)
genotype <- cbind(rep(x = c("A", "C"), each = 15),
                  replicate(n = 99, expr = sample(x = c("G", "T"),
                                                  size = 30, replace = T)))
trait.1 <- c(rnorm(n = 15, mean = 0, sd = 1),
                 rnorm(n = 15, mean = 0.5, sd = 1))
trait.2 <- c(rbinom(n = 15, size = 1, prob = 0.3),
                 rbinom(n = 15, size = 1, prob = 0.5))
traits <- cbind(trait.1, trait.2)
strains <- rep(x = c(1, 2, 3), each = 10)
rm(trait.1, trait.2)


stan.model <- rstan::stan_model(file = "src/stan_files/M0.stan")


out <- runMLgwas(genotype = genotype,
                 traits = traits,
                 trait.type = c("Q", "D"),
                 strains = strains,
                 model = "M0",
                 mcmc.chains = 4,
                 mcmc.steps = 2500,
                 mcmc.warmup = 500,
                 cores = 8,
                 hdi.level = 0.95,
                 stat.learn.method = "svm",
                 cv.steps = 100)
