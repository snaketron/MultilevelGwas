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



genotype <- cbind(c(sample(x = c("A", "C"), size = 40, prob = c(0, 1), replace = T),
                    sample(x = c("A", "C"), size = 40, prob = c(1, 0), replace = T)),
                  c(sample(x = c("H", "C"), size = 40, prob = c(0.15, 0.85), replace = T),
                    sample(x = c("H", "C"), size = 40, prob = c(0.85, 0.15), replace = T)),
                  c(sample(x = c("J", "G"), size = 40, prob = c(0.15, 0.85), replace = T),
                    sample(x = c("J", "G"), size = 40, prob = c(0.85, 0.15), replace = T)),
                  replicate(n = 22, expr = sample(x = sample(x = LETTERS, size = 2, replace = F),
                                                  size = 80, replace = T)))
trait.1 <- c(rnorm(n = 20, mean = 0, sd = 0.3),
             rnorm(n = 20, mean = 0.25, sd = 0.5),
             rnorm(n = 20, mean = 0.75, sd = 0.5),
             rnorm(n = 20, mean = 1, sd = 0.5))
trait.2 <- c(rbinom(n = 20, size = 1, prob = 0.3),
             rbinom(n = 20, size = 1, prob = 0.4),
             rbinom(n = 20, size = 1, prob = 0.6),
             rbinom(n = 20, size = 1, prob = 0.7))
traits <- cbind(trait.1, trait.2)
strains <- rep(x = c(1, 2, 3, 4), each = 20)
rm(trait.1, trait.2)

plot(x = strains, y = traits[, 1])
plot(x = strains, y = traits[, 1], col = as.factor(genotype[, 1]))



# out <- runMLgwas(genotype = genotype,
#                  traits = traits,
#                  trait.type = c("Q", "D"),
#                  strains = strains,
#                  model = "M0",
#                  mcmc.chains = 4,
#                  mcmc.steps = 2500,
#                  mcmc.warmup = 500,
#                  cores = 8,
#                  hdi.level = 0.95,
#                  stat.learn.method = "svm",
#                  cv.steps = 100)



mc <- runModelComparison(genotype = genotype,
                         traits = traits,
                         trait.type = c("Q", "D"),
                         strains = strains,
                         models = c("M0", "M1"),#c("M0", "M1"),
                         mcmc.chains = 4,
                         mcmc.steps = 1500,
                         mcmc.warmup = 500,
                         cores = 4,
                         hdi.level = 0.95)
