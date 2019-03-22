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
N <- 7

genotype <- cbind(rep(x = c("A", "C"), each = N*10/2),
                  replicate(n = 19, rep(x = sample(x = c("A", "C"),
                            size = 10, replace = T), each = N)))


sigmas.k <- c(0.3, 0.4, 0.5, 0.4, 0.5,
              0.3, 0.4, 0.5, 0.4, 0.8)

trait.1 <- c(rnorm(n = N, mean = -1.5, sd = sigmas.k[1]),
             rnorm(n = N, mean = -2, sd = sigmas.k[2]),
             rnorm(n = N, mean = -1, sd = sigmas.k[3]),
             rnorm(n = N, mean = 0.2, sd = sigmas.k[4]),
             rnorm(n = N, mean = 0.25, sd = sigmas.k[5]), #5
             rnorm(n = N, mean = 1.1, sd = sigmas.k[6]),
             rnorm(n = N, mean = 1.25, sd = sigmas.k[7]),
             rnorm(n = N, mean = 1.5, sd = sigmas.k[8]),
             rnorm(n = N, mean = 2, sd = sigmas.k[9]),
             rnorm(n = N, mean = 2, sd = sigmas.k[10]))

trait.2 <- c(rnorm(n = N, mean = -1.5, sd = sigmas.k[1]),
             rnorm(n = N, mean = -1, sd = sigmas.k[2]),
             rnorm(n = N, mean = -0.75, sd = sigmas.k[3]),
             rnorm(n = N, mean = 0.2, sd = sigmas.k[4]),
             rnorm(n = N, mean = 0.25, sd = sigmas.k[5]), #5
             rnorm(n = N, mean = 1.1, sd = sigmas.k[6]),
             rnorm(n = N, mean = 1.25, sd = sigmas.k[7]),
             rnorm(n = N, mean = 1.2, sd = sigmas.k[8]),
             rnorm(n = N, mean = 2, sd = sigmas.k[9]),
             rnorm(n = N, mean = 2, sd = sigmas.k[10]))

trait.3 <- c(rbinom(n = N, size = 1, prob = 0.05),
             rbinom(n = N, size = 1, prob = 0.2),
             rbinom(n = N, size = 1, prob = 0.3),
             rbinom(n = N, size = 1, prob = 0.2),
             rbinom(n = N, size = 1, prob = 0.3), #5
             rbinom(n = N, size = 1, prob = 0.65),
             rbinom(n = N, size = 1, prob = 0.7),
             rbinom(n = N, size = 1, prob = 0.75),
             rbinom(n = N, size = 1, prob = 0.75),
             rbinom(n = N, size = 1, prob = 0.8))

traits <- cbind(trait.1, trait.2, trait.3)
strains <- rep(x = 1:10, each = N)
rm(trait.1, trait.2, trait.3)

plot(x = strains, y = traits[, 1],
     col = as.factor(genotype[, 1]))
plot(x = strains, y = traits[, 2],
     col = as.factor(genotype[, 1]))
plot(x = strains, y = traits[, 3],
     col = as.factor(genotype[, 1]))




# rc <- runMLgwas(genotype = genotype,
#                 traits = traits,
#                 trait.type = c("Q", "Q", "D"),
#                 strains = strains,
#                 model = "M2",
#                 mcmc.chains = 4,
#                 mcmc.steps = 1500,
#                 mcmc.warmup = 500,
#                 cores = 4,
#                 hdi.level = 0.95,
#                 stat.learn.method = "svm",
#                 cv.steps = 100)




mc <- runModelComparison(genotype = genotype,
                         traits = traits,
                         trait.type = c("Q", "Q", "D"),
                         strains = strains,
                         models = c("M0", "M1", "M2"),
                         mcmc.chains = 4,
                         mcmc.steps = 2500,
                         mcmc.warmup = 1000,
                         cores = 4,
                         hdi.level = 0.999)



