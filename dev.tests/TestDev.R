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
genotype <- cbind(rep(x = c("A", "C"), each = 40),
                  replicate(n = 19, rep(x = sample(x = c("A", "C"),
                            size = 16, replace = T), each = 5)))

# mixture within strains
# genotype <- cbind(c(sample(x = c("A", "C"), size = 40, prob = c(0, 1), replace = T),
#                     sample(x = c("A", "C"), size = 40, prob = c(1, 0), replace = T)),
#                   c(sample(x = c("H", "C"), size = 40, prob = c(0.15, 0.85), replace = T),
#                     sample(x = c("H", "C"), size = 40, prob = c(0.85, 0.15), replace = T)),
#                   c(sample(x = c("J", "G"), size = 40, prob = c(0.15, 0.85), replace = T),
#                     sample(x = c("J", "G"), size = 40, prob = c(0.85, 0.15), replace = T)),
#                   replicate(n = 12, expr = sample(x = sample(x = LETTERS,
#                             size = 2, replace = F), size = 80, replace = T)))

trait.1 <- c(rnorm(n = 5, mean = -2.5, sd = 0.5),
             rnorm(n = 5, mean = -2, sd = 0.5),
             rnorm(n = 5, mean = -1, sd = 0.4),
             rnorm(n = 5, mean = 0, sd = 0.5),
             rnorm(n = 5, mean = 0, sd = 0.2),
             rnorm(n = 5, mean = 0.1, sd = 0.3),
             rnorm(n = 5, mean = 0.2, sd = 0.3),
             rnorm(n = 5, mean = 0.25, sd = 0.2), #8
             rnorm(n = 5, mean = 0.9, sd = 0.1),
             rnorm(n = 5, mean = 1, sd = 0.2),
             rnorm(n = 5, mean = 1.1, sd = 0.2),
             rnorm(n = 5, mean = 1.25, sd = 0.3),
             rnorm(n = 5, mean = 0.9, sd = 0.2),
             rnorm(n = 5, mean = 2, sd = 0.3),
             rnorm(n = 5, mean = 2, sd = 0.4),
             rnorm(n = 5, mean = 3, sd = 0.5))

trait.2 <- c(rnorm(n = 5, mean = -1, sd = 0.3),
             rnorm(n = 5, mean = -2, sd = 0.2),
             rnorm(n = 5, mean = -0.5, sd = 0.2),
             rnorm(n = 5, mean = 0, sd = 0.2),
             rnorm(n = 5, mean = 0, sd = 0.1),
             rnorm(n = 5, mean = 0.1, sd = 0.2),
             rnorm(n = 5, mean = 0.2, sd = 0.1),
             rnorm(n = 5, mean = 0.25, sd = 0.2), #8
             rnorm(n = 5, mean = 0.9, sd = 0.1),
             rnorm(n = 5, mean = 1, sd = 0.2),
             rnorm(n = 5, mean = 1.1, sd = 0.2),
             rnorm(n = 5, mean = 1.25, sd = 0.2),
             rnorm(n = 5, mean = 0.9, sd = 0.3),
             rnorm(n = 5, mean = 1, sd = 0.2),
             rnorm(n = 5, mean = 2, sd = 0.5),
             rnorm(n = 5, mean = 3, sd = 0.5))

trait.3 <- c(rbinom(n = 5, size = 1, prob = 0.01),
             rbinom(n = 5, size = 1, prob = 0.01),
             rbinom(n = 5, size = 1, prob = 0.2),
             rbinom(n = 5, size = 1, prob = 0.3),
             rbinom(n = 5, size = 1, prob = 0.2),
             rbinom(n = 5, size = 1, prob = 0.3),
             rbinom(n = 5, size = 1, prob = 0.2),
             rbinom(n = 5, size = 1, prob = 0.3),
             rbinom(n = 5, size = 1, prob = 0.65),
             rbinom(n = 5, size = 1, prob = 0.7),
             rbinom(n = 5, size = 1, prob = 0.75),
             rbinom(n = 5, size = 1, prob = 0.8),
             rbinom(n = 5, size = 1, prob = 0.65),
             rbinom(n = 5, size = 1, prob = 0.7),
             rbinom(n = 5, size = 1, prob = 0.75),
             rbinom(n = 5, size = 1, prob = 0.8))

traits <- cbind(trait.1, trait.2, trait.3)
strains <- rep(x = 1:16, each = 5)
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
                         mcmc.steps = 1500,
                         mcmc.warmup = 500,
                         cores = 4,
                         hdi.level = 0.95)




