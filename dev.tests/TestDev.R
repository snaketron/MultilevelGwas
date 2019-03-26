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

genotype <- cbind(rep(x = c("A", "C"), each = N*10/2),
                  replicate(n = 49, rep(x = sample(x = c("A", "C"),
                            size = 10, replace = T), each = N)))


mu.k <- c(rnorm(n = 5, mean = 2, sd = 1),
          rnorm(n = 5, mean = 0, sd = 1))

sigmas.k <- c(0.3, 0.6, 1.5, 0.4, 0.5,
              0.3, 0.4, 0.5, 0.4, 1.2)

trait.1 <- c(rnorm(n = N, mean = mu.k[1], sd = sigmas.k[1]),
             rnorm(n = N, mean = mu.k[2], sd = sigmas.k[2]),
             rnorm(n = N, mean = mu.k[3], sd = sigmas.k[3]),
             rnorm(n = N, mean = mu.k[4], sd = sigmas.k[4]),
             rnorm(n = N, mean = mu.k[5], sd = sigmas.k[5]), #5
             rnorm(n = N, mean = mu.k[6], sd = sigmas.k[6]),
             rnorm(n = N, mean = mu.k[7], sd = sigmas.k[7]),
             rnorm(n = N, mean = mu.k[8], sd = sigmas.k[8]),
             rnorm(n = N, mean = mu.k[9], sd = sigmas.k[9]),
             rnorm(n = N, mean = mu.k[10], sd = sigmas.k[10]))

trait.2 <- c(rnorm(n = N, mean = mu.k[1], sd = sigmas.k[1]),
             rnorm(n = N, mean = mu.k[2], sd = sigmas.k[2]),
             rnorm(n = N, mean = mu.k[3], sd = sigmas.k[3]),
             rnorm(n = N, mean = mu.k[4], sd = sigmas.k[4]),
             rnorm(n = N, mean = mu.k[5], sd = sigmas.k[5]), #5
             rnorm(n = N, mean = mu.k[6], sd = sigmas.k[6]),
             rnorm(n = N, mean = mu.k[7], sd = sigmas.k[7]),
             rnorm(n = N, mean = mu.k[8], sd = sigmas.k[8]),
             rnorm(n = N, mean = mu.k[9], sd = sigmas.k[9]),
             rnorm(n = N, mean = mu.k[10], sd = sigmas.k[10]))

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

par(mfrow = c(3, 1))
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



s <- data.frame(summary(mc$ps$M2, par = "mu_beta")$summary)
s$par <- rownames(s)
s <- s[regexpr(pattern = "mu_beta\\[3,", text = s$par) != -1, ]
plot(density(s$mean))
points(x = s$mean, y = rep(x = 0, times = length(s$mean)))
hist(s$mean, breaks = 100)
