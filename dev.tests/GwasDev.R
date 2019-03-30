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



data.list <- get(load(file = "~/MiceGwas/output/posterior/vsv/data.list.RData"))
genotype <- data.list$X[, 1:300]
genotype[1, 1] <- as.character(genotype[1, 1])
traits <- cbind(data.list$Yc, data.list$Yd)

strains <- data.list$K



mc <- runModelComparison(genotype = genotype,
                         traits = traits,
                         trait.type = c("Q", "Q"),
                         # trait.type = c("Q", "Q", "Q", "D"),
                         strains = strains,
                         models = c("M0c", "M1c", "M2c"),
                         mcmc.chains = 4,
                         mcmc.steps = 2500,
                         mcmc.warmup = 1000,
                         cores = 4,
                         hdi.level = 0.95,
                         adapt_delta = 0.99)

save(mc, file = "dev.tests/mc.vsv.cov.300.RData")


