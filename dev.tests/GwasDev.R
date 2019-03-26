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



data.list <- get(load(file = "~/MiceGwas/output/posterior/tnf/data.list.RData"))
genotype <- data.list$X[, 1:100]
genotype[1, 1] <- as.character(genotype[1, 1])
traits <- cbind(data.list$Yc, data.list$Yd)

strains <- data.list$K



mc <- runModelComparison(genotype = genotype,
                         traits = traits,
                         trait.type = c("Q", "Q", "Q", "D"),
                         strains = strains,
                         models = c("M0", "M1", "M2"),
                         mcmc.chains = 4,
                         mcmc.steps = 1500,
                         mcmc.warmup = 500,
                         cores = 4,
                         hdi.level = 0.95,
                         adapt_delta = 0.99)




s <- data.frame(summary(mc$ps$M2, par = "mu_beta")$summary)
s$par <- rownames(s)
s <- s[regexpr(pattern = "mu_beta\\[3,", text = s$par) != -1, ]
plot(density(s$mean))
points(x = s$mean, y = rep(x = 0, times = length(s$mean)))
hist(s$mean, breaks = 100)
save(mc, file = "dev.tests/mc.vsv.RData")
