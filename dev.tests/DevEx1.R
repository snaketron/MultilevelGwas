# Example 1:
# Data:
# SNPs : 30
# Traits: 1
# Strains: 14
# N per strain: 8
# Homozygous within strain
# Remark: 5 strong effects, 25 noise


require(rstan)
require(loo)
require(parallel)
require(foreach)
require(doParallel)
require(Biostrings)


source("R/main.R")
source("R/util.R")
source("R/bayesian.R")
source("R/statlearn.R")
source("R/modelcomparison.R")
source("R/ppc.R")



# # simulate data
# N <- 8
# K <- 14
# P <- 1
#
# genotype <- cbind(rep(x = c("A", "C"), each = N*K/2),
#                   rep(x = c("C", "A", "A", "C"), times = c(N, N*(K/2-1), N, N*(K/2-1))),
#                   rep(x = c("C", "A", "A", "C"), times = c(N*2, N*(K/2-2), N*2, N*(K/2-2))),
#                   rep(x = c("C", "A", "A", "C"), times = c(N*3, N*(K/2-3), N*3, N*(K/2-3))),
#                   rep(x = c("C", "A", "A", "C"), times = c(N*4, N*(K/2-4), N*4, N*(K/2-4))),
#                   replicate(n = 25, rep(x = sample(x = c("A", "C"),
#                             size = K, replace = T), each = N)))
#
#
# mu1.k <- c(rnorm(n = K/2, mean = 2, sd = 1),
#           rnorm(n = K/2, mean = -2, sd = 1))
# sigma.k <- abs(rnorm(n = K, mean = 0, sd = 1))
#
# traits <- matrix(data = 0, nrow = N*K, ncol = P)
# for(i in 1:K) {
#   traits[(N*(i-1)+1):(N*i), 1] <- rnorm(n = N, mean = mu1.k[i], sd = sigma.k[i])
# }
# rm(i)
# strains <- rep(x = 1:K, each = N)
#
#
#
# d <- rbind(data.frame(Y = traits[, 1], K = strains, trait = "1"))
#
#
# plot(traits[, 1], x = strains, col = as.factor(genotype[, 1]))
#
# d <- list(genotype = genotype,
#           traits = traits,
#           strains = strains)
# save(d, file = "dev.tests/ex1/d.RData")





d <- get(load(file = "dev.tests/ex1/d.RData"))

mc.ex1 <- runModelComparison(genotype = d$genotype,
                          traits = d$traits,
                          trait.type = c("Q"),
                          strains = d$strains,
                          models = c("M0", "M1"),
                          mcmc.chains = 4,
                          mcmc.steps = 1500,
                          mcmc.warmup = 500,
                          cores = 4,
                          hdi.level = 0.95,
                          adapt_delta = 0.95,
                          max_treedepth = 10)

save(mc.ex1, file = "dev.tests/ex1/mc.ex1.RData")



# posteriors
summary(mc.ex1$ps$M0, par = c("nu", "nu_help", "mu_beta"))$summary
summary(mc.ex1$ps$M1, par = c("nu", "nu_help", "grand_mu_beta"))$summary


getTMoments(nu = summary(mc.ex1$ps$M0, par = c("nu"))$summary[1, "mean"],
            mu = summary(mc.ex1$ps$M0, par = c("mu_beta"))$summary[1, "mean"])
getTMoments(nu = summary(mc.ex1$ps$M1, par = c("nu"))$summary[1, "mean"],
            mu = summary(mc.ex1$ps$M1, par = c("grand_mu_beta"))$summary[1, "mean"])


# correlation between snp slopes
plot(summary(mc.ex1$ps$M0, par = c("beta"))$summary[, "mean"],
     summary(mc.ex1$ps$M1, par = c("mu_beta"))$summary[, "mean"])
abline(0, 1)



# check loo.ic
mc.ex1$ic$loo
mc.ex1$ic$waic




# ppc
ppc <- do.call(rbind, mc.ex1$ppc)

ggplot(data = ppc[ppc$prediction.level == "individual-level", ])+
  facet_grid(facets = t~model)+
  geom_density2d(aes(y = y.ppc.mean, x = y.real.mean), col = "orange")+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()

ggplot(data = ppc[ppc$prediction.level == "strain-level", ])+
  facet_grid(facets = t~model)+
  geom_density2d(aes(y = y.ppc.mean, x = y.real.mean), col = "orange")+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()

ggplot(data = ppc[ppc$prediction.level == "snp-level" &
                    ppc$parameter.level == "lowest", ])+
  facet_grid(facets = t~model)+
  geom_density2d(aes(y = y.ppc.mean, x = y.real.mean), col = "orange")+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()

ggplot(data = ppc[ppc$prediction.level == "snp-level", ])+
  facet_grid(facets = t~model+parameter.level)+
  geom_errorbar(aes(ymin = y.ppc.L, ymax = y.ppc.H, x = y.real.mean), col = "darkgray")+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(aes(y = y.ppc.mean, x = y.real.mean), shape = 21, fill = "white")+
  theme_bw()
