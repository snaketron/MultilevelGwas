# Example 1:
# Data:
# SNPs : 30
# Traits: 3 (correlated)
# Strains: 14
# N per strain: 8
# Homozygous within strain
# Remark: 3 strong effects, 27 noise


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
#
# genotype <- cbind(rep(x = c("A", "C"), each = N*K/2),
#                   rep(x = c("C", "A", "A", "C"), times = c(N, N*(K/2-1), N, N*(K/2-1))),
#                   rep(x = c("C", "A", "A", "C"), times = c(N*2, N*(K/2-2), N*2, N*(K/2-2))),
#                   replicate(n = 27, rep(x = sample(x = c("A", "C"),
#                             size = K, replace = T), each = N)))
# unique(apply(X = genotype, MARGIN = 2, FUN = base::paste, collapse = ''))
#
#
# rho.12 <- rnorm(n = K, mean = 0.8, sd = 0.1)
# rho.13 <- rnorm(n = K, mean = 0.8, sd = 0.1)
# rho.12 <- sapply(X = rho.12, FUN = getMax <- function(x) {min(x, 1)})
# rho.13 <- sapply(X = rho.13, FUN = getMax <- function(x) {min(x, 1)})
#
#
# mu1.k <- c(rnorm(n = K/2, mean = 2.5, sd = 1),
#           rnorm(n = K/2, mean = 0, sd = 1))
# mu2.k <- mu1.k * rho.12
# mu3.k <- mu1.k * rho.13
#
# library("Hmisc")
# library("corrplot")
# corrplot(rcorr(cbind(mu1.k, mu2.k, mu3.k))$r,
#          type = "upper", order = "hclust",
#          tl.col = "black", tl.srt = 45)
#
# sigma.k <- abs(rnorm(n = K, mean = 0.5, sd = 1))
#
# traits <- matrix(data = 0, nrow = N*K, ncol = 3)
# for(i in 1:K) {
#   traits[(N*(i-1)+1):(N*i), 1] <- rnorm(n = N, mean = mu1.k[i], sd = sigma.k[i])
#   traits[(N*(i-1)+1):(N*i), 2] <- rnorm(n = N, mean = mu2.k[i], sd = sigma.k[i])
#   traits[(N*(i-1)+1):(N*i), 3] <- rnorm(n = N, mean = mu3.k[i], sd = sigma.k[i])
# }
# rm(i)
# strains <- rep(x = 1:K, each = N)
#
#
#
# d <- rbind(data.frame(Y = traits[, 1], K = strains, trait = "1"),
#            data.frame(Y = traits[, 2], K = strains, trait = "2"),
#            data.frame(Y = traits[, 3], K = strains, trait = "3"))
#
# d <- list(genotype = genotype,
#           traits = traits,
#           strains = strains)
# save(d, file = "dev.tests/ex3/d.RData")


d <- get(load(file = "dev.tests/ex3/d.RData"))

mc.ex3 <- runModelComparison(genotype = d$genotype,
                             traits = d$traits,
                             trait.type = c("Q", "Q", "Q"),
                             strains = d$strains,
                             models = c("M0", "M0c",
                                        "M1", "M1c"),
                             mcmc.chains = 4,
                             mcmc.steps = 1500,
                             mcmc.warmup = 500,
                             cores = 4,
                             hdi.level = 0.95,
                             adapt_delta = 0.95,
                             max_treedepth = 10)




x <- data.frame(rstan::extract(object = mc.ex3$ps$M0, pars = "log_lik2"))
for(t in 1:3) {
  print(loo::loo(x = as.matrix(x[, which(regexpr(pattern = paste("log_lik2\\.", t, sep = ''),
                                                 text = colnames(x)) != -1)])))
}
x.loo <- loo::extract_log_lik(stanfit = mc.ex3$ps$M0, parameter_name = "log_lik2")
loo::loo(x = x.loo)

m1 <- as.matrix(x)
which(regexpr(pattern = "log_lik2\\.1", text = colnames(x)) != -1)




save(mc.ex3, file = "dev.tests/ex3/mc.ex3.RData")



# posterior
mc.ex3 <- get(load(file = "dev.tests/ex3/mc.ex3.RData"))
s <- data.frame(summary(mc.ex3$ps$M1)$summary)
s$par <- rownames(s)



s <- summary(mc.ex3$ps$M1, par = c("nu", "nu_help", "grand_mu_beta", "mu_beta"))$summary
s["nu", "mean"]

getTMoments(nu = s["nu", "mean"], mu = s["grand_mu_beta[1]", "mean"])
plot(density(summary(mc.ex3$ps$M1, par = "mu_beta")$summary[, "mean"]))
s["nu", "mean"]
s["nu_help", "mean"]




# check loo.ic
mc.ex3$ic$loo
mc.ex3$ic$waic
loo::compare(mc.ex3$ic$loo$M0, mc.ex3$ic$loo$M0c)



# ppc
ppc <- do.call(rbind, mc.ex3$ppc)

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
