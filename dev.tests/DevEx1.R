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


rho.12 <- rnorm(n = K, mean = 0.75, sd = 0.50)
rho.13 <- rnorm(n = K, mean = 0.66, sd = 0.50)
rho.14 <- rnorm(n = K, mean = 0.85, sd = 0.50)
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



d <- rbind(data.frame(Y = traits[, 1], K = strains, trait = "1"),
           data.frame(Y = traits[, 2], K = strains, trait = "2"),
           data.frame(Y = traits[, 3], K = strains, trait = "3"),
           data.frame(Y = traits[, 4], K = strains, trait = "4"))

ggplot(data = d)+
  facet_wrap(facets = ~trait, ncol = 1)+
  geom_jitter(aes(x = K, y = Y),
              width = 0.1, height = 0)+
  theme_bw()



# create data list
data.list <- list(Yq = traits,
                  Yd = matrix(data = 0,
                              nrow = nrow(traits),
                              ncol = 0),
                  K = strains,
                  Nk = max(strains),
                  Ntd = 0,
                  Ntq = ncol(traits),
                  N = nrow(traits))






mc <- runModelComparison(genotype = genotype,
                         traits = traits,
                         trait.type = c("Q", "Q", "Q", "Q"),
                         strains = strains,
                         models = c("M0", "M0c", "M1", "M1c", "M2", "M2c"),
                         mcmc.chains = 4,
                         mcmc.steps = 3000,
                         mcmc.warmup = 1500,
                         cores = 4,
                         hdi.level = 0.95,
                         adapt_delta = 0.99,
                         max_treedepth = 12)

save(mc, file = "dev.tests/dev.ex1.RData")


# check loo.ic
mc$ic$loo
mc$ic$waic



mc$ppc$M0$model <- "M0"
mc$ppc$M1$model <- "M1"
mc$ppc$M2$model <- "M2"
mc$ppc$M0c$model <- "M0c"
mc$ppc$M1c$model <- "M1c"
mc$ppc$M2c$model <- "M2c"
ppc <- rbind(mc$ppc$M0, mc$ppc$M1, mc$ppc$M2,
             mc$ppc$M0c, mc$ppc$M1c, mc$ppc$M2c)

ggplot(ppc[ppc$level == "snp", ])+
  facet_grid(facets = t~model)+
  geom_errorbar(aes(x = s, ymin = error.L, ymax = error.H),
                col = "black")+
  geom_point(aes(x = s, y = error), shape = 21,
             fill = "white", col = "black")+
  theme_bw()

ggplot(ppc[ppc$level == "strain" & ppc$s == 1, ])+
  facet_grid(facets = t~model)+
  geom_errorbar(aes(x = k, y = error, ymin = error.L,
                    ymax = error.H), col = "black")+
  geom_point(aes(x = k, y = error), shape = 21,
             fill = "white", col = "black")+
  theme_bw()



x <- data.frame(extract(object = mc$ps$M2,
                pars = c("alpha", "beta", "mu_beta")))
hist(x$mu_beta.1.1, breaks = 100)
hist(x$mu_beta.4.1, breaks = 100)
hist(x$mu_beta.1.25, breaks = 100)
