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



# homozygous strains
N <- 5
K <- 10

genotype <- cbind(rep(x = c("A", "C"), each = N*K/2),
                  replicate(n = 24, rep(x = sample(x = c("A", "C"),
                            size = K, replace = T), each = N)))
unique(apply(X = genotype, MARGIN = 2, FUN = base::paste, collapse = ''))


rho.12 <- rnorm(n = K, mean = 1, sd = 0.1)
rho.13 <- rnorm(n = K, mean = 1, sd = 0.1)
rho.14 <- rnorm(n = K, mean = 1, sd = 0.1)
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




mc <- runModelComparison(genotype = genotype,
                         traits = traits,
                         trait.type = c("Q", "Q", "Q", "Q"),
                         strains = strains,
                         models = c("M0", "M1"),
                         mcmc.chains = 4,
                         mcmc.steps = 2000,
                         mcmc.warmup = 500,
                         cores = 4,
                         hdi.level = 0.95,
                         adapt_delta = 0.95,
                         max_treedepth = 12)

save(mc, file = "dev.tests/dev.ex2.RData")




# check loo.ic
mc$ic$loo
mc$ic$waic

mc$ppc <- x



ppc <- rbind(mc$ppc$M0, mc$ppc$M1, mc$ppc$M2,
             mc$ppc$M0c, mc$ppc$M1c, mc$ppc$M2c)


#### Parameters: lowest ####
# Prediction: individual-level
ggplot(data = ppc[ppc$prediction.level == "individual-level", ])+
  facet_grid(facets = t~model)+
  geom_errorbar(aes(ymin = y.ppc.L, ymax = y.ppc.H, x = y.real.mean), col = "darkgray")+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(aes(y = y.ppc.mean, x = y.real.mean), shape = 21, fill = "white")+
  theme_bw()

ggplot(data = ppc[ppc$prediction.level == "individual-level", ])+
  facet_grid(facets = t~model)+
  geom_density2d(aes(y = y.ppc.mean, x = y.real.mean), col = "orange")+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()



# Prediction: strain-level
ggplot(data = ppc[ppc$prediction.level == "strain-level", ])+
  facet_grid(facets = t~model)+
  geom_errorbar(aes(ymin = y.ppc.L, ymax = y.ppc.H, x = y.real.mean), col = "darkgray")+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(aes(y = y.ppc.mean, x = y.real.mean), shape = 21, fill = "white")+
  theme_bw()

ggplot(data = ppc[ppc$prediction.level == "strain-level", ])+
  facet_grid(facets = t~model)+
  geom_density2d(aes(y = y.ppc.mean, x = y.real.mean), col = "orange")+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()



# Prediction: snp-level
ggplot(data = ppc[ppc$prediction.level == "snp-level" &
                    ppc$parameters.level == "lowest", ])+
  facet_grid(facets = t~model)+
  geom_errorbar(aes(ymin = y.ppc.L, ymax = y.ppc.H, x = y.real.mean), col = "darkgray")+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(aes(y = y.ppc.mean, x = y.real.mean), shape = 21, fill = "white")+
  theme_bw()

ggplot(data = ppc[ppc$prediction.level == "snp-level" &
                    ppc$parameters.level == "lowest", ])+
  facet_grid(facets = t~model)+
  geom_density2d(aes(y = y.ppc.mean, x = y.real.mean), col = "orange")+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()



#### Parameters: lowest-plus-one ####

ggplot(data = ppc[ppc$prediction.level == "snp-level" &
                    ppc$parameters.level == "lowest-plus-one", ])+
  facet_grid(facets = t~model)+
  geom_errorbar(aes(ymin = y.ppc.L, ymax = y.ppc.H, x = y.real.mean), col = "darkgray")+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(aes(y = y.ppc.mean, x = y.real.mean), shape = 21, fill = "white")+
  theme_bw()



ggplot(data = ppc[ppc$prediction.level == "snp-level", ])+
  facet_grid(facets = t~model+parameters.level)+
  geom_errorbar(aes(ymin = y.ppc.L, ymax = y.ppc.H, x = y.real.mean), col = "darkgray")+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(aes(y = y.ppc.mean, x = y.real.mean), shape = 21, fill = "white")+
  theme_bw()




# Function:
# summarize errors
getErrorSummary <- function(ppc, hdi.level) {
  source("R/util.R")
  error.summary <- c()

  key <- ppc[, c("model", "t", "level")]
  key <- key[duplicated(key) == F, ]

  for(i in 1:nrow(key)) {
    error <- ppc$error.mean[ppc$level == key$level[i] &
                              ppc$t == key$t[i] &
                              ppc$model == key$model[i]]
    error <- error[is.finite(error)]

    # error.summary
    error.hdi <- getHdi(vec = error, hdi.level = hdi.level)
    row <- data.frame(model = key$model[i],
                      t = key$t[i],
                      level = key$level[i],
                      error.mean = mean(error),
                      error.L = error.hdi[1],
                      error.H = error.hdi[2])
    error.summary <- rbind(error.summary, row)

    # rm
    rm(row, error.hdi, error)
  }


  error.summary$model <- factor(x = error.summary$model,
                                levels = c("M0", "M0c", "M1",
                                           "M1c", "M2", "M2c"))
  return (error.summary)
}



error.summary <- getErrorSummary(ppc = ppc, hdi.level = 0.95)
error.summary$cov <- F
error.summary$cov[regexpr(pattern = "c", text = error.summary$model) != -1] <- T
error.summary$model <- gsub(pattern = 'c', replacement = '', x = error.summary$model)

ggplot(data = error.summary)+
  facet_grid(facets = t~level)+
  geom_errorbar(aes(x = model, ymin = error.L, ymax = error.H, shape = cov),
                col = "black", width = 0.35, position = position_dodge(width = 0.5))+
  geom_point(aes(x = model, y = error.mean, fill = model, shape = cov, col = model),
             size = 3, position = position_dodge(width = 0.5))+
  theme_bw()+
  theme(legend.position = "top")



x <- summary(mc$ps$M0, pars = "beta")$summary
y <- summary(mc$ps$M1c, pars = "mu_beta")$summary
plot(x = x[301:600, "mean"], y = y[301:600, "mean"])
abline(0, 1)

z <- data.frame(summary(mc$ps$M1c)$summary)
z$par <- rownames(z)
