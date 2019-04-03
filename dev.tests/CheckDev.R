require(rstan)
require(loo)
require(foreach)
require(doParallel)
require(parallel)
require(ggplot2)

source("R/main.R")
source("R/util.R")
source("R/bayesian.R")
source("R/statlearn.R")
source("R/modelcomparison.R")
source("R/ppc.R")

# mc.vsv <- get(load("~/MLgwas/dev.tests/vsv.300.RData"))
mc.lcmv <- get(load("~/MLgwas/dev.tests/lcmv.300.RData"))
# mc.tnf.1 <- get(load("~/MLgwas/dev.tests/tnf.300.RData"))
# mc.tnf.2 <- get(load("~/MLgwas/dev.tests/tnf.cov.300.RData"))
rm(mc)

x <- getPpcMc(ps = mc.tnf.1$ps,
              gt.data = mc.tnf.1$gt.data,
              models = names(mc.tnf.1$ps),
              hdi.level = 0.95, cores = 6)
mc.tnf.1$ppc <- x
save(mc.tnf.1, file = "~/MLgwas/dev.tests/tnf.300.RData")


x <- getPpcMc(ps = mc.tnf.2$ps,
              gt.data = mc.tnf.2$gt.data,
              models = names(mc.tnf.2$ps),
              hdi.level = 0.95, cores = 6)
mc.tnf.2$ppc <- x
save(mc.tnf.2, file = "~/MLgwas/dev.tests/tnf.cov.300.RData")



x <- getPpcMc(ps = mc.vsv$ps,
              gt.data = mc.vsv$gt.data,
              models = names(mc.vsv$ps),
              hdi.level = 0.95, cores = 6)
mc.vsv$ppc <- x
save(mc.vsv, file = "~/MLgwas/dev.tests/vsv.300.RData")



x <- getPpcMc(ps = mc.lcmv$ps,
              gt.data = mc.lcmv$gt.data,
              models = names(mc.lcmv$ps),
              hdi.level = 0.95, cores = 6)
mc.lcmv$ppc <- x
save(mc.lcmv, file = "~/MLgwas/dev.tests/lcmv.300.RData")






# vsv
ppc <- rbind(mc.vsv$ppc$M0, mc.vsv$ppc$M1, mc.vsv$ppc$M2,
             mc.vsv$ppc$M0c, mc.vsv$ppc$M1c, mc.vsv$ppc$M2c)


# tnf
ppc <- rbind(mc.tnf.1$ppc$M0, mc.tnf.1$ppc$M1, mc.tnf.1$ppc$M2,
             mc.tnf.2$M0c, mc.tnf.2$M1c, mc.tnf.2$M2c)




# check loo.ic
mc.vsv$ic$loo
mc.lcmv$ic$loo
mc.tnf.1$ic$loo
mc.tnf.2$ic$loo










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
                    ppc$parameter.level == "lowest", ])+
  facet_grid(facets = t~model)+
  geom_errorbar(aes(ymin = y.ppc.L, ymax = y.ppc.H, x = y.real.mean), col = "darkgray")+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(aes(y = y.ppc.mean, x = y.real.mean), shape = 21, fill = "white")+
  theme_bw()

ggplot(data = ppc[ppc$prediction.level == "snp-level" &
                    ppc$parameter.level == "lowest", ])+
  facet_grid(facets = t~model)+
  geom_density2d(aes(y = y.ppc.mean, x = y.real.mean), col = "orange")+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()



#### Parameters: lowest ####

ggplot(data = ppc[ppc$prediction.level == "snp-level" &
                    ppc$parameter.level == "lowest-plus-one", ])+
  facet_grid(facets = t~model)+
  geom_errorbar(aes(ymin = y.ppc.L, ymax = y.ppc.H, x = y.real.mean), col = "darkgray")+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(aes(y = y.ppc.mean, x = y.real.mean), shape = 21, fill = "white")+
  theme_bw()

ggplot(data = ppc[ppc$prediction.level == "snp-level" &
                    ppc$parameter.level == "lowest-plus-one", ])+
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

ggplot(data = ppc[ppc$prediction.level == "snp-level", ])+
  facet_grid(facets = t~model+parameter.level)+
  geom_density2d(aes(y = y.ppc.median, x = y.real.mean), col = "orange")+
  geom_abline(slope = 1, intercept = 0)+
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
