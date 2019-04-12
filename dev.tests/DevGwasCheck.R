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

mc.vsv <- get(load("~/MLgwas/dev.tests/vsv.300.RData"))
mc.lcmv <- get(load("~/MLgwas/dev.tests/lcmv.300.RData"))
mc.tnf <- get(load("~/MLgwas/dev.tests/tnf.300.RData"))
rm(mc)


# x <- getPpc(ps = mc.vsv$ps,
#             gt.data = mc.vsv$gt.data,
#             models = names(mc.vsv$ps),
#             hdi.level = 0.95)
# mc.vsv$ppc <- x
# save(mc.vsv, file = "~/MLgwas/dev.tests/vsv.300.RData")
#
# x <- getPpc(ps = mc.lcmv$ps,
#             gt.data = mc.lcmv$gt.data,
#             models = names(mc.lcmv$ps),
#             hdi.level = 0.95)
# mc.lcmv$ppc <- x
# save(mc.lcmv, file = "~/MLgwas/dev.tests/lcmv.300.RData")


x <- getPpc(ps = mc.tnf$ps,
            gt.data = mc.tnf$gt.data,
            models = names(mc.tnf$ps),
            hdi.level = 0.95)
mc.tnf$ppc <- x
save(mc.tnf, file = "~/MLgwas/dev.tests/tnf.300.RData")


ggplot(data = ppc[ppc$parameter.level == "mid" &
                    ppc$prediction.level == "snp", ])+
  facet_grid(facets = model~t)+
  geom_errorbar(aes(x = s, ymin = y.error.L, ymax = y.error.H))+
  geom_point(aes(x = s, y = y.error.mean))


ggplot(data = ppc[ppc$parameter.level == "mid" &
                  ppc$prediction.level == "snp", ])+
  facet_grid(facets = model~t)+
  geom_errorbar(aes(x = i, ymin = y.error.L, ymax = y.error.H))+
  geom_point(aes(x = i, y = y.error.mean))



ggplot(data = mc.vsv$ppc[mc.vsv$ppc$prediction.level == "individual", ])+
  facet_grid(facets = model~t)+
  geom_errorbar(aes(x = i, ymin = y.error.L, ymax = y.error.H))+
  geom_point(aes(x = i, y = y.error.mean))



ggplot(data = mc.vsv$ppc[mc.vsv$ppc$prediction.level == "strain", ])+
  facet_grid(facets = model~t)+
  geom_errorbar(aes(x = k, ymin = y.error.L, ymax = y.error.H))+
  geom_point(aes(x = k, y = y.error.mean))

ggplot(data = mc.vsv$ppc[mc.vsv$ppc$prediction.level == "snp", ])+
  facet_grid(facets = model~t)+
  geom_errorbar(aes(x = s, ymin = y.error.L, ymax = y.error.H))+
  geom_point(aes(x = s, y = y.error.mean))




















# violin
ggplot(data = mc.vsv$ppc)+
  facet_grid(facets = t~parameter.level+prediction.level, scales = "free")+
  geom_violin(aes(x = model, y = y.error.mean, fill = model))+
  theme_bw()+
  theme(legend.position = "top")+
  xlab(label = '')+
  ylab(label = "Error")

ggplot(data = mc.lcmv$ppc)+
  facet_grid(facets = t~parameter.level+prediction.level, scales = "free")+
  geom_violin(aes(x = model, y = y.error.mean, fill = model))+
  theme_bw()+
  theme(legend.position = "top")+
  xlab(label = '')+
  ylab(label = "Error")

ggplot(data = mc.tnf$ppc)+
  facet_grid(facets = t~parameter.level+prediction.level, scales = "free")+
  geom_violin(aes(x = model, y = y.error.mean, fill = model))+
  theme_bw()+
  theme(legend.position = "top")+
  xlab(label = '')+
  ylab(label = "Error")




# check loo.ic
mc.vsv$ic$loo
mc.lcmv$ic$loo
mc.tnf$ic$loo


