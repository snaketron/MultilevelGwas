default.values <- list(adapt_delta = 0.95,
                       max_treedepth = 10,
                       ntree = 1000,
                       cv.fold = 0.66,
                       refresh = 250,
                       verbose = TRUE,
                       with.sample.file = FALSE)

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
source("R/pareto.R")


mc.ex0 <- get(load(file = "dev.tests/ex0/mc.ex0.RData"))
s <- runStatLearn(gt.data = mc.ex0$gt.data,
                  method = "rf",
                  cv.steps = 100,
                  hdi.level = 0.95,
                  cores = 8,
                  dot.param = default.values)





mc.ex3 <- get(load(file = "dev.tests/ex3/mc.ex3.RData"))
s <- runStatLearn(gt.data = mc.ex3$gt.data,
                  method = "svm",
                  cv.steps = 100,
                  hdi.level = 0.95,
                  cores = 8,
                  dot.param = default.values)





o <- runMLgwas(genotype = mc.ex3$gt.data$genotype,
               traits = mc.ex3$gt.data$Y,
               trait.type = mc.ex3$gt.data$trait.type,
               strains = mc.ex3$gt.data$strains,
               model = "M0",
               mcmc.chains = 2,
               mcmc.steps = 2500,
               mcmc.warmup = 500,
               cores = 8,
               hdi.level = 0.95,
               stat.learn.method = "rf",
               cv.steps = 100)


ggplot(data = o$scores)+
  facet_wrap(facets = ~trait)+
  geom_point(aes(x = k.mean, y = abs(beta.mean), col = front.median))+
  scale_color_gradientn(colours = terrain.colors(n = 10))+
  theme_bw()
