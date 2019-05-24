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



# VSV
mc.vsv <- get(load("~/MLgwas/dev.tests/vsv.300.RData"))

# simplify
gt.data <- mc.vsv$gt.data
ps <- vector(mode = "list", length = 4)
names(ps) <- names(mc.vsv$ps)
for(i in 1:4) {
  ps[[i]] <- rstan::extract(object = mc.vsv$ps[[i]],
                            include = FALSE,
                            pars = c("z", "log_lik",
                                     "log_lik2"))
}
rm(mc.vsv, i)
gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();

# ppc
ppc.mc.vsv <- getPpcMc(ps = ps,
                       gt.data = gt.data,
                       models = names(ps),
                       hdi.level = 0.95,
                       cores = 4,
                       debug = T)
save(ppc.mc.vsv, file = "~/MLgwas/dev.tests/ppc.vsv.300.RData")
rm(mc.vsv, ppc.mc.vsv)




# LCMV
mc.lcmv <- get(load("~/MLgwas/dev.tests/lcmv.300.RData"))

# simplify
gt.data <- mc.lcmv$gt.data
ps <- vector(mode = "list", length = 4)
names(ps) <- names(mc.lcmv$ps)
for(i in 1:4) {
  ps[[i]] <- rstan::extract(object = mc.lcmv$ps[[i]],
                            include = FALSE,
                            pars = c("z", "log_lik",
                                     "log_lik2"))
}
rm(mc.lcmv, i)
gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();

ppc.mc.lcmv <- getPpcMc(ps = ps,
                        gt.data = gt.data,
                        models = names(ps),
                        hdi.level = 0.95,
                        cores = 4,
                        debug = T)
save(ppc.mc.lcmv, file = "~/MLgwas/dev.tests/ppc.lcmv.300.RData")
rm(mc.lcmv, ppc.mc.lcmv)




# TNF
mc.tnf <- get(load("~/MLgwas/dev.tests/tnf.300.RData"))

# simplify
gt.data <- mc.tnf$gt.data
ps <- vector(mode = "list", length = 4)
names(ps) <- names(mc.tnf$ps)
for(i in 1:4) {
  ps[[i]] <- rstan::extract(object = mc.tnf$ps[[i]],
                            include = FALSE,
                            pars = c("z", "log_lik",
                                     "log_lik2"))
}
rm(mc.tnf, i)
gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();

ppc.mc.tnf <- getPpcMc(ps = ps,
                       gt.data = gt.data,
                       models = names(ps),
                       hdi.level = 0.95,
                       cores = 4,
                       debug = T)
save(ppc.mc.tnf, file = "~/MLgwas/dev.tests/ppc.tnf.300.RData")
rm(mc.tnf, ppc.mc.tnf)




ppc.mc.vsv
# violin (low)
ppc.mc.vsv$model <- gsub(pattern = 'c', replacement = '\\*', x = ppc.mc.vsv$model)
vsv.errors.low <- ggplot(data = ppc.mc.vsv[ppc.mc.vsv$parameter.level == "low", ])+
  facet_grid(facets = t~prediction.level, scales = "free")+
  geom_violin(aes(x = model, y = y.error.mean, fill = model), size = 0.3)+
  theme_bw(base_size = 9)+
  theme(legend.position = "top")+
  xlab(label = '')+
  ylab(label = expression(epsilon))
vsv.errors.low <- ggplot(data = ppc.mc.vsv[ppc.mc.vsv$parameter.level == "low", ])+
  facet_grid(facets = t~prediction.level, scales = "free")+
  geom_jitter(aes(x = model, y = y.error.mean, fill = model),  width = 0.2, height = 0.01)+
  theme_bw(base_size = 9)+
  theme(legend.position = "top")+
  xlab(label = '')+
  ylab(label = expression(epsilon))



mc.lcmv$ppc$model <- gsub(pattern = 'c', replacement = '\\*', x = mc.lcmv$ppc$model)
lcmv.errors.low <- ggplot(data = mc.lcmv$ppc[mc.lcmv$ppc$parameter.level == "low", ])+
  facet_grid(facets = t~prediction.level, scales = "free")+
  geom_violin(aes(x = model, y = y.error.mean, fill = model), size = 0.3)+
  theme_bw(base_size = 9)+
  theme(legend.position = "top")+
  xlab(label = '')+
  ylab(label = expression(epsilon))


mc.tnf$ppc$model <- gsub(pattern = 'c', replacement = '\\*', x = mc.tnf$ppc$model)
tnf.errors.low <- ggplot(data = mc.tnf$ppc[mc.tnf$ppc$parameter.level == "low", ])+
  facet_grid(facets = t~prediction.level, scales = "free")+
  geom_violin(aes(x = model, y = y.error.mean, fill = model), size = 0.3)+
  theme_bw(base_size = 9)+
  theme(legend.position = "top")+
  xlab(label = '')+
  ylab(label = expression(epsilon))

ggsave(filename = "out.manuscript/vsv.errors.low.eps",
       plot = vsv.errors.low, device = "eps",
       width = 5, height = 3.75, dpi = 600)
ggsave(filename = "out.manuscript/lcmv.errors.low.eps",
       plot = lcmv.errors.low, device = "eps",
       width = 5, height = 3.75, dpi = 600)
ggsave(filename = "out.manuscript/tnf.errors.low.eps",
       plot = tnf.errors.low, device = "eps",
       width = 5, height = 7, dpi = 600)





# violin (mid)
vsv.errors.mid <- ggplot(data = mc.vsv$ppc[mc.vsv$ppc$parameter.level == "mid" &
                                             mc.vsv$ppc$prediction.level == "snp", ])+
  facet_grid(facets = t~prediction.level, scales = "free")+
  geom_violin(aes(x = model, y = y.error.mean, fill = model), size = 0.3)+
  theme_bw(base_size = 9)+
  theme(legend.position = "top")+
  xlab(label = '')+
  ylab(label = expression(epsilon))

lcmv.errors.mid <- ggplot(data = mc.lcmv$ppc[mc.lcmv$ppc$parameter.level == "mid" &
                                               mc.lcmv$ppc$prediction.level == "snp", ])+
  facet_grid(facets = t~prediction.level, scales = "free")+
  geom_violin(aes(x = model, y = y.error.mean, fill = model), size = 0.3)+
  theme_bw(base_size = 9)+
  theme(legend.position = "top")+
  xlab(label = '')+
  ylab(label = expression(epsilon))

tnf.errors.mid <- ggplot(data = mc.tnf$ppc[mc.tnf$ppc$parameter.level == "mid" &
                                             mc.tnf$ppc$prediction.level == "snp", ])+
  facet_grid(facets = t~prediction.level, scales = "free")+
  geom_violin(aes(x = model, y = y.error.mean, fill = model), size = 0.3)+
  theme_bw(base_size = 9)+
  theme(legend.position = "top")+
  xlab(label = '')+
  ylab(label = expression(epsilon))


pdf(file = "out.manuscript/errors.mid.pdf", width = 10, height = 6.5)
gridExtra::grid.arrange(vsv.errors.mid,
                        lcmv.errors.mid,
                        tnf.errors.mid, nrow = 1)
dev.off()








#### PPC: Vertical -> Individual ####
x.vsv <- mc.vsv$ppc
x.lcmv <- mc.lcmv$ppc
x.tnf <- mc.tnf$ppc
g <- gridExtra::arrangeGrob(ncol = 2,
                            ggplot(data = x.vsv[x.vsv$parameter.level == "low" &
                                                  x.vsv$prediction.level == "individual" &
                                                  x.vsv$t == 1, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L, ymax = y.ppc.H),
                                            width = 0, col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D1: IFN"),
                            ggplot(data = x.vsv[x.vsv$parameter.level == "low" &
                                                  x.vsv$prediction.level == "individual" &
                                                  x.vsv$t == 2, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L, ymax = y.ppc.H),
                                            width = 0, col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D1: Replication"),
                            ggplot(data = x.lcmv[x.lcmv$parameter.level == "low" &
                                                   x.lcmv$prediction.level == "individual" &
                                                   x.lcmv$t == 1, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L, ymax = y.ppc.H),
                                            width = 0, col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D2: IFN"),
                            ggplot(data = x.lcmv[x.lcmv$parameter.level == "low" &
                                                   x.lcmv$prediction.level == "individual" &
                                                   x.lcmv$t == 2, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L, ymax = y.ppc.H),
                                            width = 0, col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D2: Replication"),
                            ggplot(data = x.tnf[x.tnf$parameter.level == "low" &
                                                  x.tnf$prediction.level == "individual" &
                                                  x.tnf$t == 1, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L, ymax = y.ppc.H),
                                            width = 0, col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D3: ALT"),
                            ggplot(data = x.tnf[x.tnf$parameter.level == "low" &
                                                  x.tnf$prediction.level == "individual" &
                                                  x.tnf$t == 2, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L, ymax = y.ppc.H),
                                            width = 0, col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D3: AST"),
                            ggplot(data = x.tnf[x.tnf$parameter.level == "low" &
                                                  x.tnf$prediction.level == "individual" &
                                                  x.tnf$t == 3, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L, ymax = y.ppc.H),
                                            width = 0, col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D3: LDH"),
                            ggplot(data = x.tnf[x.tnf$parameter.level == "low" &
                                                  x.tnf$prediction.level == "individual" &
                                                  x.tnf$t == 4, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L, ymax = y.ppc.H),
                                            width = 0, col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              scale_x_continuous(breaks = c(0, 0.5, 1),
                                                 labels = c(0, 0.5, 1),
                                                 limits = c(0, 1))+
                              scale_y_continuous(breaks = c(0, 0.5, 1),
                                                 labels = c(0, 0.5, 1),
                                                 limits = c(0, 1))+
                              ggtitle(label = "D3: Survival"))
plot(g)
ggsave(filename = "out.manuscript/VerticalIndividual.jpg", plot = g,
       device = "jpg", width = 4.5, height = 9, dpi = 600)
ggsave(filename = "out.manuscript/VerticalIndividual.eps", plot = g,
       device = "eps", width = 4.5, height = 9, dpi = 600)





#### PPC: Vertical -> Strain ####
x.vsv <- mc.vsv$ppc
x.vsv$model <- gsub(pattern = 'c', replacement = '\\*', x = x.vsv$model)

x.lcmv <- mc.lcmv$ppc
x.lcmv$model <- gsub(pattern = 'c', replacement = '\\*', x = x.lcmv$model)

x.tnf <- mc.tnf$ppc
x.tnf$model <- gsub(pattern = 'c', replacement = '\\*', x = x.tnf$model)


g <- gridExtra::arrangeGrob(ncol = 2,
                            ggplot(data = x.vsv[x.vsv$parameter.level == "low" &
                                                  x.vsv$prediction.level == "strain" &
                                                  x.vsv$t == 1, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L, ymax = y.ppc.H),
                                            width = 0, col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D1: IFN"),
                            ggplot(data = x.vsv[x.vsv$parameter.level == "low" &
                                                  x.vsv$prediction.level == "strain" &
                                                  x.vsv$t == 2, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L, ymax = y.ppc.H),
                                            width = 0, col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D1: Replication"),
                            ggplot(data = x.lcmv[x.lcmv$parameter.level == "low" &
                                                   x.lcmv$prediction.level == "strain" &
                                                   x.lcmv$t == 1, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L, ymax = y.ppc.H),
                                            width = 0, col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D2: IFN"),
                            ggplot(data = x.lcmv[x.lcmv$parameter.level == "low" &
                                                   x.lcmv$prediction.level == "strain" &
                                                   x.lcmv$t == 2, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L, ymax = y.ppc.H),
                                            width = 0, col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D2: Replication"),
                            ggplot(data = x.tnf[x.tnf$parameter.level == "low" &
                                                  x.tnf$prediction.level == "strain" &
                                                  x.tnf$t == 1, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L, ymax = y.ppc.H),
                                            width = 0, col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D3: ALT"),
                            ggplot(data = x.tnf[x.tnf$parameter.level == "low" &
                                                  x.tnf$prediction.level == "strain" &
                                                  x.tnf$t == 2, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L, ymax = y.ppc.H),
                                            width = 0, col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D3: AST"),
                            ggplot(data = x.tnf[x.tnf$parameter.level == "low" &
                                                  x.tnf$prediction.level == "strain" &
                                                  x.tnf$t == 3, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L, ymax = y.ppc.H),
                                            width = 0, col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D3: LDH"),
                            ggplot(data = x.tnf[x.tnf$parameter.level == "low" &
                                                  x.tnf$prediction.level == "strain" &
                                                  x.tnf$t == 4, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L, ymax = y.ppc.H),
                                            width = 0, col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              scale_x_continuous(breaks = c(0, 0.5, 1),
                                                 labels = c(0, 0.5, 1),
                                                 limits = c(0, 1))+
                              scale_y_continuous(breaks = c(0, 0.5, 1),
                                                 labels = c(0, 0.5, 1),
                                                 limits = c(0, 1))+
                              ggtitle(label = "D3: Survival"))
plot(g)
ggsave(filename = "out.manuscript/VerticalStrain.jpg", plot = g,
       device = "jpg", width = 4.5, height = 9, dpi = 600)
ggsave(filename = "out.manuscript/VerticalStrain.eps", plot = g,
       device = "eps", width = 4.5, height = 9, dpi = 600)





#### PPC: Vertical -> SNP ####
x.vsv <- mc.vsv$ppc
x.lcmv <- mc.lcmv$ppc
x.tnf <- mc.tnf$ppc
g <- gridExtra::arrangeGrob(ncol = 2,
                            ggplot(data = x.vsv[x.vsv$parameter.level == "low" &
                                                  x.vsv$prediction.level == "snp" &
                                                  x.vsv$t == 1, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L,
                                                ymax = y.ppc.H), col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              geom_density2d(aes(x = y.real.mean, y = y.ppc.mean),
                                             col = "orange", bins = 5, size = 0.15)+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              scale_x_continuous(breaks = c(3, 3.5, 4.0, 4.5),
                                                 labels = c(3, 3.5, 4.0, 4.5),
                                                 limits = c(3, 4.5))+
                              scale_y_continuous(breaks = c(3, 3.5, 4.0, 4.5),
                                                 labels = c(3, 3.5, 4.0, 4.5),
                                                 limits = c(3, 4.5))+
                              ggtitle(label = "D1: IFN"),
                            ggplot(data = x.vsv[x.vsv$parameter.level == "low" &
                                                  x.vsv$prediction.level == "snp" &
                                                  x.vsv$t == 2, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L,
                                                ymax = y.ppc.H), col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              geom_density2d(aes(x = y.real.mean, y = y.ppc.mean),
                                             col = "orange", bins = 5, size = 0.15)+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              scale_x_continuous(breaks = c(5, 6, 7, 8),
                                                 labels = c(5, 6, 7, 8))+
                              scale_y_continuous(breaks = c(5, 6, 7, 8),
                                                 labels = c(5, 6, 7, 8))+
                              ggtitle(label = "D1: Replication"),
                            ggplot(data = x.lcmv[x.lcmv$parameter.level == "low" &
                                                   x.lcmv$prediction.level == "snp" &
                                                   x.lcmv$t == 1, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L,
                                                ymax = y.ppc.H), col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              geom_density2d(aes(x = y.real.mean, y = y.ppc.mean),
                                             col = "orange", bins = 5, size = 0.15)+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D2: IFN"),
                            ggplot(data = x.lcmv[x.lcmv$parameter.level == "low" &
                                                   x.lcmv$prediction.level == "snp" &
                                                   x.lcmv$t == 2, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L,
                                                ymax = y.ppc.H), col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              geom_density2d(aes(x = y.real.mean, y = y.ppc.mean),
                                             col = "orange", bins = 5, size = 0.15)+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D2: Replication"),
                            ggplot(data = x.tnf[x.tnf$parameter.level == "low" &
                                                  x.tnf$prediction.level == "snp" &
                                                  x.tnf$t == 1, ])+
                               facet_wrap(facets = ~model, nrow = 2)+
                               geom_abline(intercept = 0, slope = 1)+
                               geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L,
                                                 ymax = y.ppc.H), col = "darkgray", size = 0.3)+
                               geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                          shape = 21, size = 0.5, stroke = 0.3,
                                          col = "black", fill = "white")+
                               geom_density2d(aes(x = y.real.mean, y = y.ppc.mean),
                                              col = "orange", bins = 5, size = 0.15)+
                               theme_bw(base_size = 9)+
                               xlab(label = "Y")+
                               ylab(label = "Yhat")+
                              scale_x_continuous(breaks = c(2, 2.4, 2.8, 3.2),
                                                 labels = c(2, 2.4, 2.8, 3.2))+
                              scale_y_continuous(breaks = c(2, 2.4, 2.8, 3.2),
                                                 labels = c(2, 2.4, 2.8, 3.2))+
                               ggtitle(label = "D3: ALT"),
                            ggplot(data = x.tnf[x.tnf$parameter.level == "low" &
                                                  x.tnf$prediction.level == "snp" &
                                                  x.tnf$t == 2, ])+
                               facet_wrap(facets = ~model, nrow = 2)+
                               geom_abline(intercept = 0, slope = 1)+
                               geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L,
                                                 ymax = y.ppc.H), col = "darkgray", size = 0.3)+
                               geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                          shape = 21, size = 0.5, stroke = 0.3,
                                          col = "black", fill = "white")+
                               geom_density2d(aes(x = y.real.mean, y = y.ppc.mean),
                                              col = "orange", bins = 5, size = 0.15)+
                               theme_bw(base_size = 9)+
                               xlab(label = "Y")+
                               ylab(label = "Yhat")+
                               ggtitle(label = "D3: AST"),
                            ggplot(data = x.tnf[x.tnf$parameter.level == "low" &
                                                  x.tnf$prediction.level == "snp" &
                                                  x.tnf$t == 3, ])+
                               facet_wrap(facets = ~model, nrow = 2)+
                               geom_abline(intercept = 0, slope = 1)+
                               geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L,
                                                 ymax = y.ppc.H), col = "darkgray", size = 0.3)+
                               geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                          shape = 21, size = 0.5, stroke = 0.3,
                                          col = "black", fill = "white")+
                               geom_density2d(aes(x = y.real.mean, y = y.ppc.mean),
                                              col = "orange", bins = 5, size = 0.15)+
                               theme_bw(base_size = 9)+
                               xlab(label = "Y")+
                               ylab(label = "Yhat")+
                              scale_x_continuous(breaks = c(3.5, 3.75, 4.0, 4.25),
                                                 labels = c(3.5, 3.75, 4.0, 4.25))+
                              scale_y_continuous(breaks = c(3.5, 3.75, 4.0, 4.25),
                                                 labels = c(3.5, 3.75, 4.0, 4.25))+
                               ggtitle(label = "D3: LDH"),
                            ggplot(data = x.tnf[x.tnf$parameter.level == "low" &
                                                  x.tnf$prediction.level == "snp" &
                                                  x.tnf$t == 4, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L,
                                                ymax = y.ppc.H), col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              geom_density2d(aes(x = y.real.mean, y = y.ppc.mean),
                                             col = "orange", bins = 5, size = 0.15)+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              scale_x_continuous(breaks = c(0, 0.5, 1),
                                                 labels = c(0, 0.5, 1),
                                                 limits = c(0, 1))+
                              scale_y_continuous(breaks = c(0, 0.5, 1),
                                                 labels = c(0, 0.5, 1),
                                                 limits = c(0, 1))+
                              ggtitle(label = "D3: Survival"))
plot(g)
ggsave(filename = "out.manuscript/VerticalSnp.jpg", plot = g,
       device = "jpg", width = 4.5, height = 9, dpi = 600)
ggsave(filename = "out.manuscript/VerticalSnp.eps", plot = g,
       device = "eps", width = 4.5, height = 9, dpi = 600)




#### PPC: Horizontal -> SNP ####
x.vsv <- mc.vsv$ppc
x.lcmv <- mc.lcmv$ppc
x.tnf <- mc.tnf$ppc
g <- gridExtra::arrangeGrob(ncol = 2,
                            ggplot(data = x.vsv[x.vsv$parameter.level == "mid" &
                                                  x.vsv$prediction.level == "snp" &
                                                  x.vsv$t == 1, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L,
                                                ymax = y.ppc.H), col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              geom_density2d(aes(x = y.real.mean, y = y.ppc.mean),
                                             col = "orange", bins = 5, size = 0.15)+
                              scale_x_continuous(breaks = c(3, 3.5, 4.0, 4.5, 5),
                                                 labels = c(3, 3.5, 4.0, 4.5, 5),
                                                 limits = c(3, 5))+
                              scale_y_continuous(breaks = c(3, 3.5, 4.0, 4.5, 5),
                                                 labels = c(3, 3.5, 4.0, 4.5, 5),
                                                 limits = c(3, 5))+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D1: IFN"),
                            ggplot(data = x.vsv[x.vsv$parameter.level == "mid" &
                                                  x.vsv$prediction.level == "snp" &
                                                  x.vsv$t == 2, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L,
                                                ymax = y.ppc.H), col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              geom_density2d(aes(x = y.real.mean, y = y.ppc.mean),
                                             col = "orange", bins = 5, size = 0.15)+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D1: Replication"),
                            ggplot(data = x.lcmv[x.lcmv$parameter.level == "mid" &
                                                   x.lcmv$prediction.level == "snp" &
                                                   x.lcmv$t == 1, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L,
                                                ymax = y.ppc.H), col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              geom_density2d(aes(x = y.real.mean, y = y.ppc.mean),
                                             col = "orange", bins = 5, size = 0.15)+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D2: IFN"),
                            ggplot(data = x.lcmv[x.lcmv$parameter.level == "mid" &
                                                   x.lcmv$prediction.level == "snp" &
                                                   x.lcmv$t == 2, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L,
                                                ymax = y.ppc.H), col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              geom_density2d(aes(x = y.real.mean, y = y.ppc.mean),
                                             col = "orange", bins = 5, size = 0.15)+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D2: Replication"),
                            ggplot(data = x.tnf[x.tnf$parameter.level == "mid" &
                                                  x.tnf$prediction.level == "snp" &
                                                  x.tnf$t == 1, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L,
                                                ymax = y.ppc.H), col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              geom_density2d(aes(x = y.real.mean, y = y.ppc.mean),
                                             col = "orange", bins = 5, size = 0.15)+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D3: ALT"),
                            ggplot(data = x.tnf[x.tnf$parameter.level == "mid" &
                                                  x.tnf$prediction.level == "snp" &
                                                  x.tnf$t == 2, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L,
                                                ymax = y.ppc.H), col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              geom_density2d(aes(x = y.real.mean, y = y.ppc.mean),
                                             col = "orange", bins = 5, size = 0.15)+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D3: AST"),
                            ggplot(data = x.tnf[x.tnf$parameter.level == "mid" &
                                                  x.tnf$prediction.level == "snp" &
                                                  x.tnf$t == 3, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L,
                                                ymax = y.ppc.H), col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              geom_density2d(aes(x = y.real.mean, y = y.ppc.mean),
                                             col = "orange", bins = 5, size = 0.15)+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D3: LDH"),
                            ggplot(data = x.tnf[x.tnf$parameter.level == "mid" &
                                                  x.tnf$prediction.level == "snp" &
                                                  x.tnf$t == 4, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L,
                                                ymax = y.ppc.H), col = "darkgray", size = 0.3)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 0.5, stroke = 0.3,
                                         col = "black", fill = "white")+
                              geom_density2d(aes(x = y.real.mean, y = y.ppc.mean),
                                             col = "orange", bins = 5, size = 0.15)+
                              scale_x_continuous(breaks = c(0, 0.5, 1),
                                                 labels = c(0, 0.5, 1),
                                                 limits = c(0, 1))+
                              scale_y_continuous(breaks = c(0, 0.5, 1),
                                                 labels = c(0, 0.5, 1),
                                                 limits = c(0, 1))+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D3: Survival"))
plot(g)
ggsave(filename = "out.manuscript/HorizontalSnp.jpg", plot = g,
       device = "jpg", width = 4.5, height = 9, dpi = 600)
ggsave(filename = "out.manuscript/HorizontalSnp.eps", plot = g,
       device = "eps", width = 4.5, height = 9, dpi = 600)








#### PPC: Horizontal -> Strain ####
x.vsv <- mc.vsv$ppc
x.vsv <- x.vsv[x.vsv$parameter.level == "mid" &
                 x.vsv$prediction.level == "strain", ]
a <- merge(x = aggregate(y.ppc.mean~t+k+model, data = x.vsv, FUN = mean),
           y = aggregate(y.real.mean~t+k+model, data = x.vsv, FUN = mean),
           by = c("t", "k", "model"))
b <- merge(x = aggregate(y.ppc.L~t+k+model, data = x.vsv, FUN = mean),
           y = aggregate(y.ppc.H~t+k+model, data = x.vsv, FUN = mean),
           by = c("t", "k", "model"))
x.vsv  <- merge(x = a, y = b, by = c("t", "k", "model"))
rm(a, b)

x.lcmv <- mc.lcmv$ppc
x.lcmv <- x.lcmv[x.lcmv$parameter.level == "mid" &
                 x.lcmv$prediction.level == "strain", ]
a <- merge(x = aggregate(y.ppc.mean~t+k+model, data = x.lcmv, FUN = mean),
           y = aggregate(y.real.mean~t+k+model, data = x.lcmv, FUN = mean),
           by = c("t", "k", "model"))
b <- merge(x = aggregate(y.ppc.L~t+k+model, data = x.lcmv, FUN = mean),
           y = aggregate(y.ppc.H~t+k+model, data = x.lcmv, FUN = mean),
           by = c("t", "k", "model"))
x.lcmv  <- merge(x = a, y = b, by = c("t", "k", "model"))
rm(a, b)


x.tnf <- mc.tnf$ppc
x.tnf <- x.tnf[x.tnf$parameter.level == "mid" &
                 x.tnf$prediction.level == "strain", ]
a <- merge(x = aggregate(y.ppc.mean~t+k+model, data = x.tnf, FUN = mean),
           y = aggregate(y.real.mean~t+k+model, data = x.tnf, FUN = mean),
           by = c("t", "k", "model"))
b <- merge(x = aggregate(y.ppc.L~t+k+model, data = x.tnf, FUN = mean),
           y = aggregate(y.ppc.H~t+k+model, data = x.tnf, FUN = mean),
           by = c("t", "k", "model"))
x.tnf  <- merge(x = a, y = b, by = c("t", "k", "model"))
rm(a, b)



g <- gridExtra::arrangeGrob(ncol = 2,
                            ggplot(data = x.vsv[x.vsv$t == 1, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L, ymax = y.ppc.H),
                                            col = "darkgray", size = 0.5, width = 0)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 1, stroke = 0.5,
                                         col = "black", fill = "white")+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D1: IFN"),
                            ggplot(data = x.vsv[x.vsv$t == 2, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L, ymax = y.ppc.H),
                                            col = "darkgray", size = 0.5, width = 0)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 1, stroke = 0.5,
                                         col = "black", fill = "white")+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D1: Replication"),
                            ggplot(data = x.lcmv[x.lcmv$t == 1, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L, ymax = y.ppc.H),
                                            col = "darkgray", size = 0.5, width = 0)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 1, stroke = 0.5,
                                         col = "black", fill = "white")+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D2: IFN"),
                            ggplot(data = x.lcmv[x.lcmv$t == 2, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L, ymax = y.ppc.H),
                                            col = "darkgray", size = 0.5, width = 0)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 1, stroke = 0.5,
                                         col = "black", fill = "white")+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D2: Replication"),
                            ggplot(data = x.tnf[x.tnf$t == 1, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L, ymax = y.ppc.H),
                                            col = "darkgray", size = 0.5, width = 0)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 1, stroke = 0.5,
                                         col = "black", fill = "white")+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D3: ALT"),
                            ggplot(data = x.tnf[x.tnf$t == 2, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L, ymax = y.ppc.H),
                                            col = "darkgray", size = 0.5, width = 0)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 1, stroke = 0.5,
                                         col = "black", fill = "white")+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D3: AST"),
                            ggplot(data = x.tnf[x.tnf$t == 3, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L, ymax = y.ppc.H),
                                            col = "darkgray", size = 0.5, width = 0)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 1, stroke = 0.5,
                                         col = "black", fill = "white")+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D3: LDH"),
                            ggplot(data = x.tnf[x.tnf$t == 4, ])+
                              facet_wrap(facets = ~model, nrow = 2)+
                              geom_abline(intercept = 0, slope = 1)+
                              geom_errorbar(aes(x = y.real.mean, ymin = y.ppc.L, ymax = y.ppc.H),
                                            col = "darkgray", size = 0.5, width = 0)+
                              geom_point(aes(x = y.real.mean, y = y.ppc.mean),
                                         shape = 21, size = 1, stroke = 0.5,
                                         col = "black", fill = "white")+
                              scale_x_continuous(breaks = c(0, 0.5, 1),
                                                 labels = c(0, 0.5, 1),
                                                 limits = c(0, 1))+
                              scale_y_continuous(breaks = c(0, 0.5, 1),
                                                 labels = c(0, 0.5, 1),
                                                 limits = c(0, 1))+
                              theme_bw(base_size = 9)+
                              xlab(label = "Y")+
                              ylab(label = "Yhat")+
                              ggtitle(label = "D3: Survival"))
plot(g)
ggsave(filename = "out.manuscript/HorizontalStrain.jpg", plot = g,
       device = "jpg", width = 4.5, height = 9, dpi = 600)
ggsave(filename = "out.manuscript/HorizontalStrain.eps", plot = g,
       device = "eps", width = 4.5, height = 9, dpi = 600)














#### PPC: Horizontal -> SNP 95% HDI length ####
x.vsv <- mc.vsv$ppc
x.lcmv <- mc.lcmv$ppc
x.tnf <- mc.tnf$ppc

g <- gridExtra::arrangeGrob(ncol = 2,
                            ggplot(data = x.vsv[x.vsv$parameter.level == "mid" &
                                                  x.vsv$prediction.level == "snp" &
                                                  x.vsv$t == 1, ])+
                              geom_density(aes(y.ppc.H-y.ppc.L, col = model))+
                              theme_bw(base_size = 9)+
                              theme(legend.position =  c(1.25, 1.25), legend.direction = "horizontal")+
                              xlab(label = "Predictive uncertainty (Length of 95% HDI)")+
                              ylab(label = "Density")+
                              ggtitle(label = "D1: IFN"),
                            ggplot(data = x.vsv[x.vsv$parameter.level == "mid" &
                                                  x.vsv$prediction.level == "snp" &
                                                  x.vsv$t == 2, ])+
                              geom_density(aes(y.ppc.H-y.ppc.L, col = model))+
                              theme_bw(base_size = 9)+
                              theme(legend.position = "none")+
                              xlab(label = "Predictive uncertainty (Length of 95% HDI)")+
                              ylab(label = "Density")+
                              ggtitle(label = "D1: Replication"),
                            ggplot(data = x.lcmv[x.lcmv$parameter.level == "mid" &
                                                  x.lcmv$prediction.level == "snp" &
                                                  x.lcmv$t == 1, ])+
                              geom_density(aes(y.ppc.H-y.ppc.L, col = model))+
                              theme_bw(base_size = 9)+
                              theme(legend.position = "none")+
                              xlab(label = "Predictive uncertainty (Length of 95% HDI)")+
                              ylab(label = "Density")+
                              ggtitle(label = "D2: IFN"),
                            ggplot(data = x.lcmv[x.lcmv$parameter.level == "mid" &
                                                   x.lcmv$prediction.level == "snp" &
                                                   x.lcmv$t == 1, ])+
                              geom_density(aes(y.ppc.H-y.ppc.L, col = model))+
                              theme_bw(base_size = 9)+
                              theme(legend.position = "none")+
                              xlab(label = "Predictive uncertainty (Length of 95% HDI)")+
                              ylab(label = "Density")+
                              ggtitle(label = "D2: Replication"),
                            ggplot(data = x.lcmv[x.lcmv$parameter.level == "mid" &
                                                  x.lcmv$prediction.level == "snp" &
                                                  x.lcmv$t == 1, ])+
                              geom_density(aes(y.ppc.H-y.ppc.L, col = model))+
                              theme_bw(base_size = 9)+
                              theme(legend.position = "none")+
                              xlab(label = "Predictive uncertainty (Length of 95% HDI)")+
                              ylab(label = "Density")+
                              ggtitle(label = "D3: ALT"),
                            ggplot(data = x.tnf[x.tnf$parameter.level == "mid" &
                                                   x.tnf$prediction.level == "snp" &
                                                   x.tnf$t == 2, ])+
                              geom_density(aes(y.ppc.H-y.ppc.L, col = model))+
                              theme_bw(base_size = 9)+
                              theme(legend.position = "none")+
                              xlab(label = "Predictive uncertainty (Length of 95% HDI)")+
                              ylab(label = "Density")+
                              ggtitle(label = "D3: AST"),
                            ggplot(data = x.tnf[x.tnf$parameter.level == "mid" &
                                                   x.tnf$prediction.level == "snp" &
                                                   x.tnf$t == 3, ])+
                              geom_density(aes(y.ppc.H-y.ppc.L, col = model))+
                              theme_bw(base_size = 9)+
                              theme(legend.position = "none")+
                              xlab(label = "Predictive uncertainty (Length of 95% HDI)")+
                              ylab(label = "Density")+
                              ggtitle(label = "D3: LDH"),
                            ggplot(data = x.tnf[x.tnf$parameter.level == "mid" &
                                                   x.tnf$prediction.level == "snp" &
                                                   x.tnf$t == 4, ])+
                              geom_density(aes(y.ppc.H-y.ppc.L, col = model))+
                              theme_bw(base_size = 9)+
                              theme(legend.position = "none")+
                              xlab(label = "Predictive uncertainty (Length of 95% HDI)")+
                              ylab(label = "Density")+
                              ggtitle(label = "D3: Survival"))
plot(g)

ggsave(filename = "out.manuscript/HorizontalSnpHDI.jpg", plot = g,
       device = "jpg", width = 5.5, height = 6, dpi = 600)
ggsave(filename = "out.manuscript/HorizontalSnpHDI.eps", plot = g,
       device = "eps", width = 5.5, height = 6, dpi = 600)

