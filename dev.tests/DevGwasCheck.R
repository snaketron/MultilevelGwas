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



# model comparison

# VSV
loo::compare(mc.vsv$ic$loo$M0, mc.vsv$ic$loo$M0c)
loo::compare(mc.vsv$ic$loo$M0, mc.vsv$ic$loo$M1)
loo::compare(mc.vsv$ic$loo$M0, mc.vsv$ic$loo$M1c)
loo::compare(mc.vsv$ic$loo$M0c, mc.vsv$ic$loo$M1)
loo::compare(mc.vsv$ic$loo$M0c, mc.vsv$ic$loo$M1c)
loo::compare(mc.vsv$ic$loo$M1, mc.vsv$ic$loo$M1c)

# LCMV
loo::compare(mc.lcmv$ic$loo$M0, mc.lcmv$ic$loo$M0c)
loo::compare(mc.lcmv$ic$loo$M0, mc.lcmv$ic$loo$M1)
loo::compare(mc.lcmv$ic$loo$M0, mc.lcmv$ic$loo$M1c)
loo::compare(mc.lcmv$ic$loo$M0c, mc.lcmv$ic$loo$M1)
loo::compare(mc.lcmv$ic$loo$M0c, mc.lcmv$ic$loo$M1c)
loo::compare(mc.lcmv$ic$loo$M1, mc.lcmv$ic$loo$M1c)

# TNF
loo::compare(mc.tnf$ic$loo$M0, mc.tnf$ic$loo$M0c)
loo::compare(mc.tnf$ic$loo$M0, mc.tnf$ic$loo$M1)
loo::compare(mc.tnf$ic$loo$M0, mc.tnf$ic$loo$M1c)
loo::compare(mc.tnf$ic$loo$M0c, mc.tnf$ic$loo$M1)
loo::compare(mc.tnf$ic$loo$M0c, mc.tnf$ic$loo$M1c)
loo::compare(mc.tnf$ic$loo$M1, mc.tnf$ic$loo$M1c)


loo::extract_log_lik(mc.vsv$ps$M0, parameter_name = "log_lik2")


# x <- getPpcMc(ps = mc.tnf.1$ps,
#               gt.data = mc.tnf.1$gt.data,
#               models = names(mc.tnf.1$ps),
#               hdi.level = 0.95, cores = 3)
# mc.tnf.1$ppc.alt <- x
# save(mc.tnf.1, file = "~/MLgwas/dev.tests/tnf.300.RData")



# x <- getPpcMc(ps = mc.tnf.2$ps,
#               gt.data = mc.tnf.2$gt.data,
#               models = names(mc.tnf.2$ps),
#               hdi.level = 0.95, cores = 6)
# mc.tnf.2$ppc.alt <- x
# save(mc.tnf.2, file = "~/MLgwas/dev.tests/tnf.cov.300.RData")



# x <- getPpcMc(ps = mc.vsv$ps,
#               gt.data = mc.vsv$gt.data,
#               models = names(mc.vsv$ps),
#               hdi.level = 0.95, cores = 6)
# mc.vsv$ppc.alt <- x
# save(mc.vsv, file = "~/MLgwas/dev.tests/vsv.300.RData")



# x <- getPpcMc(ps = mc.lcmv$ps,
#               gt.data = mc.lcmv$gt.data,
#               models = names(mc.lcmv$ps),
#               hdi.level = 0.95, cores = 6)
# mc.lcmv$ppc.alt <- x
# save(mc.lcmv, file = "~/MLgwas/dev.tests/lcmv.300.RData")




# form PPC data
# vsv
ppc.vsv <- rbind(mc.vsv$ppc$M0, mc.vsv$ppc$M1,
                 mc.vsv$ppc$M0c, mc.vsv$ppc$M1c)
# lcmv
ppc.lcmv <- rbind(mc.lcmv$ppc$M0, mc.lcmv$ppc$M1,
                  mc.lcmv$ppc$M0c, mc.lcmv$ppc$M1c)
# tnf
ppc.tnf <- rbind(mc.tnf$ppc$M0, mc.tnf$ppc$M1,
                 mc.tnf$ppc$M0c, mc.tnf$ppc$M1c)




# check loo.ic
mc.vsv$ic$loo
mc.lcmv$ic$loo
mc.tnf.1$ic$loo
mc.tnf.2$ic$loo


# compute error
ppc.vsv$e <- abs(ppc.vsv$y.real.mean-ppc.vsv$y.ppc.mean)
ppc.lcmv$e <- abs(ppc.lcmv$y.real.mean-ppc.lcmv$y.ppc.mean)
ppc.tnf$e <- abs(ppc.tnf$y.real.mean-ppc.tnf$y.ppc.mean)





# PPC density


getPpcViz <- function(ppc, nrow) {
  # # Prediction: individual-level
  # ppc.gi <- ppc[ppc$s == 1, ]
  #
  # g1 <- ggplot(data = ppc.gi[ppc.gi$prediction.level == "individual-level", ])+
  #   facet_wrap(facets = t~model, nrow = nrow, scales = "free")+
  #   geom_errorbar(aes(ymin = y.ppc.L, ymax = y.ppc.H, x = y.real.mean), col = "darkgray")+
  #   geom_abline(slope = 1, intercept = 0)+
  #   geom_point(aes(y = y.ppc.mean, x = y.real.mean), shape = 21, fill = "white")+
  #   theme_bw()
  #
  # g1 <- ggplot(data = ppc[ppc$prediction.level == "individual-level", ])+
  #   facet_wrap(facets = t~model, nrow = nrow, scales = "free")+
  #   geom_density2d(aes(y = y.ppc.mean, x = y.real.mean), col = "orange")+
  #   geom_abline(slope = 1, intercept = 0)+
  #   theme_bw()
  # ggsave(filename = "dev.tests/plots/lcmv.individual.lowest.pdf", plot = g1,
  #        device = "pdf", width = 10, height = 6)


  # Prediction: strain-level
  # g2 <- ggplot(data = ppc[ppc$prediction.level == "strain-level", ])+
  #   facet_wrap(facets = t~model, nrow = nrow, scales = "free")+
  #   geom_errorbar(aes(ymin = y.ppc.L, ymax = y.ppc.H, x = y.real.mean), col = "darkgray")+
  #   geom_abline(slope = 1, intercept = 0)+
  #   geom_point(aes(y = y.ppc.mean, x = y.real.mean), shape = 21, fill = "white")+
  #   theme_bw()+
  #   ggtitle(label = "Prediction:  Trait (strain-level), Param (lowest-level)")

  g2 <- ggplot(data = ppc[ppc$prediction.level == "strain-level", ])+
    facet_wrap(facets = t~model, nrow = nrow, scales = "free")+
    geom_density2d(aes(y = y.ppc.mean, x = y.real.mean), col = "orange")+
    geom_abline(slope = 1, intercept = 0)+
    theme_bw()+
    ggtitle(label = "Prediction:  Trait (strain-level), Param (lowest-level)")



  # Prediction: snp-level
  # g3 <- ggplot(data = ppc[ppc$prediction.level == "snp-level" &
  #                           ppc$parameter.level == "lowest", ])+
  #   facet_wrap(facets = t~model, nrow = nrow, scales = "free")+
  #   geom_errorbar(aes(ymin = y.ppc.L, ymax = y.ppc.H, x = y.real.mean), col = "darkgray")+
  #   geom_abline(slope = 1, intercept = 0)+
  #   geom_point(aes(y = y.ppc.mean, x = y.real.mean), shape = 21, fill = "white")+
  #   theme_bw()

  g3 <- ggplot(data = ppc[ppc$prediction.level == "snp-level" &
                            ppc$parameter.level == "lowest", ])+
    facet_wrap(facets = t~model, nrow = nrow, scales = "free")+
    geom_density2d(aes(y = y.ppc.mean, x = y.real.mean), col = "orange")+
    geom_abline(slope = 1, intercept = 0)+
    theme_bw()+
    ggtitle(label = "Prediction:  Trait (snp-level), Param (lowest-level)")


  # Prediction: snp-level (mid level)
  # g4 <- ggplot(data = ppc[ppc$prediction.level == "snp-level" &
  #                           ppc$parameter.level == "lowest-plus-one", ])+
  #   facet_wrap(facets = t~model, nrow = nrow, scales = "free")+
  #   geom_errorbar(aes(ymin = y.ppc.L, ymax = y.ppc.H, x = y.real.mean), col = "darkgray")+
  #   geom_abline(slope = 1, intercept = 0)+
  #   geom_point(aes(y = y.ppc.mean, x = y.real.mean), shape = 21, fill = "white")+
  #   theme_bw()

  g4 <-  ggplot(data = ppc[ppc$prediction.level == "snp-level" &
                             ppc$parameter.level == "lowest-plus-one", ])+
    facet_wrap(facets = t~model, nrow = nrow, scales = "free")+
    geom_density2d(aes(y = y.ppc.mean, x = y.real.mean), col = "orange")+
    geom_abline(slope = 1, intercept = 0)+
    theme_bw()+
    ggtitle(label = "Prediction:  Trait (snp-level), Param (mid-level)")

  return(list(strain = g2, snp.lowest = g3, snp.mid = g4))
}



# VSV
plots.vsv <- getPpcViz(ppc = ppc.vsv, nrow = 2)
ggsave(filename = "dev.tests/plots/vsv.strain.lowest.pdf",
       plot = plots.vsv$strain, device = "pdf", width = 10, height = 6)
ggsave(filename = "dev.tests/plots/vsv.snp.lowest.pdf",
       plot = plots.vsv$snp.lowest, device = "pdf", width = 10, height = 6)
ggsave(filename = "dev.tests/plots/vsv.snp.mid.pdf",
       plot = plots.vsv$snp.mid, device = "pdf", width = 10, height = 6)




# LCMV
plots.lcmv <- getPpcViz(ppc = ppc.lcmv, nrow = 2)
ggsave(filename = "dev.tests/plots/lcmv.strain.lowest.pdf",
       plot = plots.lcmv$strain, device = "pdf", width = 12, height = 6)
ggsave(filename = "dev.tests/plots/lcmv.snp.lowest.pdf",
       plot = plots.lcmv$snp.lowest, device = "pdf", width = 12, height = 6)
ggsave(filename = "dev.tests/plots/lcmv.snp.mid.pdf",
       plot = plots.lcmv$snp.mid, device = "pdf", width = 12, height = 6)




# TNF
plots.tnf <- getPpcViz(ppc = ppc.tnf, nrow = 4)
ggsave(filename = "dev.tests/plots/tnf.strain.lowest.pdf",
       plot = plots.tnf$strain, device = "pdf", width = 12, height = 9)
ggsave(filename = "dev.tests/plots/tnf.snp.lowest.pdf",
       plot = plots.tnf$snp.lowest, device = "pdf", width = 12, height = 9)
ggsave(filename = "dev.tests/plots/tnf.snp.mid.pdf",
       plot = plots.tnf$snp.mid, device = "pdf", width = 12, height = 9)








# Errors #

getViolinSummary <- function(ppc) {
  ppc <- ppc[ppc$parameter.level == "lowest"
             & ppc$prediction.level == "individual-level", ]
  ppc$e <- abs(ppc$y.ppc.mean-ppc$y.real.mean)

  summary.i <- aggregate(e~i+k+t+model, data = ppc, FUN = mean)
  summary.i$level <- "individual"
  summary.k <- aggregate(e~k+t+model, data = ppc, FUN = mean)
  summary.k$level <- "strain"
  summary.s <- aggregate(e~s+t+model, data = ppc, FUN = mean)
  summary.s$level <- "snp"

  summary.t <- rbind(summary.i[, c("e", "t", "model", "level")],
                     summary.k[, c("e", "t", "model", "level")],
                     summary.s[, c("e", "t", "model", "level")])
  summary.t$model <- gsub(pattern = "c", replacement = '*', x = summary.t$model)
  summary.t$level <- factor(x = summary.t$level,
                            levels = c("individual", "strain", "snp"))

  o <- aggregate(e~t+model+level, data = summary.t, FUN = mean)

  ggplot(data = summary.t)+
    facet_grid(facets = t~level, scales = "free")+
    # geom_boxplot(aes(x = model, y = e, fill = model))+
    geom_violin(aes(x = model, y = e, fill = model))+
    # geom_text(data = o, aes(x = model, y = 2,
    #                         label = round(x = e, digits = 3)),
    #           size = 3.5)+
    theme_bw()+
    theme(legend.position = "top")+
    xlab(label = '')+
    ylab(label = "Error")

}


getErrorSummary <- function(ppc, hdi.level) {

  ppc <- ppc[ppc$parameter.level == "lowest"
             & ppc$prediction.level == "individual-level", ]
  ppc$e <- abs(ppc$y.ppc.mean-ppc$y.real.mean)


  # At I
  ppc.key <- ppc[, c("i", "k", "t", "model")]
  ppc.key <- ppc.key[duplicated(ppc.key) == F, ]
  out <- c()
  for(i in 1:nrow(ppc.key)) {
    e <- ppc$e[ppc$model == ppc.key$model[i]
               & ppc$t == ppc.key$t[i]
               & ppc$k == ppc.key$k[i]
               & ppc$i == ppc.key$i[i]]

    e.hdi <- getHdi(vec = e, hdi.level = hdi.level)
    row <- data.frame(model = ppc.key$model[i],
                      t = ppc.key$t[i],
                      k = ppc.key$k[i],
                      i = ppc.key$i[i],
                      e.mean = mean(e),
                      e.L = e.hdi[1],
                      e.H = e.hdi[2])
    out <- rbind(out, row)
    cat(i, "\n")
  }
  ggplot(data = out)+
    facet_grid(facets = t~model, scales = "free_y")+
    geom_errorbar(aes(x = k, ymin = e.L, ymax = e.H),
                  width = 0.1, col = "darkgray")+
    geom_point(aes(x = k, y = e.mean, col = as.factor(k)))+
    theme_bw()+
    theme(legend.position = "top")+
    scale_x_continuous(breaks = 1:17, labels = 1:17)+
    scale_color_discrete(name = "Strain")+
    xlab(label = '')+
    ylab(label = "Mean Error (95% HDI)")+
    guides(colour = guide_legend(nrow = 2))




  # At K
  ppc.key <- ppc[, c("k", "t", "model")]
  ppc.key <- ppc.key[duplicated(ppc.key) == F, ]
  out <- c()
  for(i in 1:nrow(ppc.key)) {
    e <- ppc$e[ppc$model == ppc.key$model[i]
               & ppc$t == ppc.key$t[i]
               & ppc$k == ppc.key$k[i]]

    e.hdi <- getHdi(vec = e, hdi.level = hdi.level)
    row <- data.frame(model = ppc.key$model[i],
                      t = ppc.key$t[i],
                      k = ppc.key$k[i],
                      e.mean = mean(e),
                      e.L = e.hdi[1],
                      e.H = e.hdi[2])
    out <- rbind(out, row)
    cat(i, "\n")
  }
  ggplot(data = out)+
    facet_grid(facets = t~model, scales = "free_y")+
    geom_errorbar(aes(x = k, ymin = e.L, ymax = e.H),
                  width = 0.1, col = "darkgray")+
    geom_point(aes(x = k, y = e.mean, col = as.factor(k)))+
    theme_bw()+
    theme(legend.position = "top")+
    scale_x_continuous(breaks = 1:17, labels = 1:17)+
    scale_color_discrete(name = "Strain")+
    xlab(label = '')+
    ylab(label = "Mean Error (95% HDI)")+
    guides(colour = guide_legend(nrow = 2))




  # At T
  ppc.key <- ppc[, c("t", "model")]
  ppc.key <- ppc.key[duplicated(ppc.key) == F, ]
  out <- c()
  for(i in 1:nrow(ppc.key)) {
    e <- ppc$e[ppc$model == ppc.key$model[i]
               & ppc$t == ppc.key$t[i]]

    e.hdi <- getHdi(vec = e, hdi.level = hdi.level)
    row <- data.frame(model = ppc.key$model[i],
                      t = ppc.key$t[i],
                      e.mean = mean(e),
                      e.L = e.hdi[1],
                      e.H = e.hdi[2])
    out <- rbind(out, row)
    cat(i, "\n")
  }
  out$model <- gsub(pattern = "c", replacement = '*', x = out$model)
  out$model <- factor(x = out$model,
                      levels = c("M0", "M0*", "M1", "M1*", "M2", "M2*"))
  ggplot(data = out)+
    facet_grid(facets = ~t, scales = "free_y")+
    geom_errorbar(aes(x = model, ymin = e.L, ymax = e.H),
                  width = 0.1, col = "darkgray")+
    geom_point(aes(x = model, y = e.mean, col = as.factor(model)))+
    theme_bw()+
    theme(legend.position = "top")+
    scale_color_discrete(name = "Strain")+
    xlab(label = '')+
    ylab(label = "Mean Error (95% HDI)")+
    guides(colour = guide_legend(nrow = 2))



}





getViolinSummary(ppc = ppc.vsv)
getViolinSummary(ppc = ppc.lcmv)


