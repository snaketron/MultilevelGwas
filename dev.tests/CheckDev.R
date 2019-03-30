ppc <- getPpc(ps = mc$ps, gt.data = mc$gt.data,
              models = names(mc$ps), hdi.level = 0.95)
mc$ppc <- ppc
save(mc, file = "dev.tests/mc.tnf.cov.300.RData")

# check loo.ic
mc$loo.ic


mc$ppc$M0$model <- "M0"
mc$ppc$M1$model <- "M1"
mc$ppc$M2$model <- "M2"

mc$ppc$M0c$model <- "M0c"
mc$ppc$M1c$model <- "M1c"
mc$ppc$M2c$model <- "M2c"

ppc <- rbind(mc$ppc$M0, mc$ppc$M1, mc$ppc$M2)

ppc <- rbind(mc$ppc$M0c, mc$ppc$M1c, mc$ppc$M2c)

ggplot(data = ppc[ppc$level == "snp", ])+
  facet_grid(facets = model~t)+
  geom_errorbar(aes(ymin = y.ppc.L, ymax = y.ppc.H, x = y.real), col = "darkgray")+
  geom_point(aes(y = y.ppc, x = y.real))+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()


ggplot(data = ppc[ppc$level == "snp", ])+
  facet_wrap(facets = ~t)+
  geom_density(aes(y.real-y.ppc, col = model))+
  theme_bw()

ggplot(data = ppc[ppc$level == "snp", ])+
  facet_wrap(facets = ~t, scales = "free_x")+
  geom_density(aes(y.ppc.H-y.ppc.L, col = model))+
  theme_bw()


ggplot(data = ppc[ppc$level == "strain", ])+
  facet_grid(facets = model~t)+
  geom_errorbar(aes(ymin = y.ppc.L, ymax = y.ppc.H, x = y.real), col = "darkgray")+
  geom_point(aes(y = y.ppc, x = y.real))+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()
