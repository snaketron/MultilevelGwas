
# check loo.ic
mc$loo.ic

mc$ppc$M0$model <- "M0"
mc$ppc$M1$model <- "M1"
mc$ppc$M2$model <- "M2"
ppc <- rbind(mc$ppc$M0, mc$ppc$M1, mc$ppc$M2)

ggplot(data = ppc[ppc$level == "snp", ])+
  facet_grid(facets = model~t)+
  geom_errorbar(aes(ymin = y.ppc.L, ymax = y.ppc.H, x = y.real), col = "darkgray")+
  geom_point(aes(y = y.ppc, x = y.real))+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()

ggplot(data = ppc[ppc$level == "strain", ])+
  facet_grid(facets = model~t)+
  geom_errorbar(aes(ymin = y.ppc.L, ymax = y.ppc.H, x = y.real), col = "darkgray")+
  geom_point(aes(y = y.ppc, x = y.real))+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()
