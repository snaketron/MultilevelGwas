getStrain <- function(x) {
  return (x$strain.summary)
}
getSnp <- function(x) {
  return (x$snp.summary)
}

# check loo.ic
mc$loo.ic


# ppc
m0.ppc <- do.call(rbind, mc$ppc$M0)
m0.ppc$model <- "M0"
m1.ppc <- do.call(rbind, mc$ppc$M1)
m1.ppc$model <- "M1"
m2.ppc <- do.call(rbind, mc$ppc$M2)
m2.ppc$model <- "M2"

ppc <- rbind(m0.ppc, m1.ppc, m2.ppc)

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
