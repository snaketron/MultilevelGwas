getStrain <- function(x) {
  return (x$strain.summary)
}
getSnp <- function(x) {
  return (x$snp.summary)
}

# check loo.ic
mc$loo.ic


# check ppc
m0.snp <- do.call(rbind, lapply(mc$ppc$M0, getSnp))
m0.snp$model <- "M0"
m1.snp <- do.call(rbind, lapply(mc$ppc$M1, getSnp))
m1.snp$model <- "M1"
m2.snp <- do.call(rbind, lapply(mc$ppc$M2, getSnp))
m2.snp$model <- "M2"

m.snp <- rbind(m0.snp, m1.snp, m2.snp)

ggplot(data = m.snp)+
  facet_grid(facets = model~t)+
  geom_point(aes(y = yhat, x = y))+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()


m1.strain <- do.call(rbind, lapply(mc$ppc$M1, getStrain))
m1.strain$model <- "M1"
m2.strain <- do.call(rbind, lapply(mc$ppc$M2, getStrain))
m2.strain$model <- "M2"
m.strain <- rbind(m1.strain, m2.strain)

x <- m.strain[m.strain$s == 1 & m.strain$t == 1, ]
x <- x[complete.cases(x), ]
ggplot(data = x)+
  facet_grid(facets = model~t)+
  geom_point(aes(y = yhat, x = y))+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()


#
s0 <- data.frame(summary(mc$ps[[1]], pars = c("alpha", "beta",
                                              "sigma"))$summary)
s1 <- data.frame(summary(mc$ps[[2]], pars = c("alpha", "beta",
                                              "mu_beta", "sigma"))$summary)
s2 <- data.frame(summary(mc$ps[[3]], pars = c("alpha", "beta",
                                             "mu_beta", "sigma",
                                             "mean_trait", "sigma_trait"))$summary)

s0$par <- rownames(s0)
s1$par <- rownames(s1)
s2$par <- rownames(s2)

# "mu_beta\\[1"
s0 <- s0[regexpr(pattern = "beta\\[1", text = s0$par) != -1, ]
s1 <- s1[regexpr(pattern = "mu_beta\\[1", text = s1$par) != -1, ]
s2 <- s2[regexpr(pattern = "mu_beta\\[1", text = s2$par) != -1, ]

plot(s1$mean, s2$mean)
plot(s0$mean, s1$mean)
abline(0, 1)

s["alpha[1]", "mean"]+s["mu_beta[1,1]", "mean"]*1
s["alpha[1]", "mean"]+s["mu_beta[1,1]", "mean"]*-1

