require(rstan)
require(loo)
gc();gc();gc();gc();gc();gc();gc();gc();gc();
gc();gc();gc();gc();gc();gc();gc();gc();gc();
gc();gc();gc();gc();gc();gc();gc();gc();gc();
gc();gc();gc();gc();gc();gc();gc();gc();gc();


getExtCollapsed <- function(p, limit.draws) {
  ext <- rstan::extract(object = p)
  ext <- ext$log_lik2
  dims <- dim(ext)

  # take section of the log-lik posterior
  if(limit.draws == 0) {
    limit.draws <- dims[1]
  }
  if(limit.draws > dims[1]) {
    limit.draws <- dims[1]
  }
  s <- sample(x = 1:dims[1], size = limit.draws, replace = F)

  ext.loglik <- matrix(data = 0,
                       nrow = limit.draws*dims[4],
                       ncol = dims[2]*dims[3])
  count <- 1
  for(t in 1:dims[2]) {
    for(i in 1:dims[3]) {
      ext.loglik[, count] <- as.vector(ext[s, t, i, ])
      count <- count + 1
    }
  }
  rm(t, i, count)

  return (ext.loglik)
}




# mc <- get(load("~/MLgwas/dev.tests/vsv.300.RData"))
mc <- get(load("~/MLgwas/dev.tests/lcmv.300.RData"))
# mc <- get(load("~/MLgwas/dev.tests/tnf.300.RData"))
rm(mc.vsv, mc.lcmv, mc.tnf)


m0.loo <- getExtCollapsed(p = mc$ps$M0, limit.draws = 1000)
m0c.loo <- getExtCollapsed(p = mc$ps$M0c, limit.draws = 1000)
m1.loo <- getExtCollapsed(p = mc$ps$M1, limit.draws = 1000)
m1c.loo <- getExtCollapsed(p = mc$ps$M1c, limit.draws = 1000)
rm(mc)

m0.loo <- loo::loo(x = m0.loo, cores = 2)
m0c.loo <- loo::loo(x = m0c.loo, cores = 2)
m1.loo <- loo::loo(x = m1.loo, cores = 2)
m1c.loo <- loo::loo(x = m1c.loo, cores = 2)

loo.out <- loo::compare(m0.loo, m0c.loo, m1.loo, m1c.loo)

loo.out.stats <- cbind(paste(round(x = loo.out[, 1], digits = 2), " (",
                             round(x = loo.out[, 2], digits = 2), ")", sep = ''),
                       paste(round(x = loo.out[, 3], digits = 2), " (",
                             round(x = loo.out[, 4], digits = 2), ")", sep = ''),
                       paste(round(x = loo.out[, 5], digits = 2), " (",
                             round(x = loo.out[, 6], digits = 2), ")", sep = ''))

stargazer::stargazer(loo.out.stats, summary = F)


loo.compare.df <- c()
z <- loo::compare(m0.loo, m0c.loo)
loo.compare.df <- rbind(loo.compare.df,
                        data.frame(key = "M0 vs M0*",
                                   elpd_diff = as.numeric(z[1]),
                                   se = as.numeric(z[2])))

z <- loo::compare(m0.loo, m1.loo)
loo.compare.df <- rbind(loo.compare.df,
                        data.frame(key = "M0 vs M1",
                                   elpd_diff = as.numeric(z[1]),
                                   se = as.numeric(z[2])))

z <- loo::compare(m0.loo, m1c.loo)
loo.compare.df <- rbind(loo.compare.df,
                        data.frame(key = "M0 vs M1*",
                                   elpd_diff = as.numeric(z[1]),
                                   se = as.numeric(z[2])))

z <- loo::compare(m0c.loo, m1.loo)
loo.compare.df <- rbind(loo.compare.df,
                        data.frame(key = "M0* vs M1",
                                   elpd_diff = as.numeric(z[1]),
                                   se = as.numeric(z[2])))

z <- loo::compare(m0c.loo, m1c.loo)
loo.compare.df <- rbind(loo.compare.df,
                        data.frame(key = "M0* vs M1*",
                                   elpd_diff = as.numeric(z[1]),
                                   se = as.numeric(z[2])))

z <- loo::compare(m1.loo, m1c.loo)
loo.compare.df <- rbind(loo.compare.df,
                        data.frame(key = "M1 vs M1*",
                                   elpd_diff = as.numeric(z[1]),
                                   se = as.numeric(z[2])))
rm(z)
loo.compare.df$elpd_diff <- round(x = loo.compare.df$elpd_diff, digits = 2)
loo.compare.df$se <- round(x = loo.compare.df$se, digits = 2)

