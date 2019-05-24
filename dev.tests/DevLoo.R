# VSV
mc.vsv <- get(load("~/MLgwas/dev.tests/vsv.300.RData"))

ex.m0 <- loo::extract_log_lik(stanfit = mc.vsv$ps$M0, parameter_name = "log_lik2")
ex.m0 <- ex.m0[sample(x = 1:nrow(ex.m0), size = 1000, replace = TRUE), ]

ex.m0c <- loo::extract_log_lik(stanfit = mc.vsv$ps$M0c, parameter_name = "log_lik2")
ex.m0c <- ex.m0c[sample(x = 1:nrow(ex.m0c), size = 1000, replace = TRUE), ]

ex.m1 <- loo::extract_log_lik(stanfit = mc.vsv$ps$M1, parameter_name = "log_lik2")
ex.m1 <- ex.m1[sample(x = 1:nrow(ex.m1), size = 1000, replace = TRUE), ]

ex.m1c <- loo::extract_log_lik(stanfit = mc.vsv$ps$M1c, parameter_name = "log_lik2")
ex.m1c <- ex.m1c[sample(x = 1:nrow(ex.m1c), size = 1000, replace = TRUE), ]

rm(mc.vsv)
gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();

loo.m0 <- loo::loo(x = ex.m0)
loo.m0c <- loo::loo(x = ex.m0c)
loo.m1 <- loo::loo(x = ex.m1)
loo.m1c <- loo::loo(x = ex.m1c)

loo.m0
loo.m0c
loo.m1
loo.m1c


loo::compare(loo.m0, loo.m0c)
loo::compare(loo.m0, loo.m1)
loo::compare(loo.m0, loo.m1c)
loo::compare(loo.m0c, loo.m1)
loo::compare(loo.m0c, loo.m1c)
loo::compare(loo.m1, loo.m1c)




# LCMV
mc.lcmv <- get(load("~/MLgwas/dev.tests/lcmv.300.RData"))

ex.m0 <- loo::extract_log_lik(stanfit = mc.lcmv$ps$M0, parameter_name = "log_lik2")
ex.m0 <- ex.m0[sample(x = 1:nrow(ex.m0), size = 1000, replace = TRUE), ]

ex.m0c <- loo::extract_log_lik(stanfit = mc.lcmv$ps$M0c, parameter_name = "log_lik2")
ex.m0c <- ex.m0c[sample(x = 1:nrow(ex.m0c), size = 1000, replace = TRUE), ]

ex.m1 <- loo::extract_log_lik(stanfit = mc.lcmv$ps$M1, parameter_name = "log_lik2")
ex.m1 <- ex.m1[sample(x = 1:nrow(ex.m1), size = 1000, replace = TRUE), ]

ex.m1c <- loo::extract_log_lik(stanfit = mc.lcmv$ps$M1c, parameter_name = "log_lik2")
ex.m1c <- ex.m1c[sample(x = 1:nrow(ex.m1c), size = 1000, replace = TRUE), ]

rm(mc.lcmv)
gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();gc();

loo.m0 <- loo::loo(x = ex.m0)
loo.m0c <- loo::loo(x = ex.m0c)
loo.m1 <- loo::loo(x = ex.m1)
loo.m1c <- loo::loo(x = ex.m1c)

rm(ex.m0, ex.m0c, ex.m1, ex.m1c)
gc();gc();gc();gc();gc();gc();gc();

loo.m0
loo.m0c
loo.m1
loo.m1c


loo::compare(loo.m0, loo.m0c)
loo::compare(loo.m0, loo.m1)
loo::compare(loo.m0, loo.m1c)
loo::compare(loo.m0c, loo.m1)
loo::compare(loo.m0c, loo.m1c)
loo::compare(loo.m1, loo.m1c)





# TNF
mc.tnf <- get(load("~/MLgwas/dev.tests/tnf.300.RData"))

ex.m0 <- loo::extract_log_lik(stanfit = mc.tnf$ps$M0, parameter_name = "log_lik2")
ex.m0 <- ex.m0[sample(x = 1:nrow(ex.m0), size = 1000, replace = TRUE), ]

ex.m0c <- loo::extract_log_lik(stanfit = mc.tnf$ps$M0c, parameter_name = "log_lik2")
ex.m0c <- ex.m0c[sample(x = 1:nrow(ex.m0c), size = 1000, replace = TRUE), ]

ex.m1 <- loo::extract_log_lik(stanfit = mc.tnf$ps$M1, parameter_name = "log_lik2")
ex.m1 <- ex.m1[sample(x = 1:nrow(ex.m1), size = 1000, replace = TRUE), ]

ex.m1c <- loo::extract_log_lik(stanfit = mc.tnf$ps$M1c, parameter_name = "log_lik2")
ex.m1c <- ex.m1c[sample(x = 1:nrow(ex.m1c), size = 1000, replace = TRUE), ]

rm(mc.tnf)
gc();gc();gc();gc();gc();gc();gc();
gc();gc();gc();gc();gc();gc();gc();

loo.m0 <- loo::loo(x = ex.m0)
loo.m0c <- loo::loo(x = ex.m0c)
loo.m1 <- loo::loo(x = ex.m1)
loo.m1c <- loo::loo(x = ex.m1c)

loo.m0
loo.m0c
loo.m1
loo.m1c




loo::compare(loo.m0, loo.m0c)
loo::compare(loo.m0, loo.m1)
loo::compare(loo.m0, loo.m1c)
loo::compare(loo.m0c, loo.m1)
loo::compare(loo.m0c, loo.m1c)
loo::compare(loo.m1, loo.m1c)

