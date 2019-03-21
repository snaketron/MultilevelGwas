
getPpcQ <- function(ext, gt.data, model, s, hdi.level) {


  getPpcM0Q <- function(gt.data, p, hdi.level, s) {
    getMu <- function(x, y) {
      return(mean(rnorm(n = 10, mean = x[1]+x[2]*y, sd = x[3])))
    }

    snp.summary <- c()
    for(t in 1:gt.data$Ntq) {
      alpha.p <- paste("alpha", t, sep = '.')
      sigma.p <- paste("sigma", t, sep = '.')
      beta.p <- paste("beta", t, s, sep = '.')

      if(sum(gt.data$X[, s] == 1) != 0) {
        yhat <- apply(X = p[, c(alpha.p, beta.p, sigma.p)], MARGIN = 1, FUN = getMu, y = 1)
        yhat.hdi <- getHdi(vec = yhat, hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = NA, x = 1,
                          y = mean(gt.data$Yq[gt.data$X[, s] == 1, t]),
                          yhat = mean(yhat), yhat.L = yhat.hdi[1],
                          yhat.H = yhat.hdi[2], level = "snp")
        snp.summary <- rbind(snp.summary, row)
      }
      if(sum(gt.data$X[, s] == -1) != 0) {
        yhat <- apply(X = p[, c(alpha.p, beta.p, sigma.p)], MARGIN = 1, FUN = getMu, y = -1)
        yhat.hdi <- getHdi(vec = yhat, hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = NA, x = -1,
                          y = mean(gt.data$Yq[gt.data$X[, s] == -1, t]),
                          yhat = mean(yhat), yhat.L = yhat.hdi[1],
                          yhat.H = yhat.hdi[2], level = "snp")
        snp.summary <- rbind(snp.summary, row)
      }
    }

    return (list(snp.summary = snp.summary))
  }

  getPpcM1Q <- function(gt.data, p, hdi.level, s) {

    getMu <- function(x, y) {
      return(mean(rnorm(n = 10, mean = x[1]+x[2]*y, sd = x[3])))
    }


    strain.summary <- c()
    snp.summary <- c()
    for(t in 1:gt.data$Ntq) {
      alpha.p <- paste("alpha", t, sep = '.')
      sigma.p <- paste("sigma", t, sep = '.')


      # snp-level
      yhat.snp <- c()
      xs.snp <- c()

      # strain level
      for(k in 1:gt.data$Nk) {
        beta.p <- paste("beta", t, s, k, sep = '.')
        yhat <- apply(X = p[, c(alpha.p, beta.p, sigma.p)], MARGIN = 1, FUN = getMu,
                      y = unique(gt.data$X[gt.data$K == k, s]))
        yhat.hdi <- getHdi(vec = yhat, hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = k, x = unique(gt.data$X[gt.data$K == k, s]),
                          y = mean(gt.data$Yq[gt.data$K == k, t]), yhat = mean(yhat),
                          yhat.L = yhat.hdi[1], yhat.H = yhat.hdi[2], level = "strain")
        strain.summary <- rbind(strain.summary, row)

        #
        yhat.snp <- cbind(yhat.snp, yhat)
        xs.snp <- c(xs.snp, unique(gt.data$X[gt.data$K == k, s]))
      }


      # compute snp-level (X = 1)
      if(sum(xs.snp == 1) == 0) {
        row <- data.frame(t = t, s = s, k = NA, x = 1,
                          y = NA, yhat = NA, yhat.L = NA,
                          yhat.H = NA, level = "snp")
      }
      else {
        yhat.1 <- apply(X = as.matrix(yhat.snp[, xs.snp == 1]), MARGIN = 1, FUN = mean)
        yhat.hdi1 <- getHdi(vec = yhat.1, hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = NA, x = 1,
                          y = mean(gt.data$Yq[gt.data$X[, s] == 1, t]),
                          yhat = mean(yhat.1), yhat.L = yhat.hdi1[1],
                          yhat.H = yhat.hdi1[2], level = "snp")
      }
      snp.summary <- rbind(snp.summary, row)

      # compute snp-level (X = -1)
      if(sum(xs.snp == -1) == 0) {
        row <- data.frame(t = t, s = s, k = NA, x = -1,
                          y = NA, yhat = NA, yhat.L = NA,
                          yhat.H = NA, level = "snp")
      }
      else {
        yhat.0 <- apply(X = as.matrix(yhat.snp[, xs.snp == -1]), MARGIN = 1, FUN = mean)
        yhat.hdi0 <- getHdi(vec = yhat.0, hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = NA, x = -1,
                          y = mean(gt.data$Yq[gt.data$X[, s] == -1, t]),
                          yhat = mean(yhat.0), yhat.L = yhat.hdi0[1],
                          yhat.H = yhat.hdi0[2], level = "snp")
      }
      snp.summary <- rbind(snp.summary, row)


      # remove
      rm(yhat.1, yhat.hdi1, k, xs.snp, yhat, yhat.snp,
         yhat.0, yhat.hdi0, row, beta.p, sigma.p, y)
    }

    return (list(snp.summary = snp.summary, strain.summary = strain.summary))
  }

  getPpcM2Q <- function(gt.data, p, hdi.level, s) {
    getMu <- function(x, y) {
      return(mean(rnorm(n = 10, mean = x[1]+x[2]*y, sd = x[3])))
    }


    strain.summary <- c()
    snp.summary <- c()
    for(t in 1:gt.data$Ntq) {
      alpha.p <- paste("alpha", t, sep = '.')


      # snp-level
      yhat.snp <- c()
      xs.snp <- c()

      # strain level
      for(k in 1:gt.data$Nk) {
        beta.p <- paste("beta", t, s, k, sep = '.')
        sigma.p <- paste("sigma", k, sep = '.')
        yhat <- apply(X = p[, c(alpha.p, beta.p, sigma.p)], MARGIN = 1, FUN = getMu,
                      y = unique(gt.data$X[gt.data$K == k, s]))
        yhat.hdi <- getHdi(vec = yhat, hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = k, x = unique(gt.data$X[gt.data$K == k, s]),
                          y = mean(gt.data$Yq[gt.data$K == k, t]), yhat = mean(yhat),
                          yhat.L = yhat.hdi[1], yhat.H = yhat.hdi[2], level = "strain")
        strain.summary <- rbind(strain.summary, row)

        #
        yhat.snp <- cbind(yhat.snp, yhat)
        xs.snp <- c(xs.snp, unique(gt.data$X[gt.data$K == k, s]))
      }



      # compute snp-level (X = 1)
      if(sum(xs.snp == 1) == 0) {
        row <- data.frame(t = t, s = s, k = NA, x = 1,
                          y = NA, yhat = NA, yhat.L = NA,
                          yhat.H = NA, level = "snp")
      }
      else {
        yhat.1 <- apply(X = as.matrix(yhat.snp[, xs.snp == 1]), MARGIN = 1, FUN = mean)
        yhat.hdi1 <- getHdi(vec = yhat.1, hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = NA, x = 1,
                          y = mean(gt.data$Yq[gt.data$X[, s] == 1, t]),
                          yhat = mean(yhat.1), yhat.L = yhat.hdi1[1],
                          yhat.H = yhat.hdi1[2], level = "snp")
      }
      snp.summary <- rbind(snp.summary, row)

      # compute snp-level (X = -1)
      if(sum(xs.snp == -1) == 0) {
        row <- data.frame(t = t, s = s, k = NA, x = -1,
                          y = NA, yhat = NA, yhat.L = NA,
                          yhat.H = NA, level = "snp")
      }
      else {
        yhat.0 <- apply(X = as.matrix(yhat.snp[, xs.snp == -1]), MARGIN = 1, FUN = mean)
        yhat.hdi0 <- getHdi(vec = yhat.0, hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = NA, x = -1,
                          y = mean(gt.data$Yq[gt.data$X[, s] == -1, t]),
                          yhat = mean(yhat.0), yhat.L = yhat.hdi0[1],
                          yhat.H = yhat.hdi0[2], level = "snp")
      }
      snp.summary <- rbind(snp.summary, row)


      # remove
      rm(yhat.1, yhat.hdi1, k, xs.snp, yhat, yhat.snp,
         yhat.0, yhat.hdi0, row, beta.p, sigma.p, y)
    }

    return (list(snp.summary = snp.summary, strain.summary = strain.summary))
  }


  if(model == "M0") {
    return(getPpcM0Q(gt.data = gt.data, p = ext,
                     hdi.level = hdi.level, s = s))
  }
  else if(model == "M1") {
    return(getPpcM1Q(gt.data = gt.data, p = ext,
                     hdi.level = hdi.level, s = s))
  }
  else if(model == "M2") {
    return(getPpcM2Q(gt.data = gt.data, p = ext,
                     hdi.level = hdi.level, s = s))
  }

  return (ppc)
}



getPpcD <- function(ext, gt.data, model, s, hdi.level) {


  getPpcM0D <- function(gt.data, p, hdi.level, s) {
    getPr <- function(x, y) {
      1/(1 + exp(-(x[1]+x[2]*y)))
      return(mean(rbinom(n = 10, size = 1, prob = 1/(1 + exp(-(x[1]+x[2]*y))))))
    }

    snp.summary <- c()

    # D-trait
    for(d in 1:gt.data$Ntd) {
      t <- d + gt.data$Ntq
      alpha.p <- paste("alpha", t, sep = '.')
      sigma.p <- paste("sigma", t, sep = '.')
      beta.p <- paste("beta", t, s, sep = '.')

      if(sum(gt.data$X[, s] == 1) != 0) {
        yhat <- apply(X = p[, c(alpha.p, beta.p)], MARGIN = 1, FUN = getPr, y = 1)
        yhat.hdi <- getHdi(vec = yhat, hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = NA, x = 1,
                          y = mean(as.numeric(gt.data$Yd[gt.data$X[, s] == 1, d])),
                          yhat = mean(yhat), yhat.L = yhat.hdi[1],
                          yhat.H = yhat.hdi[2], level = "snp")
        snp.summary <- rbind(snp.summary, row)
      }
      if(sum(gt.data$X[, s] == -1) != 0) {
        yhat <- apply(X = p[, c(alpha.p, beta.p)], MARGIN = 1, FUN = getPr, y = -1)
        yhat.hdi <- getHdi(vec = yhat, hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = NA, x = -1,
                          y = mean(as.numeric(gt.data$Yd[gt.data$X[, s] == -1, d])),
                          yhat = mean(yhat), yhat.L = yhat.hdi[1],
                          yhat.H = yhat.hdi[2], level = "snp")
        snp.summary <- rbind(snp.summary, row)
      }
    }

    return (list(snp.summary = snp.summary))
  }

  getPpcM1D <- function(gt.data, p, hdi.level, s) {

    getPr <- function(x, y) {
      1/(1 + exp(-(x[1]+x[2]*y)))
      return(mean(rbinom(n = 10, size = 1, prob = 1/(1 + exp(-(x[1]+x[2]*y))))))
    }

    strain.summary <- c()
    snp.summary <- c()

    # D-trait
    for(d in 1:gt.data$Ntd) {
      t <- d + gt.data$Ntq

      alpha.p <- paste("alpha", t, sep = '.')

      # snp-level
      yhat.snp <- c()
      xs.snp <- c()

      # strain level
      for(k in 1:gt.data$Nk) {
        beta.p <- paste("beta", t, s, k, sep = '.')
        yhat <- apply(X = p[, c(alpha.p, beta.p)], MARGIN = 1, FUN = getPr,
                      y = unique(gt.data$X[gt.data$K == k, s]))
        yhat.hdi <- getHdi(vec = yhat, hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = k, x = unique(gt.data$X[gt.data$K == k, s]),
                          y = mean(gt.data$Yd[gt.data$K == k, d]), yhat = mean(yhat),
                          yhat.L = yhat.hdi[1], yhat.H = yhat.hdi[2], level = "strain")
        strain.summary <- rbind(strain.summary, row)

        #
        yhat.snp <- cbind(yhat.snp, yhat)
        xs.snp <- c(xs.snp, unique(gt.data$X[gt.data$K == k, s]))
      }


      # compute snp-level (X = 1)
      if(sum(xs.snp == 1) == 0) {
        row <- data.frame(t = t, s = s, k = NA, x = 1,
                          y = NA, yhat = NA, yhat.L = NA,
                          yhat.H = NA, level = "snp")
      }
      else {
        yhat.1 <- apply(X = as.matrix(yhat.snp[, xs.snp == 1]), MARGIN = 1, FUN = mean)
        yhat.hdi1 <- getHdi(vec = yhat.1, hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = NA, x = 1,
                          y = mean(as.numeric(gt.data$Yd[gt.data$X[, s] == 1, d])),
                          yhat = mean(yhat.1), yhat.L = yhat.hdi1[1],
                          yhat.H = yhat.hdi1[2], level = "snp")
      }
      snp.summary <- rbind(snp.summary, row)

      # compute snp-level (X = -1)
      if(sum(xs.snp == -1) == 0) {
        row <- data.frame(t = t, s = s, k = NA, x = -1,
                          y = NA, yhat = NA, yhat.L = NA,
                          yhat.H = NA, level = "snp")
      }
      else {
        yhat.0 <- apply(X = as.matrix(yhat.snp[, xs.snp == -1]), MARGIN = 1, FUN = mean)
        yhat.hdi0 <- getHdi(vec = yhat.0, hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = NA, x = -1,
                          y = mean(as.numeric(gt.data$Yd[gt.data$X[, s] == -1, d])),
                          yhat = mean(yhat.0), yhat.L = yhat.hdi0[1],
                          yhat.H = yhat.hdi0[2], level = "snp")
      }
      snp.summary <- rbind(snp.summary, row)


      # remove
      rm(yhat.1, yhat.hdi1, k, xs.snp, yhat, yhat.snp,
         yhat.0, yhat.hdi0, row, beta.p, y)

    }

    return (list(snp.summary = snp.summary, strain.summary = strain.summary))
  }

  getPpcM2D <- function(gt.data, p, hdi.level, s) {

    getPr <- function(x, y) {
      1/(1 + exp(-(x[1]+x[2]*y)))
      return(mean(rbinom(n = 10, size = 1, prob = 1/(1 + exp(-(x[1]+x[2]*y))))))
    }


    strain.summary <- c()
    snp.summary <- c()

    # D-trait
    for(d in 1:gt.data$Ntd) {
      t <- d + gt.data$Ntq

      alpha.p <- paste("alpha", t, sep = '.')

      # snp-level
      yhat.snp <- c()
      xs.snp <- c()

      # strain level
      for(k in 1:gt.data$Nk) {
        beta.p <- paste("beta", t, s, k, sep = '.')
        yhat <- apply(X = p[, c(alpha.p, beta.p)], MARGIN = 1, FUN = getPr,
                      y = unique(gt.data$X[gt.data$K == k, s]))
        yhat.hdi <- getHdi(vec = yhat, hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = k, x = unique(gt.data$X[gt.data$K == k, s]),
                          y = mean(gt.data$Yd[gt.data$K == k, d]), yhat = mean(yhat),
                          yhat.L = yhat.hdi[1], yhat.H = yhat.hdi[2], level = "strain")
        strain.summary <- rbind(strain.summary, row)

        #
        yhat.snp <- cbind(yhat.snp, yhat)
        xs.snp <- c(xs.snp, unique(gt.data$X[gt.data$K == k, s]))
      }


      # compute snp-level (X = 1)
      if(sum(xs.snp == 1) == 0) {
        row <- data.frame(t = t, s = s, k = NA, x = 1,
                          y = NA, yhat = NA, yhat.L = NA,
                          yhat.H = NA, level = "snp")
      }
      else {
        yhat.1 <- apply(X = as.matrix(yhat.snp[, xs.snp == 1]), MARGIN = 1, FUN = mean)
        yhat.hdi1 <- getHdi(vec = yhat.1, hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = NA, x = 1,
                          y = mean(as.numeric(gt.data$Yd[gt.data$X[, s] == 1, d])),
                          yhat = mean(yhat.1), yhat.L = yhat.hdi1[1],
                          yhat.H = yhat.hdi1[2], level = "snp")
      }
      snp.summary <- rbind(snp.summary, row)

      # compute snp-level (X = -1)
      if(sum(xs.snp == -1) == 0) {
        row <- data.frame(t = t, s = s, k = NA, x = -1,
                          y = NA, yhat = NA, yhat.L = NA,
                          yhat.H = NA, level = "snp")
      }
      else {
        yhat.0 <- apply(X = as.matrix(yhat.snp[, xs.snp == -1]), MARGIN = 1, FUN = mean)
        yhat.hdi0 <- getHdi(vec = yhat.0, hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = NA, x = -1,
                          y = mean(as.numeric(gt.data$Yd[gt.data$X[, s] == -1, d])),
                          yhat = mean(yhat.0), yhat.L = yhat.hdi0[1],
                          yhat.H = yhat.hdi0[2], level = "snp")
      }
      snp.summary <- rbind(snp.summary, row)


      # remove
      rm(yhat.1, yhat.hdi1, k, xs.snp, yhat, yhat.snp,
         yhat.0, yhat.hdi0, row, beta.p, y)
    }


    return (list(snp.summary = snp.summary,
                 strain.summary = strain.summary))
  }

  if(model == "M0") {

  }
  else if(model == "M1") {

  }
  else if(model == "M2") {

  }
}



getPpcQD <- function(ext, gt.data, model, s, hdi.level) {

  getMu <- function(x, y) {
    return(rnorm(n = 1, mean = x[1]+x[2]*y, sd = x[3]))
  }

  getPr <- function(x, y) {
    return(1/(1 + exp(-(x[1]+x[2]*y))))
  }


  getPpcM0QD <- function(gt.data, p, hdi.level, s) {
    ppc.summary <- c()

    # Q-trait
    for(t in 1:gt.data$Ntq) {
      alpha.p <- paste("alpha", t, sep = '.')
      sigma.p <- paste("sigma", t, sep = '.')
      # in case single Ntq => sigma no index
      if(gt.data$Ntq == 1) {
        sigma.p <- "sigma"
      }

      # individual level
      yhat <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)
      for(i in 1:gt.data$N) {
        beta.p <- paste("beta", t, s, sep = '.')

        yhat[, i] <- apply(X = p[, c(alpha.p, beta.p, sigma.p)], MARGIN = 1,
                           FUN = getMu, y = gt.data$X[i, s])
      }

      # snp level
      if(sum(gt.data$X[, s] == 1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == 1]),
                           hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = NA, x = 1, level = "snp",
                          y.real = mean(gt.data$Yq[gt.data$X[, s] == 1, t]),
                          y.ppc = mean(yhat[, gt.data$X[, s] == 1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
      # snp level
      if(sum(gt.data$X[, s] == -1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == -1]),
                           hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = NA, x = -1, level = "snp",
                          y.real = mean(gt.data$Yq[gt.data$X[, s] == -1, t]),
                          y.ppc = mean(yhat[, gt.data$X[, s] == -1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
    }
    # clean up
    rm(yhat, yhat.hdi, row)


    # D-trait
    for(d in 1:gt.data$Ntd) {
      t <- d + gt.data$Ntq

      alpha.p <- paste("alpha", t, sep = '.')

      # individual level
      yhat <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)
      for(i in 1:gt.data$N) {
        beta.p <- paste("beta", t, s, sep = '.')

        yhat[, i] <- apply(X = p[, c(alpha.p, beta.p)], MARGIN = 1,
                           FUN = getPr, y = gt.data$X[i, s])
      }

      # snp level
      if(sum(gt.data$X[, s] == 1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == 1]),
                           hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = NA, x = 1, level = "snp",
                          y.real=mean(as.numeric(gt.data$Yd[gt.data$X[, s] == 1,d])),
                          y.ppc = mean(yhat[, gt.data$X[, s] == 1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
      # snp level
      if(sum(gt.data$X[, s] == -1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == -1]),
                           hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = NA, x = -1, level = "snp",
                          y.real=mean(as.numeric(gt.data$Yd[gt.data$X[, s] == -1,d])),
                          y.ppc = mean(yhat[, gt.data$X[, s] == -1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
    }
    # clean up
    rm(i, yhat, yhat.hdi, row)

    return (ppc.summary)
  }

  getPpcM1QD <- function(gt.data, p, hdi.level, s) {
    ppc.summary <- c()

    # Q-trait
    for(t in 1:gt.data$Ntq) {
      alpha.p <- paste("alpha", t, sep = '.')
      sigma.p <- paste("sigma", t, sep = '.')
      # in case single Ntq => sigma no index
      if(gt.data$Ntq == 1) {
        sigma.p <- "sigma"
      }

      # individual level
      yhat <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)
      for(i in 1:gt.data$N) {
        beta.p <- paste("beta", t, s, gt.data$K[i], sep = '.')
        # TODO: check in case single Ntq => sigma no index
        if(gt.data$Nk == 1) {
          beta.p <- paste("beta", t, s, sep = '.')
        }

        yhat[, i] <- apply(X = p[, c(alpha.p, beta.p, sigma.p)], MARGIN = 1,
                           FUN = getMu, y = gt.data$X[i, s])
      }

      # strain level
      for(k in 1:gt.data$Nk) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$K == k]),
                           hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = k, x = NA, level = "strain",
                          y.real = mean(gt.data$Yq[gt.data$K == k, t]),
                          y.ppc = mean(yhat[, gt.data$K == k]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }

      # snp level
      if(sum(gt.data$X[, s] == 1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == 1]),
                           hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = NA, x = 1, level = "snp",
                          y.real = mean(gt.data$Yq[gt.data$X[, s] == 1, t]),
                          y.ppc = mean(yhat[, gt.data$X[, s] == 1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
      # snp level
      if(sum(gt.data$X[, s] == -1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == -1]),
                           hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = NA, x = -1, level = "snp",
                          y.real = mean(gt.data$Yq[gt.data$X[, s] == -1, t]),
                          y.ppc = mean(yhat[, gt.data$X[, s] == -1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
    }
    # clean up
    rm(k, yhat, yhat.hdi, row)


    # D-trait
    for(d in 1:gt.data$Ntd) {
      t <- d + gt.data$Ntq

      alpha.p <- paste("alpha", t, sep = '.')

      # individual level
      yhat <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)
      for(i in 1:gt.data$N) {
        beta.p <- paste("beta", t, s, gt.data$K[i], sep = '.')
        # TODO: check in case single Ntq => sigma no index
        if(gt.data$Nk == 1) {
          beta.p <- paste("beta", t, s, sep = '.')
        }

        yhat[, i] <- apply(X = p[, c(alpha.p, beta.p)], MARGIN = 1,
                           FUN = getPr, y = gt.data$X[i, s])
      }

      # strain level
      for(k in 1:gt.data$Nk) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$K == k]),
                           hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = k, x = NA, level = "strain",
                          y.real=mean(as.numeric(gt.data$Yq[gt.data$K == k,d])),
                          y.ppc = mean(yhat[, gt.data$K == k]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }

      # snp level
      if(sum(gt.data$X[, s] == 1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == 1]),
                           hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = NA, x = 1, level = "snp",
                          y.real=mean(as.numeric(gt.data$Yd[gt.data$X[, s]==1,d])),
                          y.ppc = mean(yhat[, gt.data$X[, s] == 1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
      # snp level
      if(sum(gt.data$X[, s] == -1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == -1]),
                           hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = NA, x = -1, level = "snp",
                          y.real=mean(as.numeric(gt.data$Yd[gt.data$X[, s]==-1,d])),
                          y.ppc = mean(yhat[, gt.data$X[, s] == -1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
    }
    # clean up
    rm(k, yhat, yhat.hdi, row)

    return (ppc.summary)
  }

  getPpcM2QD <- function(gt.data, p, hdi.level, s) {
    ppc.summary <- c()

    # Q-trait
    for(t in 1:gt.data$Ntq) {
      alpha.p <- paste("alpha", t, sep = '.')

      # individual level
      yhat <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)
      for(i in 1:gt.data$N) {
        beta.p <- paste("beta", t, s, gt.data$K[i], sep = '.')
        sigma.p <- paste("sigma", gt.data$K[i], sep = '.')
        # TODO: check in case single Ntq => sigma no index
        if(gt.data$Nk == 1) {
          beta.p <- paste("beta", t, s, sep = '.')
          sigma.p <- "sigma"
        }

        yhat[, i] <- apply(X = p[, c(alpha.p, beta.p, sigma.p)], MARGIN = 1,
                           FUN = getMu, y = gt.data$X[i, s])
      }

      # strain level
      for(k in 1:gt.data$Nk) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$K == k]),
                           hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = k, x = NA, level = "strain",
                          y.real = mean(gt.data$Yq[gt.data$K == k, t]),
                          y.ppc = mean(yhat[, gt.data$K == k]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }

      # snp level
      if(sum(gt.data$X[, s] == 1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == 1]),
                           hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = NA, x = 1, level = "snp",
                          y.real = mean(gt.data$Yq[gt.data$X[, s] == 1, t]),
                          y.ppc = mean(yhat[, gt.data$X[, s] == 1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
      # snp level
      if(sum(gt.data$X[, s] == -1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == -1]),
                           hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = NA, x = -1, level = "snp",
                          y.real = mean(gt.data$Yq[gt.data$X[, s] == -1, t]),
                          y.ppc = mean(yhat[, gt.data$X[, s] == -1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
    }
    # clean up
    rm(k, yhat, yhat.hdi, row)


    # D-trait
    for(d in 1:gt.data$Ntd) {
      t <- d + gt.data$Ntq

      alpha.p <- paste("alpha", t, sep = '.')

      # individual level
      yhat <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)
      for(i in 1:gt.data$N) {
        beta.p <- paste("beta", t, s, gt.data$K[i], sep = '.')
        # TODO: check in case single Ntq => sigma no index
        if(gt.data$Nk == 1) {
          beta.p <- paste("beta", t, s, sep = '.')
        }

        yhat[, i] <- apply(X = p[, c(alpha.p, beta.p)], MARGIN = 1,
                           FUN = getPr, y = gt.data$X[i, s])
      }

      # strain level
      for(k in 1:gt.data$Nk) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$K == k]),
                           hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = k, x = NA, level = "strain",
                          y.real=mean(as.numeric(gt.data$Yq[gt.data$K == k,d])),
                          y.ppc = mean(yhat[, gt.data$K == k]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }

      # snp level
      if(sum(gt.data$X[, s] == 1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == 1]),
                           hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = NA, x = 1, level = "snp",
                          y.real=mean(as.numeric(gt.data$Yd[gt.data$X[, s] == 1,d])),
                          y.ppc = mean(yhat[, gt.data$X[, s] == 1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
      # snp level
      if(sum(gt.data$X[, s] == -1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == -1]),
                           hdi.level = hdi.level)
        row <- data.frame(t = t, s = s, k = NA, x = -1, level = "snp",
                          y.real=mean(as.numeric(gt.data$Yd[gt.data$X[, s] == -1,d])),
                          y.ppc = mean(yhat[, gt.data$X[, s] == -1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
    }
    # clean up
    rm(k, yhat, yhat.hdi, row)

    return (ppc.summary)
  }

  if(model == "M0") {
    return(getPpcM0QD(gt.data = gt.data, p = ext,
                      hdi.level = hdi.level, s = s))
  }
  else if(model == "M1") {
    return(getPpcM1QD(gt.data = gt.data, p = ext,
                      hdi.level = hdi.level, s = s))
  }
  else if(model == "M2") {
    return(getPpcM2QD(gt.data = gt.data, p = ext,
                      hdi.level = hdi.level, s = s))
  }
}

