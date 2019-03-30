

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
      errors <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)

      for(i in 1:gt.data$N) {
        beta.p <- paste("beta", t, s, sep = '.')

        yhat[, i] <- apply(X = p[, c(alpha.p, beta.p, sigma.p)], MARGIN = 1,
                           FUN = getMu, y = gt.data$X[i, s])
        errors[, i] <- abs(gt.data$Yq[i, t] - yhat[, i])
      }

      # snp level
      if(sum(gt.data$X[, s] == 1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == 1]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$X[, s] == 1]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = NA, x = 1, level = "snp",
                          y.real = mean(gt.data$Yq[gt.data$X[, s] == 1, t]),
                          y.ppc = mean(yhat[, gt.data$X[, s] == 1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$X[, s] == 1]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)

      }
      # snp level
      if(sum(gt.data$X[, s] == -1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == -1]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$X[, s] == -1]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = NA, x = -1, level = "snp",
                          y.real = mean(gt.data$Yq[gt.data$X[, s] == -1, t]),
                          y.ppc = mean(yhat[, gt.data$X[, s] == -1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$X[, s] == -1]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
    }
    # clean up
    rm(yhat, errors, yhat.hdi, row)


    # D-trait
    for(d in 1:gt.data$Ntd) {
      t <- d + gt.data$Ntq

      alpha.p <- paste("alpha", t, sep = '.')

      # individual level
      yhat <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)
      errors <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)

      for(i in 1:gt.data$N) {
        beta.p <- paste("beta", t, s, sep = '.')

        yhat[, i] <- apply(X = p[, c(alpha.p, beta.p)], MARGIN = 1,
                           FUN = getPr, y = gt.data$X[i, s])
        errors[, i] <- abs(gt.data$Yd[i, d] - yhat[, i])
      }

      # snp level
      if(sum(gt.data$X[, s] == 1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == 1]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$X[, s] == 1]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = NA, x = 1, level = "snp",
                          y.real=mean(as.numeric(gt.data$Yd[gt.data$X[,s]==1,d])),
                          y.ppc = mean(yhat[, gt.data$X[, s] == 1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$X[, s] == 1]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
      # snp level
      if(sum(gt.data$X[, s] == -1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == -1]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$X[, s] == -1]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = NA, x = -1, level = "snp",
                          y.real=mean(as.numeric(gt.data$Yd[gt.data$X[,s]==-1,d])),
                          y.ppc = mean(yhat[, gt.data$X[, s] == -1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$X[, s] == -1]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
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
      errors <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)

      for(i in 1:gt.data$N) {
        beta.p <- paste("beta", t, s, gt.data$K[i], sep = '.')
        # TODO: check in case single Ntq => sigma no index
        if(gt.data$Nk == 1) {
          beta.p <- paste("beta", t, s, sep = '.')
        }

        yhat[, i] <- apply(X = p[, c(alpha.p, beta.p, sigma.p)], MARGIN = 1,
                           FUN = getMu, y = gt.data$X[i, s])
        errors[, i] <- abs(gt.data$Yq[i, t] - yhat[, i])
      }

      # strain level
      for(k in 1:gt.data$Nk) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$K == k]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$K == k]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = k, x = NA, level = "strain",
                          y.real = mean(gt.data$Yq[gt.data$K == k, t]),
                          y.ppc = mean(yhat[, gt.data$K == k]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$K == k]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }

      # snp level
      if(sum(gt.data$X[, s] == 1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == 1]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$X[, s] == 1]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = NA, x = 1, level = "snp",
                          y.real = mean(gt.data$Yq[gt.data$X[, s] == 1, t]),
                          y.ppc = mean(yhat[, gt.data$X[, s] == 1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$X[, s] == 1]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
      # snp level
      if(sum(gt.data$X[, s] == -1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == -1]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$X[, s] == -1]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = NA, x = -1, level = "snp",
                          y.real = mean(gt.data$Yq[gt.data$X[, s] == -1, t]),
                          y.ppc = mean(yhat[, gt.data$X[, s] == -1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$X[, s] == -1]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
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
      errors <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)

      for(i in 1:gt.data$N) {
        beta.p <- paste("beta", t, s, gt.data$K[i], sep = '.')
        # TODO: check in case single Ntq => sigma no index
        if(gt.data$Nk == 1) {
          beta.p <- paste("beta", t, s, sep = '.')
        }

        yhat[, i] <- apply(X = p[, c(alpha.p, beta.p)], MARGIN = 1,
                           FUN = getPr, y = gt.data$X[i, s])
        errors[, i] <- abs(gt.data$Yd[i, t] - yhat[, i])
      }

      # strain level
      for(k in 1:gt.data$Nk) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$K == k]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$K == k]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = k, x = NA, level = "strain",
                          y.real=mean(as.numeric(gt.data$Yd[gt.data$K == k,d])),
                          y.ppc = mean(yhat[, gt.data$K == k]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$K == k]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }

      # snp level
      if(sum(gt.data$X[, s] == 1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == 1]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$X[, s] == 1]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = NA, x = 1, level = "snp",
                          y.real=mean(as.numeric(gt.data$Yd[gt.data$X[, s]==1,d])),
                          y.ppc = mean(yhat[, gt.data$X[, s] == 1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$X[, s] == 1]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
      # snp level
      if(sum(gt.data$X[, s] == -1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == -1]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$X[, s] == -1]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = NA, x = -1, level = "snp",
                          y.real=mean(as.numeric(gt.data$Yd[gt.data$X[,s]==-1,d])),
                          y.ppc = mean(yhat[, gt.data$X[, s] == -1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$X[, s] == -1]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
    }
    # clean up
    rm(k, errors, yhat, yhat.hdi, row)

    return (ppc.summary)
  }

  getPpcM2QD <- function(gt.data, p, hdi.level, s) {
    ppc.summary <- c()

    # Q-trait
    for(t in 1:gt.data$Ntq) {
      alpha.p <- paste("alpha", t, sep = '.')

      # individual level
      yhat <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)
      errors <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)

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
        errors[, i] <- abs(gt.data$Yq[i, t] - yhat[, i])
      }

      # strain level
      for(k in 1:gt.data$Nk) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$K == k]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$K == k]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = k, x = NA, level = "strain",
                          y.real = mean(gt.data$Yq[gt.data$K == k, t]),
                          y.ppc = mean(yhat[, gt.data$K == k]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$K == k]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }

      # snp level
      if(sum(gt.data$X[, s] == 1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == 1]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$X[, s] == 1]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = NA, x = 1, level = "snp",
                          y.real = mean(gt.data$Yq[gt.data$X[, s] == 1, t]),
                          y.ppc = mean(yhat[, gt.data$X[, s] == 1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$X[, s] == 1]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
      # snp level
      if(sum(gt.data$X[, s] == -1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == -1]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$X[, s] == -1]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = NA, x = -1, level = "snp",
                          y.real = mean(gt.data$Yq[gt.data$X[, s] == -1, t]),
                          y.ppc = mean(yhat[, gt.data$X[, s] == -1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$X[, s] == -1]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
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
      errors <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)

      for(i in 1:gt.data$N) {
        beta.p <- paste("beta", t, s, gt.data$K[i], sep = '.')
        # TODO: check in case single Ntq => sigma no index
        if(gt.data$Nk == 1) {
          beta.p <- paste("beta", t, s, sep = '.')
        }

        yhat[, i] <- apply(X = p[, c(alpha.p, beta.p)], MARGIN = 1,
                           FUN = getPr, y = gt.data$X[i, s])
        errors[, i] <- abs(gt.data$Yd[i, d] - yhat[, i])
      }

      # strain level
      for(k in 1:gt.data$Nk) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$K == k]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$K == k]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = k, x = NA, level = "strain",
                          y.real=mean(as.numeric(gt.data$Yd[gt.data$K == k,d])),
                          y.ppc = mean(yhat[, gt.data$K == k]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$K == k]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }

      # snp level
      if(sum(gt.data$X[, s] == 1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == 1]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$X[, s] == 1]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = NA, x = 1, level = "snp",
                          y.real=mean(as.numeric(gt.data$Yd[gt.data$X[, s] == 1,d])),
                          y.ppc = mean(yhat[, gt.data$X[, s] == 1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$X[, s] == 1]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
      # snp level
      if(sum(gt.data$X[, s] == -1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == -1]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$X[, s] == -1]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = NA, x = -1, level = "snp",
                          y.real=mean(as.numeric(gt.data$Yd[gt.data$X[, s] == -1,d])),
                          y.ppc = mean(yhat[, gt.data$X[, s] == -1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$X[, s] == -1]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
    }
    # clean up
    rm(k, yhat, yhat.hdi, row)

    return (ppc.summary)
  }

  if(model == "M0" | model == "M0c") {
    return(getPpcM0QD(gt.data = gt.data, p = ext,
                      hdi.level = hdi.level, s = s))
  }
  else if(model == "M1" | model == "M1c") {
    return(getPpcM1QD(gt.data = gt.data, p = ext,
                      hdi.level = hdi.level, s = s))
  }
  else if(model == "M2" | model == "M2c") {
    return(getPpcM2QD(gt.data = gt.data, p = ext,
                      hdi.level = hdi.level, s = s))
  }
}



getPpcQ <- function(ext, gt.data, model, s, hdi.level) {

  getMu <- function(x, y) {
    return(rnorm(n = 1, mean = x[1]+x[2]*y, sd = x[3]))
  }

  getPpcM0Q <- function(gt.data, p, hdi.level, s) {
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
      errors <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)

      for(i in 1:gt.data$N) {
        beta.p <- paste("beta", t, s, sep = '.')

        yhat[, i] <- apply(X = p[, c(alpha.p, beta.p, sigma.p)], MARGIN = 1,
                           FUN = getMu, y = gt.data$X[i, s])
        errors[, i] <- abs(gt.data$Yq[i, t] - yhat[, i])
      }

      # snp level
      if(sum(gt.data$X[, s] == 1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == 1]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$X[, s] == 1]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = NA, x = 1, level = "snp",
                          y.real = mean(gt.data$Yq[gt.data$X[, s] == 1, t]),
                          y.ppc = mean(yhat[, gt.data$X[, s] == 1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$X[, s] == 1]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
      # snp level
      if(sum(gt.data$X[, s] == -1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == -1]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$X[, s] == -1]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = NA, x = -1, level = "snp",
                          y.real = mean(gt.data$Yq[gt.data$X[, s] == -1, t]),
                          y.ppc = mean(yhat[, gt.data$X[, s] == -1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$X[, s] == -1]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
    }
    # clean up
    rm(yhat, yhat.hdi, row)

    return (ppc.summary)
  }

  getPpcM1Q <- function(gt.data, p, hdi.level, s) {
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
      errors <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)

      for(i in 1:gt.data$N) {
        beta.p <- paste("beta", t, s, gt.data$K[i], sep = '.')
        # TODO: check in case single Ntq => sigma no index
        if(gt.data$Nk == 1) {
          beta.p <- paste("beta", t, s, sep = '.')
        }

        yhat[, i] <- apply(X = p[, c(alpha.p, beta.p, sigma.p)], MARGIN = 1,
                           FUN = getMu, y = gt.data$X[i, s])
        errors[, i] <- abs(gt.data$Yq[i, t] - yhat[, i])
      }

      # strain level
      for(k in 1:gt.data$Nk) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$K == k]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$K == k]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = k, x = NA, level = "strain",
                          y.real = mean(gt.data$Yq[gt.data$K == k, t]),
                          y.ppc = mean(yhat[, gt.data$K == k]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$K == k]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }

      # snp level
      if(sum(gt.data$X[, s] == 1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == 1]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$X[, s] == 1]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = NA, x = 1, level = "snp",
                          y.real = mean(gt.data$Yq[gt.data$X[, s] == 1, t]),
                          y.ppc = mean(yhat[, gt.data$X[, s] == 1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$X[, s] == 1]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
      # snp level
      if(sum(gt.data$X[, s] == -1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == -1]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$X[, s] == -1]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = NA, x = -1, level = "snp",
                          y.real = mean(gt.data$Yq[gt.data$X[, s] == -1, t]),
                          y.ppc = mean(yhat[, gt.data$X[, s] == -1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$X[, s] == -1]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
    }
    # clean up
    rm(k, yhat, yhat.hdi, row)

    return (ppc.summary)
  }

  getPpcM2Q <- function(gt.data, p, hdi.level, s) {
    ppc.summary <- c()

    # Q-trait
    for(t in 1:gt.data$Ntq) {
      alpha.p <- paste("alpha", t, sep = '.')

      # individual level
      yhat <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)
      errors <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)

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
        errors[, i] <- abs(gt.data$Yq[i, t] - yhat[, i])
      }

      # strain level
      for(k in 1:gt.data$Nk) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$K == k]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$K == k]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = k, x = NA, level = "strain",
                          y.real = mean(gt.data$Yq[gt.data$K == k, t]),
                          y.ppc = mean(yhat[, gt.data$K == k]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$K == k]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }

      # snp level
      if(sum(gt.data$X[, s] == 1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == 1]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$X[, s] == 1]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = NA, x = 1, level = "snp",
                          y.real = mean(gt.data$Yq[gt.data$X[, s] == 1, t]),
                          y.ppc = mean(yhat[, gt.data$X[, s] == 1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$X[, s] == 1]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
      # snp level
      if(sum(gt.data$X[, s] == -1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == -1]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$X[, s] == -1]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = NA, x = -1, level = "snp",
                          y.real = mean(gt.data$Yq[gt.data$X[, s] == -1, t]),
                          y.ppc = mean(yhat[, gt.data$X[, s] == -1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$X[, s] == -1]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
    }
    # clean up
    rm(k, yhat, yhat.hdi, row)

    return (ppc.summary)
  }

  if(model == "M0" | model == "M0c") {
    return(getPpcM0Q(gt.data = gt.data, p = ext,
                     hdi.level = hdi.level, s = s))
  }
  else if(model == "M1" | model == "M1c") {
    return(getPpcM1Q(gt.data = gt.data, p = ext,
                     hdi.level = hdi.level, s = s))
  }
  else if(model == "M2" | model == "M2c") {
    return(getPpcM2Q(gt.data = gt.data, p = ext,
                     hdi.level = hdi.level, s = s))
  }
}



getPpcD <- function(ext, gt.data, model, s, hdi.level) {

  getPr <- function(x, y) {
    return(1/(1 + exp(-(x[1]+x[2]*y))))
  }

  getPpcM0D <- function(gt.data, p, hdi.level, s) {
    ppc.summary <- c()

    # D-trait
    for(d in 1:gt.data$Ntd) {
      t <- d + gt.data$Ntq

      alpha.p <- paste("alpha", t, sep = '.')

      # individual level
      yhat <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)
      errors <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)

      for(i in 1:gt.data$N) {
        beta.p <- paste("beta", t, s, sep = '.')

        yhat[, i] <- apply(X = p[, c(alpha.p, beta.p)], MARGIN = 1,
                           FUN = getPr, y = gt.data$X[i, s])
        errors[, i] <- abs(gt.data$Yd[i, d] - yhat[, i])
      }

      # snp level
      if(sum(gt.data$X[, s] == 1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == 1]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$X[, s] == 1]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = NA, x = 1, level = "snp",
                          y.real=mean(as.numeric(gt.data$Yd[gt.data$X[, s] == 1,d])),
                          y.ppc = mean(yhat[, gt.data$X[, s] == 1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$X[, s] == 1]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
      # snp level
      if(sum(gt.data$X[, s] == -1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == -1]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$X[, s] == -1]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = NA, x = -1, level = "snp",
                          y.real=mean(as.numeric(gt.data$Yd[gt.data$X[, s] == -1,d])),
                          y.ppc = mean(yhat[, gt.data$X[, s] == -1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$X[, s] == -1]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
    }
    # clean up
    rm(i, yhat, yhat.hdi, row)

    return (ppc.summary)
  }

  getPpcM1D <- function(gt.data, p, hdi.level, s) {
    ppc.summary <- c()

    # D-trait
    for(d in 1:gt.data$Ntd) {
      t <- d + gt.data$Ntq

      alpha.p <- paste("alpha", t, sep = '.')

      # individual level
      yhat <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)
      errors <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)

      for(i in 1:gt.data$N) {
        beta.p <- paste("beta", t, s, gt.data$K[i], sep = '.')
        # TODO: check in case single Ntq => sigma no index
        if(gt.data$Nk == 1) {
          beta.p <- paste("beta", t, s, sep = '.')
        }

        yhat[, i] <- apply(X = p[, c(alpha.p, beta.p)], MARGIN = 1,
                           FUN = getPr, y = gt.data$X[i, s])
        errors[, i] <- abs(gt.data$Yd[i, d] - yhat[, i])
      }

      # strain level
      for(k in 1:gt.data$Nk) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$K == k]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$K == k]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = k, x = NA, level = "strain",
                          y.real=mean(as.numeric(gt.data$Yd[gt.data$K == k,d])),
                          y.ppc = mean(yhat[, gt.data$K == k]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$K == k]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }

      # snp level
      if(sum(gt.data$X[, s] == 1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == 1]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$X[, s] == 1]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = NA, x = 1, level = "snp",
                          y.real=mean(as.numeric(gt.data$Yd[gt.data$X[, s]==1,d])),
                          y.ppc = mean(yhat[, gt.data$X[, s] == 1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$X[, s] == 1]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
      # snp level
      if(sum(gt.data$X[, s] == -1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == -1]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$X[, s] == -1]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = NA, x = -1, level = "snp",
                          y.real=mean(as.numeric(gt.data$Yd[gt.data$X[, s]==-1,d])),
                          y.ppc = mean(yhat[, gt.data$X[, s] == -1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$X[, s] == -1]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
    }
    # clean up
    rm(k, yhat, yhat.hdi, row)

    return (ppc.summary)
  }

  getPpcM2D <- function(gt.data, p, hdi.level, s) {
    ppc.summary <- c()

    # D-trait
    for(d in 1:gt.data$Ntd) {
      t <- d + gt.data$Ntq

      alpha.p <- paste("alpha", t, sep = '.')

      # individual level
      yhat <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)
      errors <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)

      for(i in 1:gt.data$N) {
        beta.p <- paste("beta", t, s, gt.data$K[i], sep = '.')
        # TODO: check in case single Ntq => sigma no index
        if(gt.data$Nk == 1) {
          beta.p <- paste("beta", t, s, sep = '.')
        }

        yhat[, i] <- apply(X = p[, c(alpha.p, beta.p)], MARGIN = 1,
                           FUN = getPr, y = gt.data$X[i, s])
        errors[, i] <- abs(gt.data$Yd[i, d] - yhat[, i])
      }

      # strain level
      for(k in 1:gt.data$Nk) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$K == k]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$K == k]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = k, x = NA, level = "strain",
                          y.real=mean(as.numeric(gt.data$Yd[gt.data$K == k,d])),
                          y.ppc = mean(yhat[, gt.data$K == k]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$K == k]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }

      # snp level
      if(sum(gt.data$X[, s] == 1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == 1]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$X[, s] == 1]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = NA, x = 1, level = "snp",
                          y.real=mean(as.numeric(gt.data$Yd[gt.data$X[, s] == 1,d])),
                          y.ppc = mean(yhat[, gt.data$X[, s] == 1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$X[, s] == 1]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
      # snp level
      if(sum(gt.data$X[, s] == -1) != 0) {
        yhat.hdi <- getHdi(vec = as.vector(yhat[, gt.data$X[, s] == -1]),
                           hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = as.vector(errors[, gt.data$X[, s] == -1]),
                             hdi.level = hdi.level)

        row <- data.frame(t = t, s = s, k = NA, x = -1, level = "snp",
                          y.real=mean(as.numeric(gt.data$Yd[gt.data$X[, s] == -1,d])),
                          y.ppc = mean(yhat[, gt.data$X[, s] == -1]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error = mean(errors[, gt.data$X[, s] == -1]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)
      }
    }
    # clean up
    rm(k, yhat, yhat.hdi, row)

    return (ppc.summary)
  }

  if(model == "M0" | model == "M0c") {
    return(getPpcM0D(gt.data = gt.data, p = ext,
                     hdi.level = hdi.level, s = s))
  }
  else if(model == "M1" | model == "M1c") {
    return(getPpcM1D(gt.data = gt.data, p = ext,
                     hdi.level = hdi.level, s = s))
  }
  else if(model == "M2" | model == "M2c") {
    return(getPpcM2D(gt.data = gt.data, p = ext,
                     hdi.level = hdi.level, s = s))
  }
}

