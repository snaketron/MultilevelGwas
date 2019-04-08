


# Function:
# Use the most specific parameters to make predictions at different levels
getPpcLowestLevel <- function(ext, gt.data, s,
                              hdi.level, model) {

  getPar <- function(t, s, k, gt.data, model) {
    if(model == "M0" | model == "M0c") {
      # alpha
      alpha.p <- paste("alpha", t, sep = '.')
      if(gt.data$Ntq+gt.data$Ntd == 1) {
        alpha.p <- "alpha"
      }

      # beta
      beta.p <- paste("beta", t, s, sep = '.')

      # Q
      if(gt.data$trait.type[t] == "Q") {
        sigma.p <- paste("sigma", t, sep = '.')
        if(gt.data$Ntq == 1) {
          sigma.p <- "sigma"
        }
        par <- c(alpha.p, beta.p, sigma.p)
        return (par)
      }
      # D
      if(gt.data$trait.type[t] == "D") {
        par <- c(alpha.p, beta.p)
        return (par)
      }
    }

    if(model == "M1" | model == "M1c") {

      # alpha
      alpha.p <- paste("alpha", t, sep = '.')
      if(gt.data$Ntq+gt.data$Ntd == 1) {
        alpha.p <- "alpha"
      }

      # beta
      beta.p <- paste("beta", t, s, k, sep = '.')
      if(gt.data$Nk == 1) {
        beta.p <- paste("beta", t, s, sep = '.')
        if(gt.data$Ntq+gt.data$Ntd == 1) {
          beta.p <- paste("beta", s, sep = '.')
        }
      }

      # Q
      if(gt.data$trait.type[t] == "Q") {
        # sigma
        sigma.p <- paste("sigma", t, sep = '.')
        if(gt.data$Ntq == 1) {
          sigma.p <- "sigma"
        }
        par <- c(alpha.p, beta.p, sigma.p)
        return (par)
      }
      # D
      if(gt.data$trait.type[t] == "D") {
        par <- c(alpha.p, beta.p)
        return (par)
      }
    }

    if(model == "M2" | model == "M2c") {

      # alpha
      alpha.p <- paste("alpha", t, sep = '.')
      if(gt.data$Ntq+gt.data$Ntd == 1) {
        alpha.p <- "alpha"
      }

      # beta
      beta.p <- paste("beta", t, s, k, sep = '.')
      if(gt.data$Nk == 1) {
        beta.p <- paste("beta", t, s, sep = '.')
        if(gt.data$Ntq+gt.data$Ntd == 1) {
          beta.p <- paste("beta", s, sep = '.')
        }
      }

      # Q
      if(gt.data$trait.type[t] == "Q") {
        # sigma
        sigma.p <- paste("sigma", k, sep = '.')
        if(gt.data$Nk == 1) {
          sigma.p <- "sigma"
        }
        par <- c(alpha.p, beta.p, sigma.p)
        return (par)
      }
      # D
      if(gt.data$trait.type[t] == "D") {
        par <- c(alpha.p, beta.p)
        return (par)
      }
    }
  }


  getMuSnp <- function(x, y, trait.type) {
    if(trait.type == "D") {
      return(1/(1 + exp(-(x[1]+x[2]*y))))
      # return(rbinom(n = 1, size = 1, prob = 1/(1 + exp(-(x[1]+x[2]*y)))))
    }
    if(trait.type == "Q") {
      return(x[1]+x[2]*y)
      # return(rnorm(n = 1, mean = x[1]+x[2]*y, sd = x[3]))
    }
  }


  getPpc <- function(p, gt.data, s,
                     hdi.level, model) {

    # Make predictions for each individual, based on
    # the lowest level parameters in M1 and M2
    ppc.summary <- c()

    for(t in 1:(gt.data$Ntq+gt.data$Ntd)) {
      # index of D-trait
      d <- sum(gt.data$trait.type[1:t] == "D")

      # matrices of results
      yhat <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)
      errors <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)

      tt <- gt.data$trait.type[t]


      # Prediction at I
      for(i in 1:gt.data$N) {

        ps <- getPar(t = t, s = s, k = gt.data$K[i],
                     gt.data = gt.data, model = model)

        yhat[, i] <- apply(X = p[, ps], MARGIN = 1, FUN = getMuSnp,
                           y = gt.data$X[i, s], trait.type = tt)

        if(tt == "Q") {
          y.real <- gt.data$Yq[i, t]
        }
        if(tt == "D") {
          y.real <- as.numeric(gt.data$Yd[i, d])
        }

        # compute errors posterior
        errors[, i] <- abs(y.real - yhat[, i])

        # hdi's
        yhat.hdi <- getHdi(vec = yhat[, i], hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = errors[, i], hdi.level = hdi.level)

        row <- data.frame(t = t,
                          s = s,
                          k = gt.data$K[i],
                          i = i,
                          x = gt.data$X[i, s],
                          parameter.level = "lowest",
                          prediction.level = "individual-level",
                          y.real.mean = y.real,
                          y.ppc.mean = mean(yhat[, i]),
                          y.ppc.median = median(yhat[, i]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error.mean = mean(errors[, i]),
                          error.median = median(errors[, i]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)

        # cleanup
        rm(row, yhat.hdi, errors.hdi, y.real, ps)
      }



      # Prediction at: K
      for(k in 1:gt.data$Nk) {
        ks <- which(gt.data$K == k)

        if(tt == "Q") {
          y.real <- gt.data$Yq[gt.data$K == k, t]
        }
        if(tt == "D") {
          y.real <- as.numeric(gt.data$Yd[gt.data$K == k, d])
        }

        # hdi's
        yhat.hdi <- getHdi(vec = yhat[, ks], hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = errors[, ks], hdi.level = hdi.level)

        row <- data.frame(t = t,
                          s = s,
                          k = k,
                          i = NA,
                          x = NA,
                          parameter.level = "lowest",
                          prediction.level = "strain-level",
                          y.real.mean = mean(y.real),
                          y.ppc.mean = mean(yhat[, ks]),
                          y.ppc.median = median(yhat[, ks]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error.mean = mean(errors[, ks]),
                          error.median = median(errors[, ks]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)

        # cleanup
        rm(row, yhat.hdi, errors.hdi, k, ks)
      }


      # Prediction at: S
      for(x in c(1, -1)) {
        xs <- which(gt.data$X[, s] == x)
        if(tt == "Q") {
          y.real <- gt.data$Yq[xs, t]
        }
        if(tt == "D") {
          y.real <- as.numeric(gt.data$Yd[xs, d])
        }

        # hdi's
        yhat.hdi <- getHdi(vec = yhat[, xs], hdi.level = hdi.level)
        errors.hdi <- getHdi(vec = errors[, xs], hdi.level = hdi.level)

        row <- data.frame(t = t,
                          s = s,
                          k = NA,
                          i = NA,
                          x = x,
                          parameter.level = "lowest",
                          prediction.level = "snp-level",
                          y.real.mean = mean(y.real),
                          y.ppc.mean = mean(yhat[, xs]),
                          y.ppc.median = median(yhat[, xs]),
                          y.ppc.L = yhat.hdi[1],
                          y.ppc.H = yhat.hdi[2],
                          error.mean = mean(errors[, xs]),
                          error.median = median(errors[, xs]),
                          error.L = errors.hdi[1],
                          error.H = errors.hdi[2])
        ppc.summary <- rbind(ppc.summary, row)

        # cleanup
        rm(row, yhat.hdi, errors.hdi, xs)
      }
    }

    return (ppc.summary)
  }


  return(getPpc(p = ext, gt.data = gt.data,
                s = s, hdi.level = hdi.level,
                model = model))
}



# Function:
# Use the mid-level parameters (SNP effects) to make predictions
getPpcMidLevel <- function(ext, gt.data, s,
                           hdi.level, model) {

  getPar <- function(t, s, k, gt.data, model) {
    # alpha
    alpha.p <- paste("alpha", t, sep = '.')
    if(gt.data$Ntq+gt.data$Ntd == 1) {
      alpha.p <- "alpha"
    }

    # beta
    beta.p <- paste("mu_beta", t, s, sep = '.')
    if(gt.data$Ntq+gt.data$Ntd == 1) {
      beta.p <- paste("mu_beta", t, s, sep = '.')
    }

    # sigma
    sigma.p <- paste("sigma_beta", t, sep = '.')
    if(gt.data$Ntq == 1) {
      sigma.p <- "sigma_beta"
    }


    par <- c(alpha.p, beta.p, sigma.p)
    return (par)
  }


  getMuSnp <- function(x, y, trait.type) {
    # m <- rnorm(n = 1, mean = x[1]+x[2]*y, sd = x[3])
    m <- x[1]+x[2]*y
    if(trait.type == "D") {
      m <- 1/(1 + exp(-m))
    }
    return(m)
  }


  getPpc <- function(p, gt.data, s,
                     hdi.level, model) {

    # Make predictions for each individual, based on the lowest level + 1
    # parameters in M1 and M2
    ppc.summary <- c()

    # Prediction at: X
    for(t in 1:(gt.data$Ntq+gt.data$Ntd)) {
      # index of D-trait
      d <- sum(gt.data$trait.type[1:t] == "D")

      tt <- gt.data$trait.type[t]
      for(x in c(1, -1)) {
        ps <- getPar(t = t, s = s, k = gt.data$K[i],
                     gt.data = gt.data, model = model)

        yhat <- apply(X = p[, ps], MARGIN = 1, FUN = getMuSnp,
                      y = x, trait.type = tt)

        # real data
        xs <- which(gt.data$X[, s] == x)
        if(length(xs) != 0) {
          if(tt == "Q") {
            y.real <- mean(gt.data$Yq[xs, t])
          }
          if(tt == "D") {
            y.real <- mean(as.numeric(gt.data$Yd[xs, d]))
          }

          # compute errors posterior
          errors <- abs(y.real - yhat)

          # hdi's
          yhat.hdi <- getHdi(vec = yhat, hdi.level = hdi.level)
          errors.hdi <- getHdi(vec = errors, hdi.level = hdi.level)

          row <- data.frame(t = t,
                            s = s,
                            k = NA,
                            i = NA,
                            x = x,
                            parameter.level = "lowest-plus-one",
                            prediction.level = "snp-level",
                            y.real.mean = y.real,
                            y.ppc.mean = mean(x = yhat),
                            y.ppc.median = median(x = yhat),
                            y.ppc.L = yhat.hdi[1],
                            y.ppc.H = yhat.hdi[2],
                            error.mean = mean(x = errors),
                            error.median = median(x = errors),
                            error.L = errors.hdi[1],
                            error.H = errors.hdi[2])



          # cleanup
          rm(errors.hdi)
        }
        else {

          # hdi's
          yhat.hdi <- getHdi(vec = yhat, hdi.level = hdi.level)

          row <- data.frame(t = t,
                            s = s,
                            k = NA,
                            i = NA,
                            x = x,
                            parameter.level = "lowest-plus-one",
                            prediction.level = "snp-level",
                            y.real.mean = NA,
                            y.ppc.mean = NA,
                            y.ppc.median = median(x = yhat),
                            y.ppc.L = yhat.hdi[1],
                            y.ppc.H = yhat.hdi[2],
                            error.mean = NA,
                            error.median = NA,
                            error.L = NA,
                            error.H = NA)
        }
        ppc.summary <- rbind(ppc.summary, row)

        # cleanup
        rm(row, y.real, ps, yhat.hdi)
      }
    }

    return (ppc.summary)
  }


  return(getPpc(p = ext, gt.data = gt.data,
                s = s, hdi.level = hdi.level,
                model = model))
}



