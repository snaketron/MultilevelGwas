

# Function:
# Use the mid-level parameters (SNP effects) to make predictions
getPpcMid <- function(ext, gt.data, s,
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

    stats <- c()
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

          # hdi's
          yhat.hdi <- getHdi(vec = yhat, hdi.level = hdi.level)
          error <- abs(yhat-y.real)
          error.hdi <- getHdi(vec = error, hdi.level = hdi.level)

          row <- data.frame(t = t,
                            s = s,
                            x = x,
                            k = NA,
                            i = NA,
                            parameter.level = "mid",
                            prediction.level = "snp",
                            y.real.mean = y.real,
                            y.ppc.mean = mean(x = yhat),
                            y.ppc.median = median(x = yhat),
                            y.ppc.L = yhat.hdi[1],
                            y.ppc.H = yhat.hdi[2],
                            y.error.mean = mean(error),
                            y.error.median = median(error),
                            y.error.L = error.hdi[1],
                            y.error.H = error.hdi[2])

          # cleanup
          rm(y.real, yhat.hdi, xs, ps, yhat)
        }
        else {
          row <- data.frame(t = t,
                            s = s,
                            k = NA,
                            i = NA,
                            x = x,
                            parameter.level = "mid",
                            prediction.level = "snp",
                            y.real.mean = NA,
                            y.ppc.mean = NA,
                            y.ppc.median = NA,
                            y.ppc.L = NA,
                            y.ppc.H = NA,
                            y.error.mean = NA,
                            y.error.median = NA,
                            y.error.L = NA,
                            y.error.H = NA)
        }
        stats <- rbind(stats, row)

        # cleanup
        rm(row)
      }
    }

    stats$model <- model

    return (stats)
  }


  return(getPpc(p = ext, gt.data = gt.data,
                s = s, hdi.level = hdi.level,
                model = model))
}



# Function:
# Posterior prediction
getPpcLow <- function(ext, gt.data,
                      model, hdi.level) {

  # Function:
  # Use the most specific parameters to make predictions at different levels
  getPpcPosterior <- function(ext, gt.data, s, model) {

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


    getMuSnpIndividual <- function(x, y, trait.type) {
      if(trait.type == "D") {
        return(rbinom(n = 1, size = 1, prob = 1/(1 + exp(-(x[1]+x[2]*y)))))
      }
      if(trait.type == "Q") {
        return(rnorm(n = 1, mean = x[1]+x[2]*y, sd = x[3]))
      }
    }


    getPpc <- function(p, gt.data, s, model) {

      out <- vector(mode = "list", length = gt.data$Ntq+gt.data$Ntd)

      for(t in 1:(gt.data$Ntq+gt.data$Ntd)) {
        # matrices of results
        yhat <- matrix(data = NA, nrow = nrow(p), ncol = gt.data$N)

        tt <- gt.data$trait.type[t]


        # Prediction at I
        for(i in 1:gt.data$N) {

          ps <- getPar(t = t, s = s, k = gt.data$K[i],
                       gt.data = gt.data, model = model)

          yhat[, i] <- apply(X = p[, ps], MARGIN = 1,
                             FUN = getMuSnpIndividual,
                             y = gt.data$X[i, s],
                             trait.type = tt)
        }

        out[[t]] <- yhat
      }

      return (out)
    }


    return(getPpc(p = ext, gt.data = gt.data,
                  s = s, model = model))
  }



  ppc.out <- vector(mode = "list", length = gt.data$Ns)
  cat("Model:", model, "\n", sep = '')
  for(s in 1:gt.data$Ns) {
    ppc.out[[s]] <- getPpcPosterior(ext = ext,
                                    gt.data = gt.data,
                                    model = model,
                                    s = s)

    if(s %% 50 == 0) {
      cat(s, "/", gt.data$Ns, ',',  sep = '')
    }
  }
  cat("\n")



  ys <- array(data = NA, dim = c(gt.data$Ntd+gt.data$Ntq,
                                 gt.data$Ns,
                                 gt.data$N,
                                 nrow(ext)))

  yi <- array(data = NA, dim = c(gt.data$Ntd+gt.data$Ntq,
                                 gt.data$N,
                                 nrow(ext)))

  es <- array(data = NA, dim = c(gt.data$Ntd+gt.data$Ntq,
                                 gt.data$N,
                                 nrow(ext)))

  stats <- c()
  for(t in 1:(gt.data$Ntd+gt.data$Ntq)) {
    for(i in 1:gt.data$N) {
      for(s in 1:gt.data$Ns) {
        ys[t,s,i,] <- ppc.out[[s]][[t]][, i]
      }

      # stats at individual level
      yi[t,i,] <- apply(X = t(ys[t,,i,]), MARGIN = 1, FUN = mean)
      yhat.hdi <- getHdi(vec = yi[t,i,], hdi.level = hdi.level)
      es[t,i,] <- abs(yi[t,i,]-gt.data$Y[i,t])
      es.hdi <- getHdi(vec = es[t,i,], hdi.level = hdi.level)

      row <- data.frame(t = t,
                        s = NA,
                        x = NA,
                        k = NA,
                        i = i,
                        parameter.level = "low",
                        prediction.level = "individual",
                        y.real.mean = gt.data$Y[i, t],
                        y.ppc.mean = mean(yi[t,i,]),
                        y.ppc.median = median(yi[t,i,]),
                        y.ppc.L = yhat.hdi[1],
                        y.ppc.H = yhat.hdi[2],
                        y.error.mean = mean(es[t,i,]),
                        y.error.median = median(es[t,i,]),
                        y.error.L = es.hdi[1],
                        y.error.H = es.hdi[2])
      stats <- rbind(stats, row)
    }

    # strains level
    for(k in 1:gt.data$Nk) {
      ks <- which(gt.data$K == k)

      # stats at strain level
      yhat <- apply(X = t(yi[t,ks,]), MARGIN = 1, FUN = mean)
      yhat.hdi <- getHdi(vec = yhat, hdi.level = hdi.level)
      e <- apply(X = t(es[t,ks,]), MARGIN = 1, FUN = mean)
      es.hdi <- getHdi(vec = e, hdi.level = hdi.level)

      row <- data.frame(t = t,
                        s = NA,
                        x = NA,
                        k = k,
                        i = i,
                        parameter.level = "low",
                        prediction.level = "strain",
                        y.real.mean = mean(gt.data$Y[ks, t]),
                        y.ppc.mean = mean(yhat),
                        y.ppc.median = median(yhat),
                        y.ppc.L = yhat.hdi[1],
                        y.ppc.H = yhat.hdi[2],
                        y.error.mean = mean(e),
                        y.error.median = median(e),
                        y.error.L = es.hdi[1],
                        y.error.H = es.hdi[2])
      stats <- rbind(stats, row)
    }



    # stats at snp level
    for(s in 1:gt.data$Ns) {
      for(x in c(1, -1)) {
        xs <- which(gt.data$X[, s] == x)
        if(length(xs) != 0) {

          # stats at strain level
          yhat <- apply(X = t(yi[t,xs,]), MARGIN = 1, FUN = mean)
          yhat.hdi <- getHdi(vec = yhat, hdi.level = hdi.level)
          e <- apply(X = t(es[t,xs,]), MARGIN = 1, FUN = mean)
          es.hdi <- getHdi(vec = e, hdi.level = hdi.level)

          row <- data.frame(t = t,
                            s = s,
                            x = x,
                            k = NA,
                            i = NA,
                            parameter.level = "low",
                            prediction.level = "snp",
                            y.real.mean = mean(gt.data$Y[xs, t]),
                            y.ppc.mean = mean(yhat),
                            y.ppc.median = median(yhat),
                            y.ppc.L = yhat.hdi[1],
                            y.ppc.H = yhat.hdi[2],
                            y.error.mean = mean(e),
                            y.error.median = median(e),
                            y.error.L = es.hdi[1],
                            y.error.H = es.hdi[2])
          stats <- rbind(stats, row)
        }
      }
    }
  }

  stats$model <- model
  return (stats)
  # return (list(ys = ys, es = es, stats = stats))
}


