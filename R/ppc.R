

# Function:
# Posterior prediction horizontal (at the same level, except lowest)
getPpcHorizontal <- function(ext, gt.data,
                             model, hdi.level) {


  getPar <- function(t, s, k, gt.data,
                     model, level,
                     trait.type) {
    # alpha
    alpha.p <- paste("alpha", t, sep = '.')
    if(gt.data$Ntq+gt.data$Ntd == 1) {
      alpha.p <- "alpha"
    }

    # only across SNP level
    if(model %in% c("M0", "M0c")) {
      # beta
      beta.p <- paste("beta", t, s, sep = '.')

      # sigma
      sigma.p <- paste("sigma", t, sep = '.')
      if(gt.data$Ntq == 1) {
        sigma.p <- "sigma"
      }

      if(trait.type == "Q") {
        par <- c(alpha.p, beta.p, sigma.p)
      }
      else {
        par <- c(alpha.p, beta.p)
      }
    }

    if(model %in% c("M1", "M1c")) {

      # strain
      if(level == "strain") {
        # beta
        beta.p <- paste("beta", t, s, k, sep = '.')

        # sigma
        sigma.p <- paste("sigma", t, sep = '.')
        if(gt.data$Ntq == 1) {
          sigma.p <- "sigma"
        }

        if(trait.type == "Q") {
          par <- c(alpha.p, beta.p, sigma.p)
        } else {
          par <- c(alpha.p, beta.p)
        }
      }

      # snp
      if(level == "snp") {
        # beta
        beta.p <- paste("mu_beta", t, s, sep = '.')

        # sigma
        sigma.p <- paste("sigma_beta", t, sep = '.')
        if(gt.data$Ntq == 1) {
          sigma.p <- "sigma_beta"
        }

        par <- c(alpha.p, beta.p, sigma.p)
      }
    }

    return (par)
  }


  getMuSnp <- function(x, y,
                       trait.type,
                       level) {

    if(level == "strain") {
      if(trait.type == "Q") {
        m <- rnorm(n = 1, mean = x[1]+x[2]*y, sd = x[3])
      }
      if(trait.type == "D") {
        m <- rbinom(n = 1, size = 1, prob = 1/(1 + exp(-(x[1]+x[2]*y))))
      }
    }

    if(level == "snp") {
      if(model %in% c("M0", "M0c")) {
        if(trait.type == "Q") {
          m <- rnorm(n = 1, mean = x[1]+x[2]*y, sd = x[3])
        }
        if(trait.type == "D") {
          m <- rbinom(n = 1, size = 1, prob = 1/(1 + exp(-(x[1]+x[2]*y))))
        }
      }
      if(model %in% c("M1", "M1c")) {
        m <- rnorm(n = 1, mean = x[1]+x[2]*y, sd = x[3])
        if(trait.type == "D") {
          m <- rbinom(n = 1, size = 1, prob = 1/(1 + exp(-m)))
        }
      }
    }


    return(m)
  }


  # compute statistics
  ysk <- array(data = NA, dim = c(gt.data$Ntd+gt.data$Ntq,
                                 gt.data$Ns,
                                 gt.data$Nk,
                                 nrow(ext)))

  ys <- array(data = NA, dim = c(gt.data$Ntd+gt.data$Ntq,
                                  gt.data$Ns,
                                  2,
                                  nrow(ext)))


  stats <- c()
  for(t in 1:(gt.data$Ntq+gt.data$Ntd)) {

    for(s in 1:gt.data$Ns) {

      # get param from posterior
      ps <- getPar(t = t, s = s, k = NA,
                   gt.data = gt.data,
                   model = model,
                   level = "snp",
                   trait.type = gt.data$trait.type[t])

      # snp level
      ys[t,s,1,] <- apply(X = ext[, ps], MARGIN = 1,
                          FUN = getMuSnp, y = 1,
                          trait.type = gt.data$trait.type[t],
                          level = "snp")
      ys[t,s,2,] <- apply(X = ext[, ps], MARGIN = 1,
                          FUN = getMuSnp, y = -1,
                          trait.type = gt.data$trait.type[t],
                          level = "snp")


      # special case for strains
      if(model %in% c("M1", "M1c")) {
        for(k in 1:gt.data$Nk) {

          # get param from posterior
          ps <- getPar(t = t, s = s, k = k,
                       gt.data = gt.data,
                       model = model,
                       level = "strain",
                       trait.type = gt.data$trait.type[t])

          ysk[t,s,k,] <- apply(X = ext[, ps], MARGIN = 1,
                               FUN = getMuSnp,
                               y = unique(gt.data$X[gt.data$K == k, s]),
                               trait.type = gt.data$trait.type[t],
                               level = "strain")

          y.real <- mean(gt.data$Y[gt.data$K == k, t])
          yhat.hdi <- getHdi(vec = ysk[t,s,k,], hdi.level = hdi.level)
          es <- abs(ysk[t,s,k,]-y.real)
          es.hdi <- getHdi(vec = es, hdi.level = hdi.level)
          row <- data.frame(t = t,
                            s = s,
                            x = unique(gt.data$X[gt.data$K == k, s]),
                            k = k,
                            i = NA,
                            parameter.level = "mid",
                            prediction.level = "strain",
                            y.real.mean = y.real,
                            y.ppc.mean = mean(ysk[t,s,k,]),
                            y.ppc.median = median(ysk[t,s,k,]),
                            y.ppc.L = yhat.hdi[1],
                            y.ppc.H = yhat.hdi[2],
                            y.error.mean = mean(es),
                            y.error.median = median(es),
                            y.error.L = es.hdi[1],
                            y.error.H = es.hdi[2],
                            model = model)
          stats <- rbind(stats, row)
        }
      }


      # X = 1
      y.real <- mean(gt.data$Y[gt.data$X[, s] == 1, t])
      yhat.hdi <- getHdi(vec = ys[t,s,1,], hdi.level = hdi.level)
      es <- abs(ys[t,s,1,]-y.real)
      es.hdi <- getHdi(vec = es, hdi.level = hdi.level)
      row <- data.frame(t = t,
                        s = s,
                        x = 1,
                        k = NA,
                        i = NA,
                        parameter.level = "mid",
                        prediction.level = "snp",
                        y.real.mean = y.real,
                        y.ppc.mean = mean(ys[t,s,1,]),
                        y.ppc.median = median(ys[t,s,1,]),
                        y.ppc.L = yhat.hdi[1],
                        y.ppc.H = yhat.hdi[2],
                        y.error.mean = mean(es),
                        y.error.median = median(es),
                        y.error.L = es.hdi[1],
                        y.error.H = es.hdi[2],
                        model = model)
      stats <- rbind(stats, row)



      # X = -1
      y.real <- mean(gt.data$Y[gt.data$X[, s] == -1, t])
      yhat.hdi <- getHdi(vec = ys[t,s,2,], hdi.level = hdi.level)
      es <- abs(ys[t,s,2,]-y.real)
      es.hdi <- getHdi(vec = es, hdi.level = hdi.level)
      row <- data.frame(t = t,
                        s = s,
                        x = -1,
                        k = NA,
                        i = NA,
                        parameter.level = "mid",
                        prediction.level = "snp",
                        y.real.mean = y.real,
                        y.ppc.mean = mean(ys[t,s,2,]),
                        y.ppc.median = median(ys[t,s,2,]),
                        y.ppc.L = yhat.hdi[1],
                        y.ppc.H = yhat.hdi[2],
                        y.error.mean = mean(es),
                        y.error.median = median(es),
                        y.error.L = es.hdi[1],
                        y.error.H = es.hdi[2],
                        model = model)
      stats <- rbind(stats, row)


      if(s %% 50 == 0) {
        cat("trait:", t, "/", gt.data$Ntq+gt.data$Ntd,
            ", SNP:", s, "/", gt.data$Ns, '\n', sep = '')
      }
    }
  }

  return (stats)
}



# Function:
# Posterior prediction from bottom to top
getPpcHierarchical <- function(ext, gt.data,
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
  for(s in 1:gt.data$Ns) {
    ppc.out[[s]] <- getPpcPosterior(ext = ext,
                                    gt.data = gt.data,
                                    model = model,
                                    s = s)

    if(s %% 50 == 0) {
      cat("SNP:", s, "/", gt.data$Ns, ',', sep = '')
    }
  }




  # compute statistics
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
                        y.error.H = es.hdi[2],
                        model = model)
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
                        y.error.H = es.hdi[2],
                        model = model)
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
                            y.error.H = es.hdi[2],
                            model = model)
          stats <- rbind(stats, row)
        }
      }
    }
  }
  return (stats)
}



