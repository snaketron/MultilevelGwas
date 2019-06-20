

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
      if(level == "refSNP") {
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


  getMuSnp <- function(x, y, trait.type, level) {

    if(level == "strain") {
      if(trait.type == "Q") {
        m <- rnorm(n = 1, mean = x[1]+x[2]*y, sd = x[3])
      }
      if(trait.type == "D") {
        # m <- 1/(1 + exp(-(x[1]+x[2]*y)))
        m <- rbinom(n = 1, size = 1, prob = 1/(1 + exp(-(x[1]+x[2]*y))))
      }
    }

    if(level == "refSNP") {
      if(model %in% c("M0", "M0c")) {
        if(trait.type == "Q") {
          m <- rnorm(n = 1, mean = x[1]+x[2]*y, sd = x[3])
        }
        if(trait.type == "D") {
          # m <- 1/(1 + exp(-(x[1]+x[2]*y)))
          m <- mean(rbinom(n = 1, size = 1, prob = 1/(1 + exp(-(x[1]+x[2]*y)))))
        }
      }
      if(model %in% c("M1", "M1c")) {
        m <- rnorm(n = 1, mean = x[1]+x[2]*y, sd = x[3])
        if(trait.type == "D") {
          # m <- 1/(1 + exp(-m))
          m <- rbinom(n = 1, size = 1, prob = 1/(1 + exp(-m)))
        }
      }
    }

    return(m)
  }



  getMuSnpOther <- function(x, y, trait.type, level) {

    if(level == "strain") {
      if(trait.type == "Q") {
        m <- rnorm(n = 1, mean = x[1]+x[2]*y, sd = x[3])
      }
      if(trait.type == "D") {
        m <- 1/(1 + exp(-(x[1]+x[2]*y)))
        # m <- rbinom(n = 1, size = 1, prob = 1/(1 + exp(-(x[1]+x[2]*y))))
      }
    }

    if(level == "refSNP") {
      if(model %in% c("M0", "M0c")) {
        if(trait.type == "Q") {
          m <- rnorm(n = 1, mean = x[1]+x[2]*y, sd = x[3])
        }
        if(trait.type == "D") {
          m <- 1/(1 + exp(-(x[1]+x[2]*y)))
          # m <- mean(rbinom(n = 1, size = 1, prob = 1/(1 + exp(-(x[1]+x[2]*y)))))
        }
      }
      if(model %in% c("M1", "M1c")) {
        m <- rnorm(n = 1, mean = x[1]+x[2]*y, sd = x[3])
        if(trait.type == "D") {
          m <- 1/(1 + exp(-m))
          # m <- rbinom(n = 1, size = 1, prob = 1/(1 + exp(-m)))
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
                   level = "refSNP",
                   trait.type = gt.data$trait.type[t])

      # snp level
      ys[t,s,1,] <- apply(X = ext[, ps], MARGIN = 1,
                          FUN = getMuSnp, y = 1,
                          trait.type = gt.data$trait.type[t],
                          level = "refSNP")
      ys[t,s,2,] <- apply(X = ext[, ps], MARGIN = 1,
                          FUN = getMuSnp, y = -1,
                          trait.type = gt.data$trait.type[t],
                          level = "refSNP")


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
                        ppc.type = "h",
                        prediction.level = "refSNP",
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
                        ppc.type = "h",
                        prediction.level = "refSNP",
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


    # special case for strains
    if(model %in% c("M1", "M1c")) {
      for(k in 1:gt.data$Nk) {

        y.temp <- ysk[t,,k,]
        y.real <- mean(gt.data$Y[gt.data$K == k, t])
        yhat.hdi <- getHdi(vec = y.temp, hdi.level = hdi.level)
        es <- abs(y.temp-y.real)
        es.hdi <- getHdi(vec = es, hdi.level = hdi.level)
        row <- data.frame(t = t,
                          s = s,
                          x = unique(gt.data$X[gt.data$K == k, s]),
                          k = k,
                          i = NA,
                          ppc.type = "h",
                          prediction.level = "strain",
                          y.real.mean = y.real,
                          y.ppc.mean = mean(y.temp),
                          y.ppc.median = median(y.temp),
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
  }

  return (stats)
}



# Function:
# Posterior prediction from bottom to top
getPpcVertical <- function(ext,
                           gt.data,
                           model,
                           hdi.level) {

  # Function:
  # Use the most specific parameters to make predictions at different levels
  getPpcPosterior <- function(x, ext, gt.data, model) {

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
        # return(rbinom(n = 1, size = 1, prob = 1/(1 + exp(-(x[1]+x[2]*y)))))
        return(1/(1 + exp(-(x[1]+x[2]*y))))
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
                  s = x, model = model))
  }


  ppc.out <- lapply(X = 1:gt.data$Ns,
                    FUN = getPpcPosterior,
                    ext = ext,
                    gt.data = gt.data,
                    model = model)

  # compute statistics
  ys <- array(data = NA, dim = c(gt.data$Ntd+gt.data$Ntq,
                                 gt.data$Ns,
                                 gt.data$N,
                                 nrow(ext)))


  ei <- array(data = NA, dim = c(gt.data$Ntd+gt.data$Ntq,
                                 gt.data$N,
                                 nrow(ext)*gt.data$Ns))


  stats <- c()
  for(t in 1:(gt.data$Ntd+gt.data$Ntq)) {
    for(i in 1:gt.data$N) {
      for(s in 1:gt.data$Ns) {
        ys[t,s,i,] <- ppc.out[[s]][[t]][, i]
      }

      # stats at individual level
      y.ppc.hdi <- getHdi(vec = t(ys[t,,i,]), hdi.level = hdi.level)
      y.ppc.mean <- mean(t(ys[t,,i,]))
      y.ppc.median <- median(t(ys[t,,i,]))

      ei[t,i,] <- abs(as.vector(t(ys[t,,i,]))-gt.data$Y[i,t])
      y.error.mean <- mean(ei[t,i,])
      y.error.median <- median(ei[t,i,])
      y.error.hdi <- getHdi(vec = ei[t,i,], hdi.level = hdi.level)

      row <- data.frame(t = t,
                        s = NA,
                        x = NA,
                        k = NA,
                        i = i,
                        ppc.type = "v",
                        prediction.level = "individual",
                        y.real.mean = gt.data$Y[i, t],
                        y.ppc.mean = y.ppc.mean,
                        y.ppc.median = y.ppc.median,
                        y.ppc.L = y.ppc.hdi[1],
                        y.ppc.H = y.ppc.hdi[2],
                        y.error.mean = y.error.mean,
                        y.error.median = y.error.median,
                        y.error.L = y.error.hdi[1],
                        y.error.H = y.error.hdi[2],
                        model = model)
      stats <- rbind(stats, row)
    }


    # strains level
    for(k in 1:gt.data$Nk) {
      ks <- which(gt.data$K == k)

      # stats at strain level
      y.ppc.hdi <- getHdi(vec = ys[t,,ks,], hdi.level = hdi.level)
      y.ppc.mean <- mean(ys[t,,ks,])
      y.ppc.median <- median(ys[t,,ks,])

      e <- apply(X = t(ei[t,ks,]), MARGIN = 1, FUN = mean)
      y.error.mean <- mean(e)
      y.error.median <- median(e)
      y.error.hdi <- getHdi(vec = e, hdi.level = hdi.level)


      row <- data.frame(t = t,
                        s = NA,
                        x = NA,
                        k = k,
                        i = i,
                        ppc.type = "v",
                        prediction.level = "strain",
                        y.real.mean = mean(gt.data$Y[ks, t]),
                        y.ppc.mean = y.ppc.mean,
                        y.ppc.median = y.ppc.median,
                        y.ppc.L = y.ppc.hdi[1],
                        y.ppc.H = y.ppc.hdi[2],
                        y.error.mean = y.error.mean,
                        y.error.median = y.error.median,
                        y.error.L = y.error.hdi[1],
                        y.error.H = y.error.hdi[2],
                        model = model)
      stats <- rbind(stats, row)
    }



    # stats at snp level
    for(s in 1:gt.data$Ns) {
      for(x in c(1, -1)) {
        xs <- which(gt.data$X[, s] == x)
        if(length(xs) != 0) {


          # stats at SNP level
          y.temp <- ys[t,s,xs,]
          y.ppc.hdi <- getHdi(vec = y.temp, hdi.level = hdi.level)
          y.ppc.mean <- mean(y.temp)
          y.ppc.median <- median(y.temp)

          e <- apply(X = t(ei[t,xs,]), MARGIN = 1, FUN = mean)
          y.error.mean <- mean(e)
          y.error.median <- median(e)
          y.error.hdi <- getHdi(vec = e, hdi.level = hdi.level)


          row <- data.frame(t = t,
                            s = s,
                            x = x,
                            k = NA,
                            i = NA,
                            ppc.type = "v",
                            prediction.level = "refSNP",
                            y.real.mean = mean(gt.data$Y[xs, t]),
                            y.ppc.mean = y.ppc.mean,
                            y.ppc.median = y.ppc.median,
                            y.ppc.L = y.ppc.hdi[1],
                            y.ppc.H = y.ppc.hdi[2],
                            y.error.mean = y.error.mean,
                            y.error.median = y.error.median,
                            y.error.L = y.error.hdi[1],
                            y.error.H = y.error.hdi[2],
                            model = model)
          stats <- rbind(stats, row)
        }
      }

      if(s %% 50 == 0) {
        cat("SNP:", s, "/", gt.data$Ns, ',', sep = '')
      }
    }
  }
  return (stats)
}



