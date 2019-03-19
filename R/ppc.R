
# Description:
# Computes posterior prediction
getPpc <- function(posterior, genphen.data, hdi.level) {

  ppc.Q <- function(p, x) {
    return(p[1]+p[2]*x+p[3]*stats::rt(n = 1, df = p[4]))
  }

  ppc.mean.Q <- function(p, x) {
    return(mean(replicate(n = 1, expr = p[1]+p[2]*x)))
  }

  ppc.D <- function(p, x) {
    o <- p[1]+p[2]*x
    o <- exp(o)/(1+exp(o))
    return(o)
  }

  ppc.mean.D <- function(p, x) {
    o <- replicate(n = 1, expr = p[1]+p[2]*x)
    o <- 1/(1 + exp(-(o)))
    return(mean(o))
  }

  getParamKeys <- function(p, i, genphen.data) {
    alpha.key <- paste("alpha.", p, ".", i,  sep = '')
    beta.key <- paste("beta.", p, ".", i,  sep = '')

    if(genphen.data$Ntq == 1) {
      sigma.key <- "sigma"
      nu.key <- "nu"
    }
    else {
      sigma.key <- paste("sigma.", p, sep = '')
      nu.key <- paste("nu.", p, sep = '')
    }

    p <- c(alpha.key, beta.key, sigma.key, nu.key)
    return (p)
  }



  # extract posterior
  posterior <- data.frame(rstan::extract(object = posterior))
  max.nrow <- ifelse(test = 1000 > nrow(posterior),
                     yes = nrow(posterior), no = 1000)
  posterior <- posterior[sample(x = 1:nrow(posterior),
                                size = max.nrow,
                                replace = TRUE), ]

  xmap <- genphen.data$xmap
  ppc <- vector(mode = "list", length = length(genphen.data$phenotype.type))

  for(p in 1:length(ppc)) {

    ppc.out <- vector(mode = "list", length = nrow(xmap))
    ppc.vec <- numeric(length = 8)

    for(i in 1:nrow(xmap)) {
      p.keys <- getParamKeys(p = p, i = i, genphen.data = genphen.data)

      # REF
      ref <- which(genphen.data$X[, i] == 1)
      ref.p <- NA
      if(length(ref) != 0) {
        if(genphen.data$phenotype.type[p] == "Q") {
          ref.p <- genphen.data$Y[ref, p]
          ppc.vec[1] <- mean(ref.p)
          ppc.data <- apply(X = posterior[, p.keys], MARGIN = 1,
                            FUN = ppc.Q, x = 1)
          ppc.hdi <- getHdi(vec = ppc.data, hdi.level = hdi.level)
          ppc.vec[2] <- mean(ppc.data)
          ppc.vec[3] <- ppc.hdi[1]
          ppc.vec[4] <- ppc.hdi[2]
        }
        if(genphen.data$phenotype.type[p] == "D") {
          ref.p <- genphen.data$Y[ref, p]
          ppc.vec[1] <- sum(ref.p == 1)/length(ref.p)
          ppc.data <- apply(X = posterior[, p.keys[1:2]], MARGIN = 1,
                            FUN = ppc.D, x = -1)
          ppc.hdi <- getHdi(vec = ppc.data, hdi.level = hdi.level)
          ppc.vec[2] <- mean(ppc.data)
          ppc.vec[3] <- ppc.hdi[1]
          ppc.vec[4] <- ppc.hdi[2]
        }
      }

      # ALT
      alt <- which(genphen.data$X[, i] == -1)
      alt.p <- NA
      if(length(alt) != 0) {
        if(genphen.data$phenotype.type[p] == "Q") {
          alt.p <- genphen.data$Y[alt, p]
          ppc.vec[5] <- mean(alt.p)
          ppc.data <- apply(X = posterior[, p.keys], MARGIN = 1,
                            FUN = ppc.Q, x = -1)
          ppc.hdi <- getHdi(vec = ppc.data, hdi.level = hdi.level)
          ppc.vec[6] <- mean(ppc.data)
          ppc.vec[7] <- ppc.hdi[1]
          ppc.vec[8] <- ppc.hdi[2]
        }
        if(genphen.data$phenotype.type[p] == "D") {
          alt.p <- genphen.data$Y[alt, p]
          ppc.vec[5] <- sum(alt.p == 1)/length(alt.p)
          ppc.data <- apply(X = posterior[, p.keys[1:2]],
                            MARGIN = 1, FUN = ppc.D, x = -1)
          ppc.hdi <- getHdi(vec = ppc.data, hdi.level = hdi.level)
          ppc.vec[6] <- mean(ppc.data)
          ppc.vec[7] <- ppc.hdi[1]
          ppc.vec[8] <- ppc.hdi[2]
        }
      }

      ppc.out[[i]] <- list(site = xmap$site[i],
                           ref = xmap$ref[i],
                           alt = xmap$alt[i],
                           ref.Y = ref.p,
                           alt.Y = alt.p,
                           ref.real = ppc.vec[1],
                           ref.ppc = ppc.vec[2],
                           ref.ppc.hdi.low = ppc.vec[3],
                           ref.ppc.hdi.high = ppc.vec[4],
                           alt.real = ppc.vec[5],
                           alt.ppc = ppc.vec[6],
                           alt.ppc.hdi.low = ppc.vec[7],
                           alt.ppc.hdi.high = ppc.vec[8])
    }

    ppc[[p]] <- ppc.out
  }

  return (ppc)
}



