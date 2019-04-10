
# Description:
# Combine data from Bayesian inference and statistical learning
# betas (b) = matrix posterior (cols = snps, rows = posterior elements P)
# kappa = matrix posterior (cols = snps, rows = posterior elements P)
getParetoRanks <- function(p, s, hdi.level, model) {

  # Description:
  # Computes HDI given a vector, taken "Doing Bayesian Analysis"
  getHdiFronts <- function(x, hdi.level) {
    sortedPts <- sort(x)
    ciIdxInc <- floor(hdi.level * length(sortedPts))
    nCIs = length(sortedPts) - ciIdxInc
    ciWidth = rep(0 , nCIs)
    for (i in 1:nCIs) {
      ciWidth[i] = sortedPts[i + ciIdxInc] - sortedPts[i]
    }
    HDImin = sortedPts[which.min(ciWidth)]
    HDImax = sortedPts[which.min(ciWidth) + ciIdxInc]
    HDIlim = c(HDImin, HDImax)
    return(HDIlim)
  }


  getParetoRanksTrait <- function(b, k, hdi.level) {

    # check
    if(nrow(b) != nrow(k)) {
      stop("nrow b must be equal to nrow k in length")
    }
    if(ncol(b) != ncol(k)) {
      stop("ncol b must be equal to ncol k in length")
    }


    # collect front statistics
    fronts <- matrix(data = NA, nrow = nrow(k), ncol = ncol(k))
    for(i in 1:nrow(k)) {
      d <- data.frame(beta.mean = b[i, ], k = k[i, ], s = 1:ncol(k))
      p <- rPref::high(abs(d$beta.mean))*rPref::high(d$k)
      f <- rPref::psel(df = d, pref = p, top = nrow(d))
      f$rank <- f[, ".level"]
      fronts[i, ] <- f$rank[order(f$s, decreasing = FALSE)]
    }


    # statistics
    front.mean <- colMeans(x = fronts)
    front.median <- apply(X = fronts, MARGIN = 2, FUN = median)
    front.hdi <- t(apply(X = fronts, MARGIN = 2,
                         FUN = getHdiFronts,
                         hdi.level = hdi.level))


    # front summary
    front.summary <- data.frame(site = 1:length(front.mean),
                                front.mean = front.mean,
                                front.median = front.median,
                                front.L = front.hdi[, 1],
                                front.H = front.hdi[, 2])


    return (front.summary)
  }


  if(model == "M0" | model == "M0c") {
    pars <- c("beta")
  }
  if(model == "M1" | model == "M1c") {
    pars <- c("mu_beta")
  }


  b <- rstan::extract(p, pars = pars)


  if(length(s) != length(b$beta[1,,1])) {
    stop("Beta and CA not equal across traits")
  }


  # subsample
  b.N <- length(b$beta[, 1, 1])
  k.N <- nrow(s[[1]]$p.ka)
  min.N <- 500 # default minimum
  sample.N <- min(min.N, min(b.N, k.N))

  if(sample.N < k.N) {
    sample.k <- sample(x = 1:k.N,
                       size = sample.N,
                       replace = TRUE)
  } else {
    sample.k <- 1:k.N
  }

  if(sample.N < b.N) {
    sample.b <- sample(x = 1:b.N,
                       size = sample.N,
                       replace = TRUE)
  } else {
    sample.b <- 1:b.N
  }
  rm(b.N, k.N, min.N, sample.N)


  # collect ranks
  ranks <- c()
  ts <- length(s)
  for(t in 1:ts) {
    r <- getParetoRanksTrait(b = b$beta[sample.b, t, ],
                             k = s[[t]]$p.ka[sample.k, ],
                             hdi.level = hdi.level)
    r$trait <- t
    ranks <- rbind(ranks, r)
  }

  return (ranks)
}

