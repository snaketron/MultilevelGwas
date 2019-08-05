

getObservedData <- function(gt.data) {

  # Description:
  # x = snp id
  # t  = trait id
  # X = genotype matrix (cols = snps with -1 or +1)
  # Y = vector trait
  getS <- function(X, t, Xmatrix, Y) {
    i.pos <- which(Xmatrix[, X] == +1)
    i.neg <- which(Xmatrix[, X] == -1)

    mean.pos <- ifelse(test = length(i.pos) == 0, yes = NA, no = mean(Y[i.pos]))
    mean.neg <- ifelse(test = length(i.neg) == 0, yes = NA, no = mean(Y[i.neg]))
    return (rbind(data.frame(observed.mean = mean.pos, trait = t, X = 1,
                             snp = X, stringsAsFactors = FALSE),
                  data.frame(observed.mean = mean.neg, trait = t, X = -1,
                             snp = X, stringsAsFactors = FALSE)))
  }

  # Description:
  # t  = trait id
  # K = strains
  # Y = vector trait
  getK <- function(t, K, Y) {

    k.data <- data.frame(Y = Y, K = K, stringsAsFactors = FALSE)
    mean.k <- aggregate(formula = Y~K, data = k.data, FUN = mean)
    mean.k$observed.mean <- mean.k$Y
    mean.k$trait <- t
    mean.k$strain <- mean.k$K
    mean.k$K <- NULL
    mean.k$Y <- NULL

    return (mean.k)
  }

  # Description:
  # t = trait id
  # K = strains
  # Y = vector trait
  getI <- function(t, K, Y) {
    mean.i <- data.frame(individual = 1:length(Y),
                         trait = t,
                         strain = K,
                         observed.mean = Y,
                         stringsAsFactors = FALSE)

    return (mean.i)
  }


  out.s <- c()
  out.k <- c()
  out.i <- c()
  for(t in 1:(gt.data$Ntd+gt.data$Ntq)) {
    os <- lapply(X = 1:gt.data$Ns, FUN = getS, t = t,
                 Xmatrix = gt.data$X, Y = gt.data$Y[, t])
    os <- do.call(rbind, os)
    out.s <- rbind(out.s, os)
    rm(os)

    out.k <- rbind(out.k, getK(t = t, K = gt.data$K, Y = gt.data$Y[, t]))

    out.i <- rbind(out.i, getI(t = t, K = gt.data$K, Y = gt.data$Y[, t]))
  }

  return (list(out.s = out.s, out.k = out.k, out.i = out.i))
}
