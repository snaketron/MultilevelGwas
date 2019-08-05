
# Description:
# Summarizes the observed data into one used for PPC analysis.
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



# Description:
# mcmc.warmup: number of warmup steps used
# par: complete name of the parameter: Yhat_individual, Yhat_snp, or
# Yhat_strain (only in M1/M1c)
# files: list of file paths
getPpcFromSampling <- function(files, mcmc.warmup, par) {

  createCustomAwk <- function(is) {
    str <- paste('{print ', paste(paste("$", is, sep = ''),
                                  collapse = ','), "}", sep = '')
    return (str)
  }

  if(length(par) != 1) {
    stop("par must be one of 'Yhat_individual', 'Yhat_strain' or 'Yhat_snp'")
  }

  if(is.character(par) == FALSE) {
    stop("par must be one of 'Yhat_individual', 'Yhat_strain' or 'Yhat_snp'")
  }

  if(all(par %in% c("Yhat_individual", "Yhat_strain", "Yhat_snp")) == FALSE) {
    stop("par must be one of 'Yhat_individual', 'Yhat_strain' or 'Yhat_snp'")
  }

  data <- vector(mode = "list", length = length(files))

  for(f in 1:length(files)) {
    if(file.exists(files[f]) == FALSE) {
      stop(paste("File ", files[f], " does not exist \n", sep = ''))
    }

    awk.path.pipe <- pipe(paste("gawk -F , 'FNR == 26'", files[f], sep = ' '))
    cols <- readLines(con = awk.path.pipe)
    close(con = awk.path.pipe)
    cols <- unlist(strsplit(x = cols, split = ','))
    is <- which(regexpr(pattern = par, text = cols, ignore.case = TRUE) != -1)

    if(length(is) == 0) {
      stop("no such parameter found in sampling file.")
    }

    awk <- createCustomAwk(is)
    writeLines(text = awk, con = "awk_temp.awk")

    awk.path.pipe <- pipe(paste("gawk -F , -f awk_temp.awk",
                                files[f], "> temp.csv", sep = ' '))
    d <- readLines(con = awk.path.pipe)
    close(con = awk.path.pipe)

    awk.path.pipe <- pipe("gawk -F , '(NR > 25)' temp.csv > temp2.csv")
    d <- readLines(con = awk.path.pipe)
    close(con = awk.path.pipe)

    # here take only samples (-warmup)
    d <- read.table(file = "temp2.csv", header = TRUE, as.is = TRUE, sep = ' ')
    d <- d[complete.cases(d), ]
    d <- d[-(1:mcmc.warmup), ]
    data[[f]] <- d

    # do some cleanup before next file
    file.remove(c("awk_temp.awk", "temp.csv", "temp2.csv"))
  }

  # finally summarize data
  rm(d)
  data <- do.call(rbind, data)
  y.mean <- apply(X = data, MARGIN = 2, FUN = mean)
  y.median <- apply(X = data, MARGIN = 2, FUN = median)
  y.hdi <- apply(X = data, MARGIN = 2, FUN = getHdi, hdi.level = 0.95)
  y.hdi <- t(y.hdi)
  y.data <- data.frame(par = names(y.mean), ppc.mean = y.mean,
                       ppc.median = y.median, ppc.L = y.hdi[, 1],
                       ppc.H = y.hdi[, 2], stringsAsFactors = FALSE)
  par.data <- strsplit(x = y.data$par, split = "\\.")
  par.data <- do.call(rbind, par.data)


  if(par == "Yhat_individual") {
    y.data$par.type <- par.data[, 1]
    y.data$trait <- as.numeric(par.data[, 2])
    y.data$individual <- as.numeric(par.data[, 3])
    y.data$snp <- as.numeric(par.data[, 4])
  }

  if(par == "Yhat_strain") {
    y.data$par.type <- par.data[, 1]
    y.data$trait <- as.numeric(par.data[, 2])
    y.data$strain <- as.numeric(par.data[, 3])
    y.data$snp <- as.numeric(par.data[, 4])
  }

  if(par == "Yhat_snp") {
    y.data$par.type <- par.data[, 1]
    y.data$trait <- as.numeric(par.data[, 2])
    y.data$X <- ifelse(test = as.numeric(par.data[, 3]) == 1, yes = 1, no = -1)
    y.data$snp <- as.numeric(par.data[, 4])
  }

  y.data$par <- NULL

  return (y.data)
}


# Description:
# Horizontal posterior predictions
# Uses: getObservedData
# Uses: getPpcFromSampling
getPpc <- function(gt.data, sampling.files,
                   mcmc.warmup, model) {

  # 1.
  o <- getObservedData(gt.data = gt.data)

  # 2.
  yhat.snp <- getPpcFromSampling(files = sampling.files,
                                 mcmc.warmup = mcmc.warmup,
                                 par = "Yhat_snp")

  # 3.
  yhat.individual <- getPpcFromSampling(files = sampling.files,
                                        mcmc.warmup = mcmc.warmup,
                                        par = "Yhat_individual")

  # 4.
  ppc.snp <- merge(x = o$out.s, y = yhat.snp,
                   by = c("trait", "snp", "X"))
  ppc.snp$error.mean <- abs(ppc.snp$ppc.mean-ppc.snp$observed.mean)


  # 5.
  ppc.individual <- merge(x = o$out.i, y = yhat.individual,
                          by = c("trait", "individual"))
  ppc.individual$error.mean <- abs(ppc.individual$ppc.mean-
                                     ppc.individual$observed.mean)

  # 6.
  ppc.strain <- NA
  if(model %in% c("M1", "M1c")) {
    yhat.strain <- getPpcFromSampling(files = sampling.files,
                                      mcmc.warmup = mcmc.warmup,
                                      par = "Yhat_strain")

    ppc.strain <- merge(x = o$out.k, y = yhat.strain, by = c("trait", "strain"))
    ppc.strain$error.mean <- abs(ppc.strain$ppc.mean-ppc.strain$observed.mean)
  }

  # 7.
  return(list(ppc.individual = ppc.individual,
              ppc.strain = ppc.strain,
              ppc.snp = ppc.snp))
}
