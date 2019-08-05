

# Description:
# Parse input data and format it for stan and statistical learning
getStanData <- function(genotype,
                        traits,
                        trait.type,
                        strains) {


  # Description:
  # Get genotype-traits data
  getFormattedGenphen <- function(genotype,
                                  traits,
                                  trait.type,
                                  strains) {


    # convert AAMultipleAlignment to matrix if needed
    genotype <- convertMsaToGenotype(genotype = genotype)


    # if vector genotype => matrix genotype
    if(is.vector(genotype)) {
      genotype <- matrix(data = genotype, ncol = 1)
    }

    # if vector genotype => matrix genotype
    if(is.vector(traits)) {
      traits <- matrix(data = traits, ncol = 1)
    }

    traits <- data.frame(traits)


    # TODO: check
    if(sum(trait.type == "D") != 0) {
      d <- which(trait.type == "D")
      for(i in 1:length(d)) {
        if(all(traits[, d] %in% c(1, 0)) == FALSE) {
          # mapping 1st element to 1, 2nd to 0
          u <- unique(traits[, d])
          traits[traits[, d] == u[1], d] <- "1"
          traits[traits[, d] == u[2], d] <- "0"
          cat("Mapping dichotomous traits:", i,
              "(", u[1], "->1,", u[2], "->0) \n")
        }
        traits[, d] <- as.factor(as.character(traits[, d]))
      }
    }

    # return
    return (list(genotype = genotype,
                 traits = traits,
                 trait.type = trait.type,
                 strains = strains))
  }


  # Description:
  # Convert genphen data to stan data
  getStanFormat <- function(f.data) {

    # Split SNP into pairs of genotypes
    getX <- function(x, X.setup) {
      xs <- unique(x)

      # Is X setup with 1s and -1s?
      # * Yes = leave it as it is
      # * No = transform
      if(X.setup == FALSE) {
        X <- numeric(length = length(x))
        xs.dir <- sample(x = c(1, -1),
                         size = 2,
                         replace = FALSE)

        # map genotypes to 1s and -1s
        for(i in 1:length(xs)) {
          X[x == xs[i]] <- xs.dir[i]
        }
      }
      else {
        X <- x
      }



      xs <- c(xs, NA)

      ref <- unique(x[X == 1])
      alt <- unique(x[X == -1])
      if(length(ref) == 0) {
        ref <- NA
      }
      if(length(alt) == 0) {
        alt <- NA
      }

      # xmap
      xmap <- data.frame(ref = ref,
                         alt = alt,
                         refN = sum(X == 1),
                         altN = sum(X == -1))

      # return
      return(list(X = X, xmap = xmap))
    }

    # X.data
    getXData <- function(x) {
      return (x$X)
    }

    # X.map
    getXMap <- function(x) {
      return (x$xmap)
    }


    # make huge matrix with predictors
    X <- f.data$genotype


    X.setup <- all(X == 1 | X == -1)
    X <- apply(X = X, MARGIN = 2, FUN = getX, X.setup = X.setup)
    for(i in 1:length(X)) {
      X[[i]]$xmap$site <- i
    }
    xmap <- do.call(rbind, lapply(X = X, FUN = getXMap))
    X <- do.call(cbind, lapply(X = X, FUN = getXData))


    if(sum(f.data$trait.type == "Q") != 0) {
      Yq <- as.matrix(f.data$traits[, f.data$trait.type == "Q"])
    }
    else {
      Yq <- matrix(data = 0, ncol = 0, nrow = nrow(X))
    }

    if(sum(f.data$trait.type == "D") != 0) {
      Yd <- as.matrix(f.data$traits[, f.data$trait.type == "D"])
    }
    else {
      Yd <- matrix(data = 0, ncol = 0, nrow = nrow(X))
    }


    # extra conversion
    temp.Yd <- c()
    if(ncol(Yd) > 0) {
      for(i in 1:ncol(Yd)) {
        temp.Yd <- cbind(temp.Yd, as.numeric(Yd[, i]))
      }
      Yd <- temp.Yd
      rm(temp.Yd)
    }

    # set X to numeric
    class(X) <- "numeric"


    K <- as.numeric(as.factor(f.data$strains))
    Xk <- X[which(duplicated(K) == FALSE), ]


    s <- list(X = X,
              Y = cbind(Yq, Yd),
              K = K,
              Yq = Yq,
              Yd = Yd,
              N = nrow(X),
              Ns = ncol(X),
              Nk = length(unique(f.data$strains)),
              Ntq = ncol(Yq),
              Ntd = ncol(Yd),
              Xk = Xk,
              xmap = xmap,
              Korg = f.data$strains,
              strains = f.data$strains,
              genotype = f.data$genotype,
              trait.type = f.data$trait.type)

    return (s)
  }



  # genphen data
  f.data <- getFormattedGenphen(genotype = genotype,
                                traits = traits,
                                trait.type = trait.type,
                                strains = strains)


  # stan data
  stan.data <- getStanFormat(f.data = f.data)


  # return
  return (stan.data)
}




# Function:
# Given model.name, fetch appropriate stan model from
# extdata and compile
getStanModel <- function(model.name) {
  rstan::rstan_options(auto_write = TRUE)
  model.file <- system.file("extdata", paste(model.name, '.stan', sep = ''),
                            package = "MultilevelGwas")
  model <- rstan::stan_model(file = model.file, auto_write = TRUE)
  return (model)
}





# Function:
# Given model.name, fetch appropriate stan model
getStanModelDebugReal <- function(model.name,
                              comparison = TRUE) {
  if(comparison) {

    # M0
    if(model.name == "M0") {
      stan.model <- rstan::stan_model(file = "src/stan_files/M0_loglik.stan")
    }
    # M0c
    if(model.name == "M0c") {
      stan.model <- rstan::stan_model(file = "src/stan_files/M0c_loglik.stan")
    }

    # M1
    if(model.name == "M1") {
      stan.model <- rstan::stan_model(file = "src/stan_files/M1_loglik.stan")
    }
    # M1c
    if(model.name == "M1c") {
      stan.model <- rstan::stan_model(file = "src/stan_files/M1c_loglik.stan")
    }

    # M2
    if(model.name == "M2") {
      stan.model <- rstan::stan_model(file = "src/stan_files/M2_loglik.stan")
    }
    # M2c
    if(model.name == "M2c") {
      stan.model <- rstan::stan_model(file = "src/stan_files/M2c_loglik.stan")
    }

  }
  else {

    # M0
    if(model.name == "M0") {
      stan.model <- rstan::stan_model(file = "src/stan_files/M0.stan")
    }
    # M0c
    if(model.name == "M0c") {
      stan.model <- rstan::stan_model(file = "src/stan_files/M0c.stan")
    }

    # M1
    if(model.name == "M1") {
      stan.model <- rstan::stan_model(file = "src/stan_files/M1.stan")
    }
    # M1c
    if(model.name == "M1c") {
      stan.model <- rstan::stan_model(file = "src/stan_files/M1c.stan")
    }

    # M2
    if(model.name == "M2") {
      stan.model <- rstan::stan_model(file = "src/stan_files/M2.stan")
    }
    # M2c
    if(model.name == "M2c") {
      stan.model <- rstan::stan_model(file = "src/stan_files/M2c.stan")
    }

  }


  return (stan.model)
}



# Function:
# Given model.name, fetch appropriate stan model
getStanModelDebug <- function(model.name,
                              comparison = TRUE) {

  dir <- ''

  if(comparison) {
    model.name <- paste("src/stan_files/", dir, model.name, "_loglik.stan", sep = '')
    stan.model <- rstan::stan_model(file = model.name)
  }
  else {
    model.name <- paste("src/stan_files/", dir, model.name, ".stan", sep = '')
    stan.model <- rstan::stan_model(file = model.name)
  }

  return (stan.model)
}


# Description:
# These are the most necessary parameters (for the user). The remaining are
# anyhow exported in the sampling files.
getStanModelPars <- function(model.name) {

  # M0
  if(model.name == "M0") {
    pars <- c("alpha_trait", "sigma", "beta_snp",
              "nu", "nu_help")
  }
  if(model.name == "M0c") {
    pars <- c("alpha_trait", "sigma", "beta_snp",
              "nu", "nu_help", "rho")
  }
  # M1
  if(model.name == "M1") {
    pars <- c("alpha_trait", "sigma", "beta_snp",
              "sigma_snp", "nu", "nu_help")
  }
  if(model.name == "M1c") {
    pars <- c("alpha_trait", "sigma", "beta_snp",
              "sigma_snp", "nu", "nu_help", "rho")
  }

  return (pars)
}






# Description:
# Check the input arguments of the main run, stop if violated.
checkInput <- function(genotype,
                       traits,
                       trait.type,
                       strains,
                       model,
                       mcmc.chains,
                       mcmc.steps,
                       mcmc.warmup,
                       mcmc.cores,
                       hdi.level,
                       adapt.delta,
                       max.treedepth) {


  if(is.null(genotype) | missing(genotype) |
     is.null(traits) | missing(traits) |
     is.null(trait.type) | missing(trait.type) |
     is.null(strains) | missing(strains) |
     is.null(model) | missing(model) |
     is.null(mcmc.chains) | missing(mcmc.chains) |
     is.null(mcmc.steps) | missing(mcmc.steps) |
     is.null(mcmc.warmup) | missing(mcmc.warmup) |
     is.null(mcmc.cores) | missing(mcmc.cores) |
     is.null(adapt.delta) | missing(adapt.delta) |
     is.null(max.treedepth) | missing(max.treedepth) |
     is.null(hdi.level) | missing(hdi.level)) {
    stop("arguments must be non-NULL/specified")
  }

  checkGenotypeTrait(genotype = genotype, traits = traits,
                     trait.type = trait.type)
  checkTraitValidity(traits = traits, trait.type = trait.type)
  checkStrains(strains = strains, traits = traits)
  checkModel(model = model)
  checkMcmcSteps(mcmc.steps = mcmc.steps)
  checkMcmcWarmup(mcmc.warmup = mcmc.warmup)
  checkMcmcChains(mcmc.chains = mcmc.chains)
  checkMcmcCores(mcmc.cores = mcmc.cores)
  checkHdi(hdi.level = hdi.level)
  checkAdaptDelta(adapt.delta = adapt.delta)
  checkMaxTreedepth(max.treedepth = max.treedepth)
}




# Function:
# Provided the input arguments, this function checks their validity. It
# stops the execution if a problem is encountered and prints out warnings.
checkInputPhyloBias <- function(input.kinship.matrix,
                                genotype) {


  checkGenotype <- function(genotype) {
    if(is.null(attr(genotype, "class")) == FALSE) {
      if(!attr(genotype, "class") %in% c("AAMultipleAlignment",
                                         "DNAMultipleAlignment")) {
        stop("genotype can be one of the following structures: matrix or
             data.frame or DNAMultipleAlignment/AAMultipleAlignment")
      }
      else {
        temp <- as.matrix(genotype)
        if(nrow(temp) < 2 | ncol(temp) == 0) {
          stop("the genotypes cannot have less than two observations the or
               number of genotypes cannot be 0.")
        }
      }
    }
    else {
      if(is.vector(genotype)) {
        genotype <- matrix(data = genotype, ncol = 1)
      }

      if(nrow(genotype) < 2 | ncol(genotype) == 0) {
        stop("the genotypes cannot have less than two observations or the
             number of genotypes cannot be 0.")
      }

      if(!is.matrix(genotype) & !is.data.frame(genotype)) {
        stop("genotype can be one of the following structures: matrix or
             data.frame or DNAMultipleAlignment/AAMultipleAlignment")
      }

      if(typeof(genotype) != "character") {
        stop("if it is structured in matrix/data.frame the genotype must
             be of character type.")
      }
    }
  }


  checkKinship <- function(input.kinship.matrix) {

    if(is.matrix(input.kinship.matrix) == FALSE) {
      stop("precomputed kinship matrix must be a numeric matrix.")
    }

    if(is.numeric(input.kinship.matrix) == FALSE) {
      stop("precomputed kinship matrix must be a numeric matrix.")
    }

    if(nrow(input.kinship.matrix) != ncol(input.kinship.matrix)) {
      stop("precomputed kinship matrix must be NxN numeric matrix.")
    }

    if(nrow(input.kinship.matrix) <= 0) {
      stop("at least two individuals needed for the analysis.")
    }
  }


  if((is.null(input.kinship.matrix) | missing(input.kinship.matrix))
     & (is.null(genotype) | missing(genotype))) {
    stop("arguments must be non-NULL/specified")
  }

  if(is.null(input.kinship.matrix) | missing(input.kinship.matrix)) {
    if(is.null(genotype) | missing(genotype)) {
      stop("arguments must be non-NULL/specified")
    }
  }
  else {
    checkKinship(input.kinship.matrix = input.kinship.matrix)
  }
  checkGenotype(genotype = genotype)
}




# Function:
# If an object of type DNAMultipleAlignment
convertMsaToGenotype <- function(genotype) {
  if(is.null(attr(genotype, "class")) == FALSE) {
    genotype <- as.matrix(genotype)
  }
  return (genotype)
}




# Function:
# Combine data from Bayesian inference and statistical learning
# p = posterior
getScores <- function(glm, model, hdi.level, gt.data) {

  pars <- ifelse(test = model %in% c("M0", "M0c"),
                 yes = "beta", no = "mu_beta")

  d <- summary(glm, probs = c(0.5, (1-hdi.level)/2, 1-(1-hdi.level)/2),
               pars = pars)$summary
  d <- data.frame(d)
  d$par <- rownames(d)
  d$par <- gsub(pattern = "mu_beta|beta|\\[|\\]",
                replacement = '', x = d$par)
  key <- do.call(rbind, strsplit(x = d$par, split = '\\,'))
  d$trait <- as.numeric(key[, 1])
  d$site <- as.numeric(key[, 2])
  d$par <- NULL

  colnames(d) <- c("beta.mean", "beta.mean.se", "beta.sd",
                   "beta.median", "beta.L", "beta.H",
                   "neff", "rhat", "trait", "site")


  xmap <- gt.data$xmap
  scores <- merge(x = xmap, y = d, by = "site")
  scores <- scores[order(scores$trait, scores$site, decreasing = FALSE), ]

  # pmax
  pmax <- getPmax(glm)
  scores <- merge(x = scores, y = pmax, by = c("site", "trait"))

  return (scores)
}



getPmax <- function(glm) {
  pmax <- c()

  if(glm@model_name %in% c("M0", "M0c")) {
    e <- rstan::extract(object = glm, pars = c("beta"))
    e.dim <- dim(e$beta)
    Ns <- e.dim[3]
    Nt <- e.dim[2]
    Np <- e.dim[1]

    for(t in 1:Nt) {
      for(s in 1:Ns) {
        row <- data.frame(trait = t, site = s,
                          pmax = max(sum(e$beta[,t,s] > 0)/Np,
                                     sum(e$beta[,t,s] < 0)/Np),
                          stringsAsFactors = FALSE)
        pmax <- rbind(pmax, row)
      }
    }
    return (pmax)
  }

  if(glm@model_name %in% c("M1", "M1c")) {
    e <- rstan::extract(object = glm, pars = c("mu_beta"))
    e.dim <- dim(e$mu_beta)
    Ns <- e.dim[3]
    Nt <- e.dim[2]
    Np <- e.dim[1]

    for(t in 1:Nt) {
      for(s in 1:Ns) {
        row <- data.frame(trait = t, site = s,
                          pmax = max(sum(e$mu_beta[,t,s] > 0)/Np,
                                     sum(e$mu_beta[,t,s] < 0)/Np),
                          stringsAsFactors = FALSE)
        pmax <- rbind(pmax, row)
      }
    }
    return (pmax)
  }
}



# Function:
# Computes HDI given a vector, taken "Doing Bayesian Analysis"
getHdi <- function(vec, hdi.level) {
  sortedPts <- sort(vec)
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

    awk.path <- paste("gawk -F , 'FNR == 26'", files[f], sep = ' ')
    cols <- readLines(con = pipe(awk.path))
    cols <- unlist(strsplit(x = cols, split = ','))
    is <- which(regexpr(pattern = par, text = cols, ignore.case = TRUE) != -1)

    awk <- createCustomAwk(is)
    writeLines(text = awk, con = "awk_temp.awk")

    awk.path <- paste("gawk -F , -f awk_temp.awk",
                      files[f], "> temp.csv", sep = ' ')
    d <- readLines(con = pipe(awk.path))
    d <- readLines(con = pipe("gawk -F , '(NR > 25)' temp.csv > temp2.csv"))

    # here take only samples (-warmup)
    d <- read.table(file = "temp2.csv", header = TRUE, as.is = TRUE, sep = ' ')
    d <- d[complete.cases(d), ]
    d <- d[-(1:mcmc.warmup), ]
    data[[f]] <- d

    # close active connections
    # close(con = pipe(awk.path))
    # close(con = pipe("gawk -F , '(NR > 25)' temp.csv > temp2.csv"))
    closeAllConnections()

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
getPpc <- function() {

}


