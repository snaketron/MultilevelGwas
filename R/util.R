# Function:
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
    getX <- function(x) {
      xs <- unique(x)
      X <- numeric(length = length(x))
      xs.dir <- sample(x = c(1, -1), size = 2, replace = FALSE)
      # xs.dir <- c(1, -1)

      # map genotypes to 1s and -1s
      for(i in 1:length(xs)) {
        X[x == xs[i]] <- xs.dir[i]
      }

      xs <- c(xs, NA)

      # xmap
      xmap <- data.frame(ref = xs[1], alt = xs[2],
                         refN = sum(X == xs[1]),
                         altN = sum(X == xs[2]))

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
    X <- apply(X = X, MARGIN = 2, FUN = getX)
    for(i in 1:length(X)) {
      X[[i]]$xmap$site <- i
    }
    xmap <- do.call(rbind, lapply(X = X, FUN = getXMap))
    X <- do.call(cbind, lapply(X = X, FUN = getXData))


    if(sum(trait.type == "Q") != 0) {
      Yq <- as.matrix(f.data$traits[, f.data$trait.type == "Q"])
    }
    else {
      Yq <- matrix(data = 0, ncol = 0, nrow = nrow(X))
    }

    if(sum(trait.type == "D") != 0) {
      Yd <- as.matrix(f.data$traits[, f.data$trait.type == "D"])
    }
    else {
      Yd <- matrix(data = 0, ncol = 0, nrow = nrow(X))
    }

    K <- as.numeric(as.factor(strains))
    s <- list(X = X,
              Y = f.data$traits,
              K = K,
              Yq = Yq,
              Yd = Yd,
              N = nrow(X),
              Ns = ncol(X),
              Nk = length(unique(strains)),
              Ntq = ncol(Yq),
              Ntd = ncol(Yd),
              xmap = xmap,
              strains = strains,
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
# Parse input data and format it for stan and statistical learning
getPhyloData <- function(genotype) {

  # Description:
  # Convert genphen data to stan data
  getPhyloFormat <- function(genotype) {


    # Split SNP into pairs of genotypes
    getSplitX <- function(x) {
      ux <- unique(x)
      ns <- length(ux)
      nc <- max(c(choose(n = ns, k = 2), 1))
      x.split <- matrix(data = 0, nrow = length(x), ncol = nc)

      if(ns == 1) {
        x.split[, 1] <- 1
        x.map <- data.frame(ref = ux[1], alt = NA,
                            refN = sum(x == ux[1]),
                            altN = 0,
                            stringsAsFactors = FALSE)
      }
      else {
        x.map <- c()
        counter <- 1
        for(i in 1:(ns-1)) {
          for(j in (i+1):ns) {
            x.split[x == ux[i], counter] <- 1
            x.split[x == ux[j], counter] <- -1
            counter <- counter + 1
            x.map <- rbind(x.map, data.frame(ref = ux[i], alt = ux[j],
                                             refN = sum(x == ux[i]),
                                             altN = sum(x == ux[j]),
                                             stringsAsFactors = FALSE))
          }
        }
      }

      # return
      return(list(x.split = x.split,
                  x.map = x.map))
    }


    # X.data
    getXData <- function(x) {
      return (x$x.split)
    }


    # X.map
    getXMap <- function(x) {
      return (x$x.map)
    }


    # convert AAMultipleAlignment to matrix if needed
    genotype <- convertMsaToGenotype(genotype = genotype)


    # if vector genotype => matrix genotype
    if(is.vector(genotype)) {
      genotype <- matrix(data = genotype, ncol = 1)
    }

    # make huge matrix with predictors
    X <- genotype
    colnames(X) <- 1:ncol(X)
    x.data <- apply(X = X, MARGIN = 2, FUN = getSplitX)
    for(i in 1:length(x.data)) {
      x.data[[i]]$x.map$site <- i
    }
    X <- do.call(cbind, lapply(X = x.data, FUN = getXData))
    x.map <- do.call(rbind, lapply(X = x.data, FUN = getXMap))

    return (x.map)
  }

  # phylo format data
  phylo.data <- getPhyloFormat(genotype = genotype)

  # return
  return (phylo.data)
}





# Function:
# Given model.name, fetch appropriate stan model
getStanModel <- function(model.name, comparison.mode) {
  # M0
  if(model.name == "M0") {
    if(comparison.mode == TRUE) {
      model <- stanmodels$M0_loglik
    }
    else {
      model <- stanmodels$M0
    }
  }

  # M0c
  if(model.name == "M0c") {
    if(comparison.mode == TRUE) {
      model <- stanmodels$M0c_loglik
    }
    else {
      model <- stanmodels$M0c
    }
  }

  # M1
  if(model.name == "M1") {
    if(comparison.mode == TRUE) {
      model <- stanmodels$M1_loglik
    }
    else {
      model <- stanmodels$M1
    }
  }

  # M1c
  if(model.name == "M1c") {
    if(comparison.mode == TRUE) {
      model <- stanmodels$M1c_loglik
    }
    else {
      model <- stanmodels$M1c
    }
  }

  # M2
  if(model.name == "M2") {
    if(comparison.mode == TRUE) {
      model <- stanmodels$M2_loglik
    }
    else {
      model <- stanmodels$M2
    }
  }

  # M2c
  if(model.name == "M2c") {
    if(comparison.mode == TRUE) {
      model <- stanmodels$M2c_loglik
    }
    else {
      model <- stanmodels$M2c
    }
  }

  return (model)
}





# Function:
# Given model.name, fetch appropriate stan model
getStanModelDebug <- function(model.name, comparison = TRUE) {
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




getStanModelPars <- function(model.name, comparison = TRUE) {
  if(comparison) {
    # M0
    if(model.name == "M0") {
      pars <- c("alpha", "beta", "sigma",
                "mu_beta", "sigma_beta",
                "log_lik")
    }
    if(model.name == "M0c") {
      pars <- c("alpha", "beta", "sigma",
                "mu_beta", "sigma_beta",
                "log_lik", "rho")
    }

    # M1
    if(model.name == "M1") {
      pars <- c("alpha", "beta", "mu_beta", "grand_mu_beta", "sigma",
                "sigma_beta", "grand_sigma_beta", "log_lik")
    }
    if(model.name == "M1c") {
      pars <- c("alpha", "beta", "mu_beta", "grand_mu_beta", "sigma",
                "sigma_beta", "grand_sigma_beta", "rho", "log_lik")
    }

    # M2
    if(model.name == "M2") {
      pars <- c("alpha", "beta", "mu_beta", "grand_mu_beta",
                "sigma", "sigma_trait", "mean_trait", "sigma_beta",
                "grand_sigma_beta", "log_lik")
    }
    if(model.name == "M2c") {
      pars <- c("alpha", "beta", "mu_beta", "grand_mu_beta",
                "sigma", "sigma_trait", "mean_trait", "sigma_beta",
                "grand_sigma_beta", "log_lik", "rho")
    }
  }
  else {
    # M0
    if(model.name == "M0") {
      pars <- c("alpha", "beta", "sigma", "mu_beta", "sigma_beta")
    }
    if(model.name == "M0c") {
      pars <- c("alpha", "beta", "sigma", "mu_beta", "sigma_beta", "rho")
    }

    # M1
    if(model.name == "M1") {
      pars <- c("alpha", "beta", "mu_beta", "grand_mu_beta", "sigma",
                "sigma_beta", "grand_sigma_beta")
    }
    if(model.name == "M1c") {
      pars <- c("alpha", "beta", "mu_beta", "grand_mu_beta", "sigma",
                "sigma_beta", "grand_sigma_beta", "rho")
    }

    # M2
    if(model.name == "M2") {
      pars <- c("alpha", "beta", "mu_beta", "grand_mu_beta",
                "sigma", "sigma_trait", "mean_trait", "sigma_beta",
                "grand_sigma_beta")
    }
    if(model.name == "M2c") {
      pars <- c("alpha", "beta", "mu_beta", "grand_mu_beta",
                "sigma", "sigma_trait", "mean_trait", "sigma_beta",
                "grand_sigma_beta", "rho")
    }
  }

  return (pars)
}




# Function:
# Check the (...)  parameters for the model comparison function
checkDotParModelComparison <- function(...) {

  checkAdaptDelta <- function(adapt_delta) {
    if(length(adapt_delta) != 1) {
      stop("adapt_delta must be in range (0, 1) (default = 0.95).")
    }

    if(is.numeric(adapt_delta) == FALSE) {
      stop("adapt_delta must be in range (0, 1)")
    }

    if(adapt_delta >= 1 | adapt_delta <= 0) {
      stop("adapt_delta must be in range (0, 1)")
    }
  }

  checkMaxTreedepth <- function(max_treedepth) {
    if(length(max_treedepth) != 1) {
      stop("max_treedepth is numeric parameter.")
    }

    if(is.numeric(max_treedepth) == FALSE) {
      stop("max_treedepth is numeric parameter.")
    }

    if(max_treedepth < 5) {
      stop("max_treedepth >= 5 (default = 10).")
    }
  }

  checkVerbose <- function(verbose) {
    if(length(verbose) != 1) {
      stop("verbose is a logical parameter.")
    }

    if(is.logical(verbose) == FALSE) {
      stop("verbose is a logical parameter.")
    }
  }

  checkRefresh <- function(refresh) {
    if(length(refresh) != 1) {
      stop("refresh is numeric parameter.")
    }

    if(is.numeric(refresh) == FALSE) {
      stop("refresh is a numeric parameter.")
    }

    return (refresh)
  }

  checkSampleFile <- function(with.sample.file) {
    if(length(with.sample.file) != 1) {
      stop("with.sample.file is a logical parameter.")
    }

    if(is.logical(with.sample.file) == FALSE) {
      stop("with.sample.file is a logical parameter.")
    }
  }

  available.names <- c("adapt_delta",
                       "max_treedepth",
                       "refresh",
                       "verbose",
                       "with.sample.file")
  default.values <- list(adapt_delta = 0.95,
                         max_treedepth = 10,
                         refresh = 250,
                         verbose = TRUE,
                         with.sample.file = FALSE)

  # get the optional parameters
  dot.names <- names(list(...))

  if(length(dot.names) > 0) {
    if(any(dot.names %in% available.names) == FALSE) {
      wrong.names <- dot.names[!dot.names %in% available.names]
      stop("Unknown optional parameter were provided")
    }
  }

  # check each parameter
  for(p in dot.names) {
    if(is.null(list(...)[[p]]) || is.na(list(...)[[p]])) {
      stop(paste("optional parameter ", p, " can't be NULL", sep = ''))
    }
    if(p == "adapt_delta") {
      checkAdaptDelta(adapt_delta = list(...)[[p]])
      default.values[["adapt_delta"]] <- list(...)[[p]]
    }
    if(p == "max_treedepth") {
      checkMaxTreedepth(max_treedepth = list(...)[[p]])
      default.values[["max_treedepth"]] <- list(...)[[p]]
    }
    if(p == "refresh") {
      checkRefresh(refresh = list(...)[[p]])
      default.values[["refresh"]] <- list(...)[[p]]
    }
    if(p == "verbose") {
      checkVerbose(verbose = list(...)[[p]])
      default.values[["verbose"]] <- list(...)[[p]]
    }
    if(p == "sample.file") {
      checkSampleFile(verbose = list(...)[[p]])
      default.values[["sample.file"]] <- list(...)[[p]]
    }
  }

  return (default.values)
}




# Function:
# Check the (...) parameters for the main function
checkDotParRun <- function(...) {

  checkAdaptDelta <- function(adapt_delta) {
    if(length(adapt_delta) != 1) {
      stop("adapt_delta must be in range (0, 1) (default = 0.95).")
    }

    if(is.numeric(adapt_delta) == FALSE) {
      stop("adapt_delta must be in range (0, 1)")
    }

    if(adapt_delta >= 1 | adapt_delta <= 0) {
      stop("adapt_delta must be in range (0, 1)")
    }
  }

  checkMaxTreedepth <- function(max_treedepth) {
    if(length(max_treedepth) != 1) {
      stop("max_treedepth is numeric parameter.")
    }

    if(is.numeric(max_treedepth) == FALSE) {
      stop("max_treedepth is numeric parameter.")
    }

    if(max_treedepth < 5) {
      stop("max_treedepth >= 5 (default = 10).")
    }
  }

  checkCvFold <- function(cv.fold) {
    if(length(cv.fold) != 1) {
      stop("cv.fold must be in range (0, 1) (default = 0.66).")
    }

    if(is.numeric(cv.fold) == FALSE) {
      stop("cv.fold must be in range (0, 1)")
    }

    if(cv.fold >= 1 | cv.fold <= 0) {
      stop("cv.fold must be in range (0, 1)")
    }
  }

  checkNtree <- function(ntree) {
    if(length(ntree) != 1) {
      stop("ntree is numeric parameter.")
    }

    if(is.numeric(ntree) == FALSE) {
      stop("ntree is numeric parameter.")
    }

    if(ntree < 100) {
      stop("ntree >= 100 (default = 500).")
    }
  }

  checkVerbose <- function(verbose) {
    if(length(verbose) != 1) {
      stop("verbose is a logical parameter.")
    }

    if(is.logical(verbose) == FALSE) {
      stop("verbose is a logical parameter.")
    }
  }

  checkRefresh <- function(refresh) {
    if(length(refresh) != 1) {
      stop("refresh is numeric parameter.")
    }

    if(is.numeric(refresh) == FALSE) {
      stop("refresh is a numeric parameter.")
    }

    return (refresh)
  }

  checkSampleFile <- function(with.sample.file) {
    if(length(with.sample.file) != 1) {
      stop("with.sample.file is a logical parameter.")
    }

    if(is.logical(with.sample.file) == FALSE) {
      stop("with.sample.file is a logical parameter.")
    }
  }

  available.names <- c("adapt_delta",
                       "max_treedepth",
                       "ntree",
                       "cv.fold",
                       "refresh",
                       "verbose",
                       "with.sample.file")
  default.values <- list(adapt_delta = 0.95,
                         max_treedepth = 10,
                         ntree = 1000,
                         cv.fold = 0.66,
                         refresh = 250,
                         verbose = TRUE,
                         with.sample.file = FALSE)

  # get the optional parameters
  dot.names <- names(list(...))

  if(length(dot.names) > 0) {
    if(any(dot.names %in% available.names) == FALSE) {
      wrong.names <- dot.names[!dot.names %in% available.names]
      stop("Unknown optional parameter were provided")
    }
  }

  # check each parameter
  for(p in dot.names) {
    if(is.null(list(...)[[p]]) || is.na(list(...)[[p]])) {
      stop(paste("optional parameter ", p, " can't be NULL", sep = ''))
    }
    if(p == "adapt_delta") {
      checkAdaptDelta(adapt_delta = list(...)[[p]])
      default.values[["adapt_delta"]] <- list(...)[[p]]
    }
    if(p == "max_treedepth") {
      checkMaxTreedepth(max_treedepth = list(...)[[p]])
      default.values[["max_treedepth"]] <- list(...)[[p]]
    }
    if(p == "cv.fold") {
      checkCvFold(cv.fold = list(...)[[p]])
      default.values[["cv.fold"]] <- list(...)[[p]]
    }
    if(p == "ntree") {
      checkNtree(ntree = list(...)[[p]])
      default.values[["ntree"]] <- list(...)[[p]]
    }
    if(p == "refresh") {
      checkRefresh(refresh = list(...)[[p]])
      default.values[["refresh"]] <- list(...)[[p]]
    }
    if(p == "verbose") {
      checkVerbose(verbose = list(...)[[p]])
      default.values[["verbose"]] <- list(...)[[p]]
    }
    if(p == "sample.file") {
      checkSampleFile(verbose = list(...)[[p]])
      default.values[["sample.file"]] <- list(...)[[p]]
    }
  }

  return (default.values)
}





# Function:
# Check the input arguments of the main run, stop if violated.
checkInputMain <- function(genotype,
                           traits,
                           trait.type,
                           strains,
                           model,
                           mcmc.chains,
                           mcmc.steps,
                           mcmc.warmup,
                           cores,
                           hdi.level,
                           stat.learn.method,
                           cv.steps) {


  # Function:
  # Is a genotype (SNP) bi-allelic?
  checkBiallelic <- function(genotype) {

    isBi <- function(x) {
      return(length(unique(x)) <= 2)
    }

    is.bi <- apply(X = genotype, MARGIN = 2, FUN = isBi)

    return(all(is.bi == TRUE))
  }

  checkGenotypeTrait <- function(genotype, traits, trait.type) {
    # CHECK: genotype
    if(is.null(attr(genotype, "class")) == FALSE) {
      if(!attr(genotype, "class") %in% c("AAMultipleAlignment",
                                         "DNAMultipleAlignment")) {
        stop("genotype can be one of the following structures: vector (for a
             single SNP), matrix, data.frame or Biostrings structures such as
             DNAMultipleAlignment or AAMultipleAlignment")
      }
      else {
        temp <- as.matrix(genotype)
        if(nrow(temp) < 2 | ncol(temp) == 0) {
          stop("The genotypes cannot have less than two observations, or the
               number of genotypes cannot be 0.")
        }

        if(checkBiallelic(genotype = temp) == FALSE) {
          stop("Only bi-allelic genotypes allowed (each SNP
             must have at most 2 genotypes)")
        }
      }
    }
    else {
      if(is.vector(genotype)) {
        genotype <- matrix(data = genotype, ncol = 1)
      }

      if(nrow(genotype) < 2 | ncol(genotype) == 0) {
        stop("The genotypes cannot have less than two observations or the
             number of genotypes cannot be 0.")
      }

      if(!is.matrix(genotype) & !is.data.frame(genotype)) {
        stop("genotype can be one of the following structures: vector (for a
             single SNP), matrix, data.frame or Biostrings structures such as
             DNAMultipleAlignment or AAMultipleAlignment")
      }

      if(typeof(genotype) != "character") {
        stop("If it is structured as vector/matrix/data.frame,
             the genotype have be of character type.")
      }

      if(checkBiallelic(genotype = genotype) == FALSE) {
        stop("Only bi-allelic genotypes allowed (each SNP
             must have at most 2 genotypes)")
      }
    }



    # CHECK: traits
    if(!is.vector(traits) & !is.matrix(traits)) {
      stop("The traits must be either a vector (single traits) or matrix
           (with multiple traits = columns), where the rows match the rows
           of the genotype data")
    }

    # convert vector -> matrix
    if(is.vector(traits) == TRUE) {
      traits <- matrix(data = traits, ncol = 1)
    }

    if(!is.numeric(traits)) {
      stop("The traits must be of numeric type.")
    }

    if(nrow(traits) < 3) {
      stop("The traits must contain at least 3 data points.")
    }

    if(nrow(genotype) != nrow(traits)) {
      stop("length(genotype) != length(traits),
           they must be equal in length.")
    }



    # CHECK: trait.type
    if(is.vector(trait.type) == FALSE) {
      stop("trait.type must be vector. Each element in this vector refers
           to the type of each traits, with 'Q' (for quantitative traits)
           or 'D' (for dichotomous) included as each columns in the traits
           data. If a single trait is provided, then the trait.type should
           be a single 'Q' or 'D'.")
    }

    if(length(trait.type) == 0) {
      stop("The trait.type vector must contain at least 1 element.")
    }

    if(typeof(trait.type) != "character") {
      stop("trait.type must be character vector with elements 'Q' (for
           quantitative traits) or 'D' (for dichotomous)")
    }

    if(ncol(traits) != length(trait.type)) {
      stop("Number of traits provided differs from traits types.")
    }

    if(all(trait.type %in% c("Q", "D")) == FALSE) {
      stop("trait.type must be character vector with elements 'Q' (for
           quantitative traits) or 'D' (for dichotomous)")
    }
  }

  checkTraitValidity <- function(traits, trait.type) {

    # convert vector -> matrix
    if(is.vector(traits) == TRUE) {
      traits <- matrix(data = traits, ncol = 1)
    }

    # check traits types
    for(i in 1:length(trait.type)) {
      if(trait.type[i] == "D") {
        if(length(unique(traits[, i])) != 2) {
          stop("The dichotomous traits must contain exactly two
               categories (classes) \n")
        }
      }
      if(trait.type[i] == "Q") {
        if(length(unique(traits[, i])) <= 3) {
          warning("The quantitative trait/s contains 3 or less unique
                  elements \n")
        }
      }
    }
  }

  checkStrains <- function(strains, traits) {

    if(is.vector(strains) == FALSE) {
      stop("The strain identifiers must be provided in a vector")
    }

    if(is.character(strains) == FALSE &
       is.numeric(strains) == FALSE &
       is.factor(strains) == FALSE) {
      stop("The strain identifiers must be provided in a character or
           numeric or factor vector")
    }

    if(is.vector(traits) == TRUE) {
      if(length(traits) != length(strains)) {
        stop("One strain identifier per individual (trait measurement)
             must be provided")
      }
    }
    else {
      if(nrow(traits) != length(strains)) {
        stop("One strain identifier per individual (traits measurement)
             must be provided")
      }
    }

    if(length(strains) == 0) {
      stop("One strain identifier per individual (traits measurement)
             must be provided")
    }
  }

  checkModel <- function(model) {
    # CHECK: model
    if(length(model) != 1) {
      stop("model must be one of 'M0', 'M0c', 'M1', 'M1c', 'M2' or 'M2c'")
    }

    if(!is.character(model)) {
      stop("model must be of 'M0', 'M0c', 'M1', 'M1c', 'M2' or 'M2c'")
    }

    if(!model %in% c("M0", "M0c", "M1", "M1c", "M2", "M2c")) {
      stop("model must be of 'M0', 'M0c', 'M1', 'M1c', 'M2' or 'M2c'")
    }
  }

  checkMcmcIterations <- function(mcmc.steps) {
    # CHECK: mcmc.steps
    if(length(mcmc.steps) != 1) {
      stop("the mcmc.steps must be a number > 0 (default = 10000).")
    }

    if(!is.numeric(mcmc.steps)) {
      stop("mcmc.steps must be a numeric argument (default = 10000).")
    }

    if(mcmc.steps <= 0) {
      stop("mcmc.steps must be larger than 0 (default = 10000).")
    }
  }

  checkMcmcWarmup <- function(mcmc.warmup) {
    # CHECK: mcmc.warmup
    if(length(mcmc.warmup) != 1) {
      stop("the mcmc.warmup must be a number > 0 (default = 5000).")
    }

    if(!is.numeric(mcmc.warmup)) {
      stop("mcmc.warmup must be a numeric argument (default = 5000).")
    }

    if(mcmc.warmup <= 0) {
      stop("mcmc.warmup must be larger than 0 (default = 5000).")
    }
  }

  checkMcmcChains <- function(mcmc.chains) {
    # CHECK: mcmc.chains
    if(length(mcmc.chains) != 1) {
      stop("mcmc.chains must be a positive integer > 0 (default = 1).")
    }

    if(!is.numeric(mcmc.chains)) {
      stop("mcmc.chains must be a positive integer > 0 (default = 1).")
    }

    if(mcmc.chains <= 0) {
      stop("mcmc.chains must be a positive integer > 0 (default = 1).")
    }
  }

  checkCores <- function(cores) {
    # CHECK: cores
    if(length(cores) != 1) {
      stop("cores is numeric parameter.")
    }

    if(is.numeric(cores) == FALSE) {
      stop("cores is numeric parameter.")
    }

    if(cores <= 0) {
      stop("cores is numeric parameter >=1.")
    }
  }

  checkHdi <- function(hdi.level) {
    if(length(hdi.level) != 1) {
      stop("The HDI level must be in range (0, 1).")
    }

    if(is.numeric(hdi.level) == FALSE) {
      stop("The HDI level must be in range (0, 1).")
    }

    if(hdi.level >= 1 | hdi.level <= 0) {
      stop("The HDI level must be in range (0, 1).")
    }
  }

  checkMlMethod <- function(stat.learn.method) {
    # CHECK: trait.type
    if(length(stat.learn.method) != 1) {
      stop("stat.learn.method must be a string: 'rf' or 'svm'")
    }

    if(!is.character(stat.learn.method)) {
      stop("stat.learn.method must be a string: 'rf', 'svm' or 'none'")
    }

    if(!stat.learn.method %in% c("rf", "svm", "none")) {
      stop("stat.learn.method must be a string: 'rf', 'svm' or 'none'")
    }
  }

  checkCv <- function(stat.learn.method, cv.steps) {
    if(stat.learn.method %in% c("rf", "svm")) {
      if(length(cv.steps) != 1) {
        stop("cv.steps must be a number (default = 1,000).")
      }

      if(is.numeric(cv.steps) == FALSE) {
        stop("cv.steps must be a number (default = 1,000).")
      }

      if(cv.steps < 100) {
        stop("cv.steps >= 100 recomended (default = 1,000).")
      }
    }
  }

  if(is.null(genotype) | missing(genotype) |
     is.null(traits) | missing(traits) |
     is.null(trait.type) | missing(trait.type) |
     is.null(strains) | missing(strains) |
     is.null(model) | missing(model) |
     is.null(mcmc.chains) | missing(mcmc.chains) |
     is.null(mcmc.steps) | missing(mcmc.steps) |
     is.null(mcmc.warmup) | missing(mcmc.warmup) |
     is.null(cores) | missing(cores) |
     is.null(hdi.level) | missing(hdi.level) |
     is.null(stat.learn.method) | missing(stat.learn.method) |
     is.null(cv.steps) | missing(cv.steps)) {
    stop("arguments must be non-NULL/specified")
  }

  checkGenotypeTrait(genotype = genotype, traits = traits,
                     trait.type = trait.type)
  checkTraitValidity(traits = traits, trait.type = trait.type)
  checkStrains(strains = strains, traits = traits)
  checkModel(model = model)
  checkMcmcIterations(mcmc.steps = mcmc.steps)
  checkMcmcWarmup(mcmc.warmup = mcmc.warmup)
  checkMcmcChains(mcmc.chains = mcmc.chains)
  checkCores(cores = cores)
  checkHdi(hdi.level = hdi.level)
  checkMlMethod(stat.learn.method = stat.learn.method)
  checkCv(stat.learn.method = stat.learn.method, cv.steps = cv.steps)
}



# Function:
# Check the input arguments of the model comparison function
checkInputModelComparison <- function(genotype,
                                      traits,
                                      trait.type,
                                      strains,
                                      models,
                                      mcmc.chains,
                                      mcmc.steps,
                                      mcmc.warmup,
                                      cores,
                                      hdi.level) {


  # Function:
  # Is a genotype (SNP) bi-allelic?
  checkBiallelic <- function(genotype) {

    isBi <- function(x) {
      return(length(unique(x)) <= 2)
    }

    is.bi <- apply(X = genotype, MARGIN = 2, FUN = isBi)

    return(all(is.bi == TRUE))
  }

  checkGenotypeTrait <- function(genotype, traits, trait.type) {
    # CHECK: genotype
    if(is.null(attr(genotype, "class")) == FALSE) {
      if(!attr(genotype, "class") %in% c("AAMultipleAlignment",
                                         "DNAMultipleAlignment")) {
        stop("genotype can be one of the following structures: vector (for a
             single SNP), matrix, data.frame or Biostrings structures such as
             DNAMultipleAlignment or AAMultipleAlignment")
      }
      else {
        temp <- as.matrix(genotype)
        if(nrow(temp) < 2 | ncol(temp) == 0) {
          stop("The genotypes cannot have less than two observations, or the
               number of genotypes cannot be 0.")
        }

        if(checkBiallelic(genotype = temp) == FALSE) {
          stop("Only bi-allelic genotypes allowed (each SNP
               must have at most 2 genotypes)")
        }
      }
    }
    else {
      if(is.vector(genotype)) {
        genotype <- matrix(data = genotype, ncol = 1)
      }

      if(nrow(genotype) < 2 | ncol(genotype) == 0) {
        stop("The genotypes cannot have less than two observations or the
             number of genotypes cannot be 0.")
      }

      if(!is.matrix(genotype) & !is.data.frame(genotype)) {
        stop("genotype can be one of the following structures: vector (for a
             single SNP), matrix, data.frame or Biostrings structures such as
             DNAMultipleAlignment or AAMultipleAlignment")
      }

      if(typeof(genotype) != "character") {
        stop("If it is structured as vector/matrix/data.frame,
             the genotype have be of character type.")
      }

      if(checkBiallelic(genotype = genotype) == FALSE) {
        stop("Only bi-allelic genotypes allowed (each SNP
             must have at most 2 genotypes)")
      }
    }



    # CHECK: traits
    if(!is.vector(traits) & !is.matrix(traits)) {
      stop("The traits must be either a vector (single traits) or matrix
           (with multiple traits = columns), where the rows match the rows
           of the genotype data")
    }

    # convert vector -> matrix
    if(is.vector(traits) == TRUE) {
      traits <- matrix(data = traits, ncol = 1)
    }

    if(!is.numeric(traits)) {
      stop("The traits must be of numeric type.")
    }

    if(nrow(traits) < 3) {
      stop("The traits must contain at least 3 data points.")
    }

    if(nrow(genotype) != nrow(traits)) {
      stop("length(genotype) != length(traits),
           they must be equal in length.")
    }



    # CHECK: trait.type
    if(is.vector(trait.type) == FALSE) {
      stop("trait.type must be vector. Each element in this vector refers
           to the type of each traits, with 'Q' (for quantitative traits)
           or 'D' (for dichotomous) included as each columns in the traits
           data. If a single trait is provided, then the trait.type should
           be a single 'Q' or 'D'.")
    }

    if(length(trait.type) == 0) {
      stop("The trait.type vector must contain at least 1 element.")
    }

    if(typeof(trait.type) != "character") {
      stop("trait.type must be character vector with elements 'Q' (for
           quantitative traits) or 'D' (for dichotomous)")
    }

    if(ncol(traits) != length(trait.type)) {
      stop("Number of traits provided differs from traits types.")
    }

    if(all(trait.type %in% c("Q", "D")) == FALSE) {
      stop("trait.type must be character vector with elements 'Q' (for
           quantitative traits) or 'D' (for dichotomous)")
    }
  }

  checkTraitValidity <- function(traits, trait.type) {

    # convert vector -> matrix
    if(is.vector(traits) == TRUE) {
      traits <- matrix(data = traits, ncol = 1)
    }

    # check traits types
    for(i in 1:length(trait.type)) {
      if(trait.type[i] == "D") {
        if(length(unique(traits[, i])) != 2) {
          stop("The dichotomous traits must contain exactly two
               categories (classes) \n")
        }
      }
      if(trait.type[i] == "Q") {
        if(length(unique(traits[, i])) <= 3) {
          warning("The quantitative trait/s contains 3 or less unique
                  elements \n")
        }
      }
    }
  }

  checkStrains <- function(strains, traits) {

    if(is.vector(strains) == FALSE) {
      stop("The strain identifiers must be provided in a vector")
    }

    if(is.character(strains) == FALSE &
       is.numeric(strains) == FALSE &
       is.factor(strains) == FALSE) {
      stop("The strain identifiers must be provided in a character or
           numeric or factor vector")
    }

    if(is.vector(traits) == TRUE) {
      if(length(traits) != length(strains)) {
        stop("One strain identifier per individual (trait measurement)
             must be provided")
      }
    }
    else {
      if(nrow(traits) != length(strains)) {
        stop("One strain identifier per individual (traits measurement)
             must be provided")
      }
    }

    if(length(strains) == 0) {
      stop("One strain identifier per individual (traits measurement)
           must be provided")
    }
  }

  checkModels <- function(models) {
    # CHECK: model
    if(length(models) <= 0) {
      stop("models must contain a combination of the following models:
           'M0', 'M0c', 'M1', 'M1c', 'M2' or 'M2c'")
    }


    if(all(models %in% c("M0", "M0c", "M1", "M1c", "M2", "M2c")) == FALSE) {
      stop("models must contain a combination of the following models:
           'M0', 'M0c', 'M1', 'M1c', 'M2' or 'M2c'")
    }
  }

  checkMcmcIterations <- function(mcmc.steps) {
    # CHECK: mcmc.steps
    if(length(mcmc.steps) != 1) {
      stop("the mcmc.steps must be a number > 0 (default = 10000).")
    }

    if(!is.numeric(mcmc.steps)) {
      stop("mcmc.steps must be a numeric argument (default = 10000).")
    }

    if(mcmc.steps <= 0) {
      stop("mcmc.steps must be larger than 0 (default = 10000).")
    }
  }

  checkMcmcWarmup <- function(mcmc.warmup) {
    # CHECK: mcmc.warmup
    if(length(mcmc.warmup) != 1) {
      stop("the mcmc.warmup must be a number > 0 (default = 5000).")
    }

    if(!is.numeric(mcmc.warmup)) {
      stop("mcmc.warmup must be a numeric argument (default = 5000).")
    }

    if(mcmc.warmup <= 0) {
      stop("mcmc.warmup must be larger than 0 (default = 5000).")
    }
  }

  checkMcmcChains <- function(mcmc.chains) {
    # CHECK: mcmc.chains
    if(length(mcmc.chains) != 1) {
      stop("mcmc.chains must be a positive integer > 0 (default = 1).")
    }

    if(!is.numeric(mcmc.chains)) {
      stop("mcmc.chains must be a positive integer > 0 (default = 1).")
    }

    if(mcmc.chains <= 0) {
      stop("mcmc.chains must be a positive integer > 0 (default = 1).")
    }
  }

  checkCores <- function(cores) {
    # CHECK: cores
    if(length(cores) != 1) {
      stop("cores is numeric parameter.")
    }

    if(is.numeric(cores) == FALSE) {
      stop("cores is numeric parameter.")
    }

    if(cores <= 0) {
      stop("cores is numeric parameter >=1.")
    }
  }

  checkHdi <- function(hdi.level) {
    if(length(hdi.level) != 1) {
      stop("The HDI level must be in range (0, 1).")
    }

    if(is.numeric(hdi.level) == FALSE) {
      stop("The HDI level must be in range (0, 1).")
    }

    if(hdi.level >= 1 | hdi.level <= 0) {
      stop("The HDI level must be in range (0, 1).")
    }
  }

  if(is.null(genotype) | missing(genotype) |
     is.null(traits) | missing(traits) |
     is.null(trait.type) | missing(trait.type) |
     is.null(strains) | missing(strains) |
     is.null(models) | missing(models) |
     is.null(mcmc.chains) | missing(mcmc.chains) |
     is.null(mcmc.steps) | missing(mcmc.steps) |
     is.null(mcmc.warmup) | missing(mcmc.warmup) |
     is.null(cores) | missing(cores) |
     is.null(hdi.level) | missing(hdi.level)) {
    stop("arguments must be non-NULL/specified")
  }

  checkGenotypeTrait(genotype = genotype, traits = traits,
                     trait.type = trait.type)
  checkTraitValidity(traits = traits, trait.type = trait.type)
  checkStrains(strains = strains, traits = traits)
  checkModels(models = models)
  checkMcmcIterations(mcmc.steps = mcmc.steps)
  checkMcmcWarmup(mcmc.warmup = mcmc.warmup)
  checkMcmcChains(mcmc.chains = mcmc.chains)
  checkCores(cores = cores)
  checkHdi(hdi.level = hdi.level)
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
getScores <- function(p, s,
                      gt.data,
                      hdi.level) {
  probs <- c((1-hdi.level)/2, 1-(1-hdi.level)/2)
  d.full <- data.frame(summary(p$posterior, probs = probs)$summary,
                       stringsAsFactors = FALSE)
  d.full$par <- rownames(d.full)


  xmap <- gt.data$xmap
  scores <- vector(mode = "list", length = ncol(gt.data$Y))
  for(i in 1:ncol(gt.data$Y)) {
    s.p <- s[s$p == i, ]

    key <- paste("beta\\[", i, "\\,", sep = '')
    d <- d.full[which(regexpr(pattern = key, text = d.full$par) != -1
                      & regexpr(pattern = "alpha|sigma|nu|mu|tau|z",
                                text = d.full$par) == -1), ]
    d$i <- i
    d$i <- gsub(pattern = key, replacement = '', x = d$par)
    d$i <- gsub(pattern = "\\]", replacement = '', x = d$i)
    d$i <- as.numeric(d$i)
    scores[[i]] <- cbind(xmap, d, s[s$p == i, ])
  }

  return (scores)
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

