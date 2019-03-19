# Function:
# Parse input data and format it for stan and statistical learning
getStanData <- function(genotype,
                        phenotype,
                        phenotype.type) {



  # Description:
  # Get genotype-phenotype data
  getFormattedGenphen <- function(genotype,
                                  phenotype,
                                  phenotype.type) {


    # convert AAMultipleAlignment to matrix if needed
    genotype <- convertMsaToGenotype(genotype = genotype)


    # if vector genotype => matrix genotype
    if(is.vector(genotype)) {
      genotype <- matrix(data = genotype, ncol = 1)
    }

    # if vector genotype => matrix genotype
    if(is.vector(phenotype)) {
      phenotype <- matrix(data = phenotype, ncol = 1)
    }

    phenotype <- data.frame(phenotype)


    # TODO: check
    if(sum(phenotype.type == "D") != 0) {
      d <- which(phenotype.type == "D")
      for(i in 1:length(d)) {
        if(all(phenotype[, d] %in% c(1, 0)) == FALSE) {
          # mapping 1st element to 1, 2nd to 0
          u <- unique(phenotype[, d])
          phenotype[phenotype[, d] == u[1], d] <- "1"
          phenotype[phenotype[, d] == u[2], d] <- "0"
          cat("Mapping dichotomous phenotype:", i,
              "(", u[1], "->1,", u[2], "->0) \n")
        }
        phenotype[, d] <- as.factor(as.character(phenotype[, d]))
      }
    }

    # return
    return (list(genotype = genotype,
                 phenotype = phenotype,
                 phenotype.type = phenotype.type))
  }


  # Description:
  # Convert genphen data to stan data
  getStanFormat <- function(f.data) {


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


    # make huge matrix with predictors
    X <- f.data$genotype
    colnames(X) <- 1:ncol(X)
    x.data <- apply(X = X, MARGIN = 2, FUN = getSplitX)
    for(i in 1:length(x.data)) {
      x.data[[i]]$x.map$site <- i
    }
    X <- do.call(cbind, lapply(X = x.data, FUN = getXData))
    x.map <- do.call(rbind, lapply(X = x.data, FUN = getXMap))


    # TODO: remove
    # filtering
    # filter.i <- which(x.map$refN > 1 & x.map$altN > 1)
    # X <- X[, filter.i]
    # if(length(X) == 0) {
    #   stop("All SNPs have 1 element (Problem)")
    # }
    # if(is.vector(X) == T) {
    #   X <- matrix(data = X, ncol = 1)
    # }
    # x.map <- x.map[filter.i, ]


    if(sum(phenotype.type == "Q") != 0) {
      Yq <- as.matrix(f.data$phenotype[, f.data$phenotype.type == "Q"])
    }
    else {
      Yq <- matrix(data = 0, ncol = 0, nrow = nrow(X))
    }

    if(sum(phenotype.type == "D") != 0) {
      Yd <- as.matrix(f.data$phenotype[, f.data$phenotype.type == "D"])
    }
    else {
      Yd <- matrix(data = 0, ncol = 0, nrow = nrow(X))
    }

    s <- list(X = X,
              Y = f.data$phenotype,
              Yq = Yq,
              Yd = Yd,
              N = nrow(X),
              Ns = ncol(f.data$genotype),
              Nsk = ncol(X),
              Ntq = ncol(Yq),
              Ntd = ncol(Yd),
              xmap = x.map,
              genotype = f.data$genotype,
              phenotype.type = f.data$phenotype.type)

    return (s)
  }


  # genphen data
  f.data <- getFormattedGenphen(genotype = genotype,
                                phenotype = phenotype,
                                phenotype.type = phenotype.type)


  # stan data
  stan.data <- getStanFormat(f.data = f.data)


  # return
  return (stan.data)
}





# Description:
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
# Is a genotype (SNP) bi-allelic?
isBiallelic <- function(genotype) {

  isBi <- function(x) {
    return(length(unique(x)) <= 2)
  }

  is.bi <- apply(X = genotype, MARGIN = 2, FUN = isBi)

  return(all(is.bi == TRUE))
}




# Function:
# Check the ... parameters for the model comparison function
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

  available.names <- c("adapt_delta",
                       "max_treedepth",
                       "refresh",
                       "verbose")
  default.values <- list(adapt_delta = 0.95,
                         max_treedepth = 10,
                         refresh = 250,
                         verbose = TRUE)

  # get the optional parameters
  dot.names <- names(list(...))

  if(length(dot.names) > 0) {
    if(any(dot.names %in% available.names) == FALSE) {
      wrong.names <- dot.names[!dot.names %in% available.names]
      stop(paste("Unknown optional parameter were provided! The following
                 optional parameters are available:", dot.names, sep = ' '))
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
  }

  return (default.values)
}




# Function:
# Check the ... parameters for the main function
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

  available.names <- c("adapt_delta",
                       "max_treedepth",
                       "ntree",
                       "cv.fold",
                       "refresh",
                       "verbose")
  default.values <- list(adapt_delta = 0.95,
                         max_treedepth = 10,
                         ntree = 1000,
                         cv.fold = 0.66,
                         refresh = 250,
                         verbose = TRUE)

  # get the optional parameters
  dot.names <- names(list(...))

  if(length(dot.names) > 0) {
    if(any(dot.names %in% available.names) == FALSE) {
      wrong.names <- dot.names[!dot.names %in% available.names]
      stop(paste("Unknown optional parameter were provided! The following
                 optional parameters are available:", dot.names, sep = ' '))
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
  }

  return (default.values)
  }




# Description:
# Provided the input arguments, this function checks their validity. It
# stops the execution if a problem is encountered and prints out warnings.
checkInput <- function(genotype,
                       phenotype,
                       phenotype.type,
                       model.type,
                       mcmc.chains,
                       mcmc.steps,
                       mcmc.warmup,
                       cores,
                       hdi.level,
                       stat.learn.method,
                       cv.steps,
                       rpa.iterations,
                       with.stan.obj) {

  checkGenotypePhenotype <- function(genotype, phenotype, phenotype.type) {
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
      }



    # CHECK: phenotype
    if(!is.vector(phenotype) & !is.matrix(phenotype)) {
      stop("The phenotype must be either a vector (single phenotype) or matrix
           (with multiple phenotypes = columns), where the rows match the rows
           of the genotype data")
    }

    # convert vector -> matrix
    if(is.vector(phenotype) == TRUE) {
      phenotype <- matrix(data = phenotype, ncol = 1)
    }

    if(!is.numeric(phenotype)) {
      stop("The phenotype must be of numeric type.")
    }

    if(length(phenotype) < 3) {
      stop("The phenotype must contain at least 3 data points.")
    }

    if(nrow(genotype) != nrow(phenotype)) {
      stop("length(genotype) != length(phenotype),
           they must be equal in length.")
    }



    # CHECK: phenotype.type
    if(is.vector(phenotype.type) == FALSE) {
      stop("phenotype.type must be vector. Each element in this vector refers
           to the type of each phenotype, with 'Q' (for quantitative phenotypes)
           or 'D' (for dichotomous) included as each columns in the phenotype
           data. If a single phenotype is provided, then the phenotype type
           should be a single 'Q' or 'D'.")
    }

    if(length(phenotype.type) == 0) {
      stop("The phenotype.type vector must contain at least 1 element.")
    }

    if(typeof(phenotype.type) != "character") {
      stop("phenotype.type must be character vector with elements 'Q' (for
           quantitative phenotypes) or 'D' (for dichotomous)")
    }

    if(ncol(phenotype) != length(phenotype.type)) {
      stop("Number of phenotypes provided differs from phenotypes types.")
    }

    if(all(phenotype.type %in% c("Q", "D")) == FALSE) {
      stop("phenotype.type must be character vector with elements 'Q' (for
           quantitative phenotypes) or 'D' (for dichotomous)")
    }
    }

  checkPhenotypeValidity <- function(phenotype, phenotype.type) {

    # convert vector -> matrix
    if(is.vector(phenotype) == TRUE) {
      phenotype <- matrix(data = phenotype, ncol = 1)
    }

    # check phenotype types
    for(i in 1:length(phenotype.type)) {
      if(phenotype.type[i] == "D") {
        if(length(unique(phenotype[, i])) != 2) {
          stop("The dichotomous phenotypes must contains exactly two
               categories (classes) \n")
        }
        }
      if(phenotype.type[i] == "Q") {
        if(length(unique(phenotype[, i])) <= 3) {
          warning("The quantitative phenotype/s contains 3 or less unique
                  elements \n")
        }
        }
      }
      }

  checkModelType <- function(model.type) {
    # CHECK: model.type
    if(length(model.type) != 1) {
      stop("model.type must be a string (default = 'univariate')")
    }

    if(!is.character(model.type)) {
      stop("model.type must be a string: 'univariate' or 'hierarchical'")
    }

    if(!model.type %in% c("univariate", "hierarchical")) {
      stop("phenotype.type must be a string: 'univariate' or 'hierarchical'")
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
    # CHECK: phenotype.type
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
     is.null(phenotype) | missing(phenotype) |
     is.null(phenotype.type) | missing(phenotype.type) |
     is.null(model.type) | missing(model.type) |
     is.null(mcmc.chains) | missing(mcmc.chains) |
     is.null(mcmc.steps) | missing(mcmc.steps) |
     is.null(mcmc.warmup) | missing(mcmc.warmup) |
     is.null(cores) | missing(cores) |
     is.null(hdi.level) | missing(hdi.level) |
     is.null(stat.learn.method) | missing(stat.learn.method) |
     is.null(cv.steps) | missing(cv.steps)) {
    stop("arguments must be non-NULL/specified")
  }

  checkGenotypePhenotype(genotype = genotype,
                         phenotype = phenotype,
                         phenotype.type = phenotype.type)
  checkPhenotypeValidity(phenotype = phenotype,
                         phenotype.type = phenotype.type)
  checkModelType(model.type = model.type)
  checkMcmcIterations(mcmc.steps = mcmc.steps)
  checkMcmcWarmup(mcmc.warmup = mcmc.warmup)
  checkMcmcChains(mcmc.chains = mcmc.chains)
  checkCores(cores = cores)
  checkHdi(hdi.level = hdi.level)
  checkMlMethod(stat.learn.method = stat.learn.method)
  checkCv(stat.learn.method = stat.learn.method, cv.steps = cv.steps)
    }




# Description:
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










# Description:
# Combine data from Bayesian inference and statistical learning
getScores <- function(p, s,
                      genphen.data,
                      hdi.level) {
  probs <- c((1-hdi.level)/2, 1-(1-hdi.level)/2)
  d.full <- data.frame(summary(p$posterior, probs = probs)$summary,
                       stringsAsFactors = FALSE)
  d.full$par <- rownames(d.full)


  xmap <- genphen.data$xmap
  scores <- vector(mode = "list", length = ncol(genphen.data$Y))
  for(i in 1:ncol(genphen.data$Y)) {
    s.p <- s[s$p == i, ]
    # list(results = results,
    #      ca.list = ca.list,
    #      kappa.list = kappa.list)

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

