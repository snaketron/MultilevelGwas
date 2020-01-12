# Description:
# Parse input data and format it for stan and statistical learning
# Get genotype-trait data
getStanFormat <- function(genotype,
                          trait,
                          trait.type,
                          strains) {


  # convert AAMultipleAlignment to matrix if needed
  genotype <- convertMsaToGenotype(genotype = genotype)


  # if vector genotype => matrix genotype
  if(is.vector(genotype)) {
    genotype <- matrix(data = genotype, ncol = 1)
  }

  # TODO: check
  if(trait.type == "D") {
    if(all(trait %in% c(1, 0)) == FALSE) {
      # mapping 1st element to 1, 2nd to 0
      u <- unique(trait)
      trait[trait == u[1]] <- 1
      trait[trait == u[2]] <- 0
      cat("Mapping dichotomous trait:",
          "(", u[1], "->1,", u[2], "->0) \n")
    }
  }


  # X
  X <- apply(X = genotype, MARGIN = 2,
             FUN = function(x) {return(ifelse(
               test = x==x[1], yes = 1, no = -1))})


  d <- list(N = length(Y),
            Ns = ncol(X),
            Nk = length(unique(strains)),
            Y = trait,
            X = X,
            K = as.numeric(as.factor(strains)),
            strains = strains,
            P = ifelse(test = trait.type == "Q",
                       yes = 1, no = 0))

  return (d)
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
getStanModelDebugReal <- function(model.name) {

  # M0
  if(model.name == "M0") {
    stan.model <- rstan::stan_model(file = "src/stan_files/M0.stan")
  }

  # M1
  if(model.name == "M1") {
    stan.model <- rstan::stan_model(file = "src/stan_files/M1.stan")
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
    pars <- c("alpha_trait", "sigma", "beta_snp")
  }

  # M1
  if(model.name == "M1") {
    pars <- c("alpha_trait", "sigma",
              "beta_snp", "sigma_snp")
  }

  return (pars)
}



# Description:
# Check the input arguments of the main run, stop if violated.
checkInput <- function(genotype,
                       trait,
                       trait.type,
                       strains,
                       model,
                       mcmc.chains,
                       mcmc.steps,
                       mcmc.warmup,
                       mcmc.cores,
                       hdi.level,
                       adapt.delta,
                       max.treedepth,
                       write.samples) {


  if(is.null(genotype) | missing(genotype) |
     is.null(trait) | missing(trait) |
     is.null(trait.type) | missing(trait.type) |
     is.null(strains) | missing(strains) |
     is.null(model) | missing(model) |
     is.null(mcmc.chains) | missing(mcmc.chains) |
     is.null(mcmc.steps) | missing(mcmc.steps) |
     is.null(mcmc.warmup) | missing(mcmc.warmup) |
     is.null(mcmc.cores) | missing(mcmc.cores) |
     is.null(adapt.delta) | missing(adapt.delta) |
     is.null(max.treedepth) | missing(max.treedepth) |
     is.null(hdi.level) | missing(hdi.level) |
     is.null(write.samples) | missing(write.samples)) {
    stop("arguments must be non-NULL/specified")
  }

  checkGenotypeTrait(genotype = genotype, trait = trait,
                     trait.type = trait.type)
  checkTraitValidity(trait = trait, trait.type = trait.type)
  checkStrains(strains = strains, trait = trait)
  checkModel(model = model)
  checkMcmcSteps(mcmc.steps = mcmc.steps)
  checkMcmcWarmup(mcmc.warmup = mcmc.warmup)
  checkMcmcChains(mcmc.chains = mcmc.chains)
  checkMcmcCores(mcmc.cores = mcmc.cores)
  checkHdi(hdi.level = hdi.level)
  checkAdaptDelta(adapt.delta = adapt.delta)
  checkMaxTreedepth(max.treedepth = max.treedepth)
  checkWriteSamples(write.samples = write.samples)
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

  d <- summary(glm, probs = c(0.5, (1-hdi.level)/2,
                              1-(1-hdi.level)/2),
               pars = "beta_snp")$summary
  d <- data.frame(d)
  d$par <- rownames(d)
  d$par <- gsub(pattern = "beta_snp|\\[|\\]",
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



# Description:
# Computes max(P-, P+) for beta_snp
getPmax <- function(glm) {
  pmax <- c()

  e <- rstan::extract(object = glm, pars = c("beta_snp"))
  e.dim <- dim(e$beta_snp)
  Ns <- e.dim[3]
  Nt <- e.dim[2]
  Np <- e.dim[1]

  for(t in 1:Nt) {
    for(s in 1:Ns) {
      row <- data.frame(trait = t, site = s,
                        pmax = 2*max(sum(e$beta_snp[,t,s] > 0)/Np,
                                     sum(e$beta_snp[,t,s] < 0)/Np)-1,
                        stringsAsFactors = FALSE)
      pmax <- rbind(pmax, row)
    }
  }

  return (pmax)
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


