
# Description:
# Is the genotype (SNP) data bi-allelic?
checkBiallelic <- function(genotype) {

  # biallelic -> true, else -> false
  isBi <- function(x) {
    return(length(unique(x)) <= 2)
  }

  is.bi <- apply(X = genotype, MARGIN = 2, FUN = isBi)

  return(all(is.bi == TRUE))
}

# Description:
# Is the genotype, trait and trait.type data consistent?
checkGenotypeTrait <- function(genotype, trait, trait.type) {
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



  # CHECK: trait
  if(!is.vector(trait)) {
    stop("The trait must be provided as vector, with elements
         that match individuals")
  }

  if(!is.numeric(trait)) {
    stop("The trait must be of numeric type.")
  }

  if(length(trait) < 3) {
    stop("The trait vector must contain at least 3 data points.")
  }

  if(nrow(genotype) != length(trait)) {
    stop("nrow(genotype) != length(trait),
           they must be equal in length.")
  }



  # CHECK: trait.type
  if(is.character(trait.type) == FALSE) {
    stop("trait.type can be 'Q' or 'D'.")
  }

  if(length(trait.type) != 1) {
    stop("The trait.type vector must contain 1 element.")
  }

  if(typeof(trait.type) != "character") {
    stop("trait.type can be 'Q' or 'D'.")
  }

  if(!trait.type %in% c("Q", "D")) {
    stop("trait.type can be 'Q' or 'D'.")
  }
}

# Description:
# Is the trait valid?
checkTraitValidity <- function(trait, trait.type) {

  # check trait
  if(trait.type == "D") {
    if(length(unique(trait)) != 2) {
      stop("The dichotomous trait must contain exactly two
               categories (classes) \n")
    }
  }
  if(trait.type == "Q") {
    if(length(unique(trait)) <= 3) {
      stop("The quantitative trait must contain >= 3 unique data points")
    }
  }
}

# Description:
# Are the strains valid?
checkStrains <- function(strains, trait) {

  if(is.vector(strains) == FALSE) {
    stop("The strain identifiers must be provided in a vector")
  }

  if(is.character(strains) == FALSE &
     is.numeric(strains) == FALSE &
     is.factor(strains) == FALSE) {
    stop("The strain identifiers must be provided in a character or
           numeric or factor vector")
  }

  if(is.vector(trait) == TRUE) {
    if(length(trait) != length(strain)) {
      stop("One strain identifier per individual (trait measurement)
             must be provided")
    }
  }

  if(length(strains) == 0) {
    stop("One strain identifier per individual (trait measurement)
         must be provided")
  }
}

# Description:
# Is the specified model valid?
checkModel <- function(model) {
  # CHECK: model
  if(length(model) != 1) {
    stop("model must be one of 'M0', 'M0c', 'M1', 'M1c'")
  }

  if(!is.character(model)) {
    stop("model must be of 'M0', 'M0c', 'M1', 'M1c'")
  }

  if(!model %in% c("M0", "M0c", "M1", "M1c")) {
    stop("model must be of 'M0', 'M0c', 'M1', 'M1c'")
  }
}

# Description:
# Is the specified mcmc.step valid?
checkMcmcSteps <- function(mcmc.steps) {
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

# Description:
# Is the specified mcmc.warmup valid?
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

# Description:
# Is the specified mcmc.chain valid?
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

# Description:
# Is the specified mcmc.cores valid?
checkMcmcCores <- function(mcmc.cores) {
  # CHECK: mcmc.cores
  if(length(mcmc.cores) != 1) {
    stop("mcmc.cores is numeric parameter.")
  }

  if(is.numeric(mcmc.cores) == FALSE) {
    stop("mcmc.cores is numeric parameter.")
  }

  if(mcmc.cores <= 0) {
    stop("mcmc.cores is numeric parameter >=1.")
  }
}

# Description:
# Is the specified hdi.level valid?
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

# Description:
# Is the specified adapt.delta valid?
checkAdaptDelta <- function(adapt.delta) {
  if(length(adapt.delta) != 1) {
    stop("adapt.delta must be in range (0, 1) (default = 0.95).")
  }

  if(is.numeric(adapt.delta) == FALSE) {
    stop("adapt.delta must be in range (0, 1)")
  }

  if(adapt.delta >= 1 | adapt.delta <= 0) {
    stop("adapt.delta must be in range (0, 1)")
  }
}

# Description:
# Is the specified max.treedepth valid?
checkMaxTreedepth <- function(max.treedepth) {
  if(length(max.treedepth) != 1) {
    stop("max.treedepth is numeric parameter.")
  }

  if(is.numeric(max.treedepth) == FALSE) {
    stop("max.treedepth is numeric parameter.")
  }

  if(max.treedepth < 5) {
    stop("max.treedepth >= 5 (default = 10).")
  }
}

# Description:
# Is the specified write.samples valid
checkWriteSamples <- function(write.samples) {
  if(length(write.samples) != 1) {
    stop("write.samples is a logical parameter.")
  }

  if(is.logical(write.samples) == FALSE) {
    stop("write.samples is a logical parameter.")
  }
}
