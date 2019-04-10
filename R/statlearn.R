

# Description:
# Given two vectors, one dependent (genotype) and one independent
# (phenotype), compute the classification accuracy of classifying
# the genotype from the phenotype alone (and corresponding HDI).
# The classification is computed using random forest.
runStatLearn <- function(gt.data,
                         method,
                         cv.steps,
                         hdi.level,
                         cores,
                         dot.param) {


  # RF analysis
  runRf <- function(X, Y, cv.fold, cv.steps,
                    ntree, hdi.level, site,
                    trait.type) {


    getBoot <- function(D, cv.fold, cv.steps,
                        ntree, hdi.level) {

      # posterior output (one extra for multi-trait)
      ca.p <- matrix(data = NA, nrow = cv.steps, ncol = 1)
      k.p <- matrix(data = NA, nrow = cv.steps, ncol = 1)


      D$Y <- as.character(D$Y)
      for(i in 1:cv.steps) {
        # sample at random
        s <- sample(x = 1:nrow(D), size = ceiling(x = cv.fold*nrow(D)),
                    replace = FALSE)
        train <- D[s, ]
        test <- D[-s, ]


        # TODO: dummy (check if inference needed at all e.g. 1 class only)
        if(length(unique(train$Y)) != 1) {

          # train classification model (try condition to avoid errors in
          # case only one-level predictor is train data)
          rf.out <- try(ranger::ranger(Y~X, data = train,
                                       num.trees = ntree),
                        silent = TRUE)

          if(class(rf.out) != "try-error") {
            # test classification model
            pr <- stats::predict(object = rf.out, data = test)

            # compute classification accuracy (1 - classification error)
            test$Y <- as.character(test$Y)
            pr$predictions <- as.character(pr$predictions)
            ca <- sum(test$Y == pr$predictions)/nrow(test)


            # compute k statistics
            k <- getKappa(real = test$Y,
                          predicted = pr$predictions,
                          aas = unique(D$Y))

            ca.p[i, 1] <- ca
            k.p[i, 1] <- k
          }
        }
      }

      # compute stats
      summary <- c()
      ca <- ca.p
      ca <- ca[is.finite(ca)]
      ca.mean <- mean(x = ca)
      ca.median <- median(x = ca)
      # get HDI
      ca.hdi <- getHdi(vec = ca, hdi.level = hdi.level)
      ca.L <- as.numeric(ca.hdi[1])
      ca.H <- as.numeric(ca.hdi[2])
      ca.N <- length(ca)

      k <- k.p
      k <- k[is.finite(k)]
      k.mean <- mean(x = k)
      k.median <- median(x = k)
      # get HDI
      k.hdi <- getHdi(vec = k, hdi.level = hdi.level)
      k.L <- as.numeric(k.hdi[1])
      k.H <- as.numeric(k.hdi[2])
      k.N <- length(k)

      summary <- data.frame(ca.mean = ca.mean, ca.median = ca.median,
                            ca.L = ca.L, ca.H = ca.H, ca.N = ca.N,
                            k.mean = k.mean, k.median = k.median,
                            k.L = k.L, k.H = k.H, k.N = k.N)

      return (list(summary = summary, ca.p = ca.p, k.p = k.p))
    }


    X <- as.factor(X)
    if(trait.type == "D") {
      Y <- as.character(Y)
    }
    D <- data.frame(X = Y, Y = X)


    p <- getBoot(D = D,
                 cv.fold = cv.fold,
                 cv.steps = cv.steps,
                 hdi.level = hdi.level,
                 ntree = ntree)

    return (p)
  }


  # SVM analysis
  runSvm <- function(X, Y, cv.fold, cv.steps,
                     hdi.level, site,
                     trait.type) {


    getBoot <- function(D, cv.fold,
                        cv.steps,
                        hdi.level) {


      # posterior output (one extra for multi-trait)
      ca.p <- matrix(data = NA, nrow = cv.steps, ncol = 1)
      k.p <- matrix(data = NA, nrow = cv.steps, ncol = 1)


      D$Y <- as.character(D$Y)
      for(i in 1:cv.steps) {

        # sample at random
        s <- sample(x = 1:nrow(D), size = ceiling(x = cv.fold*nrow(D)),
                    replace = FALSE)
        train <- D[s, ]
        test <- D[-s, ]


        # TODO: dummy (check if inference needed at all e.g. 1 class only)
        if(length(unique(train$Y)) != 1) {

          # train classification model (try condition to avoid errors in case
          # only one-level predictor is train data)
          svm.out <- try(e1071::svm(as.factor(Y)~X,data = train,
                                    type = "C-classification"), silent = TRUE)

          if(class(svm.out)[1] != "try-error") {
            # test classification model
            pr <- stats::predict(object = svm.out, newdata = test)


            # compute classification accuracy (1 - classification error)
            test$Y <- as.character(test$Y)
            pr <- as.character(pr)
            ca <- sum(test$Y == pr)/nrow(test)


            # compute k statistics
            k <- getKappa(real = test$Y, predicted = pr, aas = unique(D$Y))


            # collect posteriors
            ca.p[i, 1] <- ca
            k.p[i, 1] <- k
          }
        }
      }

      # compute stats
      summary <- c()
      ca <- ca.p
      ca <- ca[is.finite(ca)]
      ca.mean <- mean(x = ca)
      ca.median <- median(x = ca)
      # get HDI
      ca.hdi <- getHdi(vec = ca, hdi.level = hdi.level)
      ca.L <- as.numeric(ca.hdi[1])
      ca.H <- as.numeric(ca.hdi[2])
      ca.N <- length(ca)

      k <- k.p
      k <- k[is.finite(k)]
      k.mean <- mean(x = k)
      k.median <- median(x = k)
      # get HDI
      k.hdi <- getHdi(vec = k, hdi.level = hdi.level)
      k.L <- as.numeric(k.hdi[1])
      k.H <- as.numeric(k.hdi[2])
      k.N <- length(k)

      summary <- data.frame(ca.mean = ca.mean, ca.median = ca.median,
                            ca.L = ca.L, ca.H = ca.H, ca.N = ca.N,
                            k.mean = k.mean, k.median = k.median,
                            k.L = k.L, k.H = k.H, k.N = k.N)

      return (list(summary = summary, ca.p = ca.p, k.p = k.p))
    }


    X <- as.factor(X)
    if(trait.type == "D") {
      Y <- as.character(Y)
    }
    D <- data.frame(X = Y, Y = X)


    p <- getBoot(D = D,
                 cv.fold = cv.fold,
                 cv.steps = cv.steps,
                 hdi.level = hdi.level)

    return (p)
  }



  # Format posterior:
  # out[[1]][[1]]
  #   * ca.p [x, y]
  #   * k.p [x, y]
  #   * summary [1, l]
  # (list(summary = summary, ca.p = ca.p, k.p = k.p))
  getPosteriorFormat <- function(out) {

    collectSummary <- function(x) {
      return(x$summary)
    }

    collectPosterior <- function(x, bit) {
      if(bit == "ca") {
        return(x$ca.p)
      }
      else {
        return(x$k.p)
      }
    }


    formatted.out <- vector(mode = "list", length = length(out))
    for(t in 1:length(out)) {
      out.summary <- do.call(rbind, lapply(X = out[[t]], FUN = collectSummary))
      out.summary$site <- 1:nrow(out.summary)
      out.summary$trait <- t

      out.ca <- do.call(cbind, lapply(X = out[[t]], FUN = collectPosterior, bit = "ca"))
      out.k <- do.call(cbind, lapply(X = out[[t]], FUN = collectPosterior, bit = "k"))

      formatted.out[[t]] <- list(summary = out.summary,
                                 p.ca = out.ca,
                                 p.ka = out.k)
    }

    return (formatted.out)
  }


  # multicore classification
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)


  # results list
  out <- vector(mode = "list", length = gt.data$Ntd+gt.data$Ntq)


  cat("Trait:")
  # loop across traits
  for(t in 1:(gt.data$Ntd+gt.data$Ntq)) {
    cat(t, "/", (gt.data$Ntd+gt.data$Ntq), ', ', sep = '')
    j <- NULL
    if(method == "rf") {
      out[[t]] <- (foreach(j = 1:ncol(gt.data$X),
                           .export = c("getHdi", "getKappa"),
                           .packages = c("ranger")) %dopar%
                     runRf(X = as.matrix(gt.data$X[, j]),
                           Y = gt.data$Y[, t],
                           cv.fold = dot.param$cv.fold,
                           cv.steps = cv.steps,
                           hdi.level = hdi.level,
                           ntree = dot.param$ntree,
                           site = j,
                           trait.type = gt.data$trait.type[t]))
    }
    else if(method == "svm") {
      out[[t]] <- (foreach(j = 1:ncol(gt.data$X),
                           .export = c("getHdi", "getKappa"),
                           .packages = c("e1071")) %dopar%
                     runSvm(X = as.matrix(gt.data$X[, j]),
                            Y = gt.data$Y[, t],
                            cv.fold = dot.param$cv.fold,
                            cv.steps = cv.steps,
                            hdi.level = hdi.level,
                            site = j,
                            trait.type = gt.data$trait.type[t]))
    }
  }
  cat("\n")

  # stop cluster
  doParallel::stopImplicitCluster()
  parallel::stopCluster(cl)

  # formatted
  formatted.out <- getPosteriorFormat(out = out)

  return (formatted.out)
}


