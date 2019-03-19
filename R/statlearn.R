

# Description:
# Given two vectors, one dependent (genotype) and one independent
# (phenotype), compute the classification accuracy of classifying
# the genotype from the phenotype alone (and corresponding HDI).
# The classification is computed using random forest.
runStatLearn <- function(genphen.data,
                         method,
                         cv.fold,
                         cv.steps,
                         ntree,
                         hdi.level,
                         cores) {


  # RF analysis
  runRf <- function(X, Y, cv.fold, cv.steps,
                    ntree, hdi.level, site) {


    getBoot <- function(D, cv.fold, cv.steps,
                        ntree, hdi.level) {

      # posterior output (one extra for multi-trait)
      posterior <- vector(mode = "list", length = ncol(D)-1)
      for(i in 1:length(posterior)) {
        posterior[[i]] <- matrix(data = NA, nrow = cv.steps, ncol = 2)
      }
      rm(i)


      D$Y <- as.character(D$Y)
      for(i in 1:cv.steps) {
        # sample at random
        s <- sample(x = 1:nrow(D), size = ceiling(x = cv.fold*nrow(D)),
                    replace = FALSE)
        train <- D[s, ]
        test <- D[-s, ]


        # TODO: dummy (check if inference needed at all e.g. 1 class only)
        if(length(unique(train$Y)) == 1) {
          for(j in 1:length(posterior)) {
            posterior[[j]][i, ] <- c(NA, NA)
          }
        }
        else {
          for(j in 1:length(posterior)) {

            # train classification model (try condition to avoid errors in
            # case only one-level predictor is train data)
            rf.out <- try(ranger::ranger(Y~., data = train[, c(1, j+1)],
                                         num.trees = ntree), silent = TRUE)

            if(class(rf.out) == "try-error") {
              posterior[[j]][i, 1:2] <- c(NA, NA)
            }
            else {
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

              posterior[[j]][i, 1:2] <- c(ca, k)
            }
          }
        }
      }

      # compute stats
      summary <- c()
      for(j in 1:length(posterior)) {
        ca <- posterior[[j]][, 1]
        ca <- ca[is.finite(ca)]
        ca.mean <- mean(x = ca)
        # get HDI
        ca.hdi <- getHdi(vec = posterior[[j]][,1], hdi.level = hdi.level)
        ca.L <- as.numeric(ca.hdi[1])
        ca.H <- as.numeric(ca.hdi[2])

        k <- posterior[[j]][, 2]
        k <- k[is.finite(k)]
        k.mean <- mean(x = k)
        # get HDI
        k.hdi <- getHdi(vec = posterior[[j]][, 2], hdi.level = hdi.level)
        k.L <- as.numeric(k.hdi[1])
        k.H <- as.numeric(k.hdi[2])

        # summary append
        row <- data.frame(ca = ca.mean, ca.L = ca.L, ca.H = ca.H,
                          k = k.mean, k.L = k.L, k.H = k.H, p = j)
        summary <- rbind(summary, row)
      }

      return (list(summary = summary,
                   posterior = posterior))
    }


    i <- which(X %in% c(1, -1))
    D <- data.frame(Y = as.character(X))
    colnames(Y) <- paste("P", 1:ncol(Y), sep = '')
    D <- cbind(D, Y)
    D <- D[i, ]


    # if c.v. can be done with given cv.fold and # of data points
    result <- c()
    if(ceiling(cv.fold*nrow(D)) == nrow(D)) {
      for(i in 1:ncol(Y)) {
        out <- data.frame(ca = NA, ca.L = NA, ca.H = NA,
                          k = NA, k.L = NA, k.H = NA,
                          p = i, stringsAsFactors = FALSE)
        result <- rbind(result, out)
      }
    }
    else {
      # run
      p <- getBoot(D = D,
                   cv.fold = cv.fold,
                   cv.steps = cv.steps,
                   hdi.level = hdi.level,
                   ntree = ntree)
      result <- rbind(result, p$summary)
    }
    return (list(result = result,
                 posterior = p$posterior))
  }


  # SVM analysis
  runSvm <- function(X, Y, cv.fold, cv.steps,
                     hdi.level, site) {


    getBoot <- function(D, cv.fold, cv.steps, hdi.level) {


      # posterior output (one extra for multi-trait)
      posterior <- vector(mode = "list", length = ncol(D)-1)
      for(i in 1:length(posterior)) {
        posterior[[i]] <- matrix(data = NA, nrow = cv.steps, ncol = 2)
      }
      rm(i)


      D$Y <- as.character(D$Y)
      for(i in 1:cv.steps) {
        # sample at random
        s <- sample(x = 1:nrow(D), size = ceiling(x = cv.fold*nrow(D)),
                    replace = FALSE)
        train <- D[s, ]
        test <- D[-s, ]


        # TODO: dummy (check if inference needed at all e.g. 1 class only)
        if(length(unique(train$Y)) == 1) {
          for(j in 1:length(posterior)) {
            posterior[[j]][i, ] <- c(NA, NA)
          }
        }
        else {
          for(j in 1:length(posterior)) {

            # train classification model (try condition to avoid errors in case
            # only one-level predictor is train data)
            svm.out <- try(e1071::svm(as.factor(Y)~.,data = train[, c(1, j+1)],
                                      type = "C-classification"), silent = TRUE)

            if(class(svm.out)[1] == "try-error") {
              posterior[[j]][i, 1:2] <- c(NA, NA)
            }
            else {
              # test classification model
              pr <- stats::predict(object = svm.out, newdata = test)


              # compute classification accuracy (1 - classification error)
              test$Y <- as.character(test$Y)
              pr <- as.character(pr)
              ca <- sum(test$Y == pr)/nrow(test)


              # compute k statistics
              k <- getKappa(real = test$Y, predicted = pr, aas = unique(D$Y))

              posterior[[j]][i, 1:2] <- c(ca, k)
            }
          }
        }
      }

      # compute stats
      summary <- c()
      for(j in 1:length(posterior)) {
        ca <- posterior[[j]][, 1]
        ca <- ca[is.finite(ca)]
        ca.mean <- mean(x = ca)
        # get HDI
        ca.hdi <- getHdi(vec = posterior[[j]][,1], hdi.level = hdi.level)
        ca.L <- as.numeric(ca.hdi[1])
        ca.H <- as.numeric(ca.hdi[2])

        k <- posterior[[j]][, 2]
        k <- k[is.finite(k)]
        k.mean <- mean(x = k)
        # get HDI
        k.hdi <- getHdi(vec = posterior[[j]][, 2], hdi.level = hdi.level)
        k.L <- as.numeric(k.hdi[1])
        k.H <- as.numeric(k.hdi[2])


        # summary append
        row <- data.frame(ca = ca.mean, ca.L = ca.L, ca.H = ca.H,
                          k = k.mean, k.L = k.L, k.H = k.H, p = j)
        summary <- rbind(summary, row)
      }

      return (list(summary = summary,
                   posterior = posterior))
    }


    i <- which(X %in% c(1, -1))
    D <- data.frame(Y = as.character(X))
    colnames(Y) <- paste("P", 1:ncol(Y), sep = '')
    D <- cbind(D, Y)
    D <- D[i, ]


    # if c.v. can be done with given cv.fold and # of data points
    result <- c()
    if(ceiling(cv.fold*nrow(D)) == nrow(D)) {
      for(i in 1:ncol(Y)) {
        out <- data.frame(ca = NA, ca.L = NA, ca.H = NA,
                          k = NA, k.L = NA, k.H = NA,
                          p = i, stringsAsFactors = FALSE)
        result <- rbind(result, out)
      }
    }
    else {
      # run
      p <- getBoot(D = D,
                   cv.fold = cv.fold,
                   cv.steps = cv.steps,
                   hdi.level = hdi.level)
      result <- rbind(result, p$summary)
    }

    return (list(result = result,
                 posterior = p$posterior))
  }


  # Format posterior
  getPosteriorFormat <- function(cas) {
    ca.list <- vector(mode = "list", length = length(cas[[1]]$posterior))
    kappa.list <- vector(mode = "list", length = length(cas[[1]]$posterior))

    # empty result matrix
    results <- c()
    for(i in 1:length(ca.list)) {
      # empty posterior matrices
      ca.m <- matrix(data = 0, nrow = nrow(cas[[1]]$posterior[[1]]),
                     ncol = length(cas))
      kappa.m <- matrix(data = 0, nrow = nrow(cas[[1]]$posterior[[1]]),
                        ncol = length(cas))


      for(j in 1:length(cas)) {

        # collect results
        if(i == 1) {
          row <- cas[[j]]$result
          row$s <- j
          results <- rbind(results, row)
        }


        # collect posteriors
        ca.m[, j] <- cas[[j]]$posterior[[i]][, 1]
        kappa.m[,j ] <- cas[[j]]$posterior[[i]][, 2]
      }

      ca.list[[i]] <- ca.m
      kappa.list[[i]] <- kappa.m
    }


    return (list(results = results,
                 ca.list = ca.list,
                 kappa.list = kappa.list))
  }


  # multicore classification
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  j <- NULL
  if(method == "rf") {
    cas <- (foreach(j = 1:ncol(genphen.data$X),
                    .export = c("getHdi", "getKappa"),
                    .packages = c("ranger")) %dopar%
              runRf(X = as.matrix(genphen.data$X[, j]),
                    Y = genphen.data$Y,
                    cv.fold = cv.fold,
                    cv.steps = cv.steps,
                    hdi.level = hdi.level,
                    ntree = ntree,
                    site = j))
  }
  else if(method == "svm") {
    cas <- (foreach(j = 1:ncol(genphen.data$X),
                    .export = c("getHdi", "getKappa"),
                    .packages = c("e1071")) %dopar%
              runSvm(X = as.matrix(genphen.data$X[, j]),
                     Y = genphen.data$Y,
                     cv.fold = cv.fold,
                     cv.steps = cv.steps,
                     hdi.level = hdi.level,
                     site = j))
  }
  # stop cluster
  parallel::stopCluster(cl = cl)
  doParallel::stopImplicitCluster()


  # format posterior
  cas <- getPosteriorFormat(cas = cas)

  return (cas)
}












