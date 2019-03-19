
# Description:
# Combine data from Bayesian inference and statistical learning
# betas = complete posterior
# cas = list with matrix (posterior) elements for traits
# kappa = list with matrix (posterior) elements for traits
getParetoScores <- function(scores) {
  phenotype.ids <- unique(scores$phenotype.id)
  result <- c()
  for(p in phenotype.ids) {
    t <- scores[scores$phenotype.id == p, ]
    p <- rPref::high(abs(t$beta.mean))*rPref::high(t$kappa.mean)
    f <- rPref::psel(t, p, top = nrow(t))
    f$rank <- f[, ".level"]
    f[, ".level"] <- NULL
    result <- rbind(result, f)
  }
  # p <- rPref::high(abs(scores$beta.mean))*
  #   rPref::high(scores$ca.mean)*
  #   rPref::high(scores$kappa.mean)
  return (result)
}




