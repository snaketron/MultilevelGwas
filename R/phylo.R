

# Description:
# Given a genotype dataset containing SNPs (columns) and N individuals (rows),
# the procedure computes a NxN kinship matrix for the individuals and estimates
# the phylogenetic bias related to each SNP.
getPhyloBias <- function(genotype,
                         k.matrix) {
  phylo.bias <- c()

  # total mean phylogenetic distance
  mean.d.t <- mean(k.matrix[upper.tri(x = k.matrix, diag = FALSE)])

  for(i in 1:ncol(genotype)) {
    gs <- unique(genotype[, i])
    for(g in gs) {
      # feature mean phylogenetic distance
      mean.d.f <- mean(k.matrix[genotype[, i] == g, genotype[, i] == g])

      row <- data.frame(site = i, genotype = g, feature.dist = mean.d.f,
                        total.dist = mean.d.t, stringsAsFactors = FALSE)
      phylo.bias <- rbind(phylo.bias, row)
    }
  }

  return(phylo.bias)
}
