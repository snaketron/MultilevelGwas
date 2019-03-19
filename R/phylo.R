
# Function:
# Phylo bias analysis
runPhyloBiasCheck <- function(input.kinship.matrix,
                              genotype) {


  # check params
  checkInputPhyloBias(input.kinship.matrix = input.kinship.matrix,
                      genotype = genotype)


  # convert input data to phylo data
  phylo.data <- getPhyloData(genotype = genotype)


  # compute kinship if needed
  if(is.null(input.kinship.matrix) | missing(input.kinship.matrix)) {
    kinship.matrix <- e1071::hamming.distance(genotype)
  }
  else {
    kinship.matrix <- input.kinship.matrix
  }

  # compute bias
  bias <- getPhyloBias(genotype = genotype,
                       k.matrix = kinship.matrix)

  # bias = 1-dist(feature)/dist(total)
  bias$bias <- 1-bias$feature.dist/bias$total.dist


  # append bias to each SNP
  phylo.data$bias.ref <- NA
  phylo.data$bias.alt <- NA
  phylo.data$bias <- NA
  for(i in 1:nrow(phylo.data)) {
    bias.ref <- bias[bias$site == phylo.data$site[i] &
                       bias$genotype == phylo.data$ref[i], ]
    bias.alt <- bias[bias$site == phylo.data$site[i] &
                       bias$genotype == phylo.data$alt[i], ]
    phylo.data$bias.ref[i] <- bias.ref$bias[1]
    phylo.data$bias.alt[i] <- bias.alt$bias[1]
    phylo.data$bias[i] <- max(bias.ref$bias[1], bias.alt$bias[1])
  }


  # sort by site
  bias <- phylo.data[, c("site", "ref", "alt", "bias.ref", "bias.alt", "bias")]
  bias <- bias[order(bias$site, decreasing = FALSE), ]

  return (list(bias = bias, kinship.matrix = kinship.matrix))
}



# Function:
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
