require(MultilevelGwas)

# VSV
data.list <- get(load(file = "/media/simo/safehouse1/backup_nfs_simo/MiceGwas/output/posterior/vsv/data.list.RData"))
genotype <- data.list$X[, 1:300]
genotype[1, 1] <- as.character(genotype[1, 1])
traits <- cbind(data.list$Yq, data.list$Yd)
strains <- data.list$K


vsv.glm <- MultilevelGwas::runComparison(genotype = genotype,
                                         traits = traits,
                                         trait.type = c("Q", "Q"),
                                         strains = strains,
                                         models = c("M0", "M0c", "M1", "M1c"),
                                         mcmc.warmup = 500,
                                         mcmc.steps = 2500,
                                         mcmc.chains = 4,
                                         mcmc.cores = 4,
                                         hdi.level = 0.95,
                                         adapt.delta = 0.95,
                                         max.treedepth = 10)
