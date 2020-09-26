library(snpR)
# import genotypes
sample_meta <- read.table("metadata.txt", header = T)
genos <- read.table("bamlist2_baits_geno.geno", stringsAsFactors = F)
sample_meta$stream <- paste0("str_", sample_meta$stream)
sample_meta$stream[sample_meta$watershed == "LHU" & sample_meta$stream == "str_INC"] <- "str_INK"
snp_meta <- genos[,1:2]
colnames(snp_meta) <- c("chr", "position")
dat <- import.snpR.data(genos[,-c(1:2)], snp_meta, sample_meta)


# import cluster data from colony
flist <- list.files("colony/")
flist <- flist[-which(flist == "colony_sibs.txt")]
clusters <- vector("list", length(flist))
for(i in 1:length(flist)){
  clusters[[i]] <- readr::read_table2(paste0("colony/", flist[i], "/colony/colony_input.BestCluster"))
  clusters[[i]]$pop <- i
}
clusters <- dplyr::bind_rows(clusters)

# filter down to one individual per family with a probability of > .5, preferentially keeping the best genotyped individual
clusters <- as.data.frame(clusters[which(clusters$Probability >= .5),])

unique.clusters <- as.data.frame(unique(clusters[,c(1, 6)]))
genotype.quality <- colSums(dat != "NN")
genotype.quality <- data.frame(OffspringID = dat@sample.meta$ID, quality = genotype.quality)
clusters <- merge(clusters, genotype.quality, by = "OffspringID")

bads <- character()
for(i in 1:nrow(unique.clusters)){
  matches <- clusters[(clusters$pop == unique.clusters$pop[i] & clusters$ClusterIndex == unique.clusters$ClusterIndex[i]),]
  if(length(matches) > 1){
    bads <- c(bads, matches$OffspringID[-which.max(matches$quality)])
  }
}
dat.clean <- subset_snpR_data(dat, samps = which(!dat@sample.meta$ID %in% bads))


dat.clean <- calc_ne(dat.clean, facets = c("stream.year", "stream", "watershed", "watershed.year"))
neS <- calc_ne(dat.clean, "stream", chr = "chr", NeEstimator_path = "D://ne_estimator/Ne2-1.exe", methods = c("LD", "het", "Coan"))
neW <- calc_ne(dat.clean, "watershed", chr = "chr", NeEstimator_path = "D://ne_estimator/Ne2-1.exe", methods = c("LD", "het", "Coan"))
neSY <- calc_ne(dat.clean, "stream.year", chr = "chr", NeEstimator_path = "D://ne_estimator/Ne2-1.exe", methods = c("LD", "het", "Coan"))
neWY <- calc_ne(dat.clean, "watershed.year", chr = "chr", NeEstimator_path = "D://ne_estimator/Ne2-1.exe", methods = c("LD", "het", "Coan"))
neS$ne$facet <- "stream"
neW$ne$facet <- "watershed"
neSY$ne$facet <- "stream.year"
neWY$ne$facet <- "watershed.year"
ne <- list(neS$ne, neW$ne, neSY$ne, neWY$ne)
ne <- dplyr::bind_rows(ne)
write.csv(ne, "ne_sib_removed.csv", row.names = F, quote = F)

