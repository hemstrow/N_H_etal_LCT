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

bd <- dat@sample.meta[which(dat@sample.meta$ID %in% bads),]
bdt <- table(paste(bd$watershed, bd$stream, bd$year, sep = "_"))
names(bdt) <- gsub("_str_", "_", names(bdt))
bdt <- as.data.frame(bdt)
colnames(bdt) <- c("Sample", "Number_of_Removed_Siblings")
write.table(bdt, "Sibling_removal_counts.txt", sep = "\t", quote = F, col.names = T, row.names = F)

dat.clean <- subset_snpR_data(dat, samps = which(!dat@sample.meta$ID %in% bads))
saveRDS(dat.clean, "sib_purged_snpRdata.RDS")
