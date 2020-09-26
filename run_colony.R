library(ggplot2); library(snpR)

sample_meta <- read.table("metadata.txt", header = T)
genos <- read.table("bamlist2_baits_geno.geno", stringsAsFactors = F)
sample_meta$stream <- paste0("str_", sample_meta$stream)
sample_meta$stream[sample_meta$watershed == "LHU" & sample_meta$stream == "str_INC"] <- "str_INK"
snp_meta <- genos[,1:2]
colnames(snp_meta) <- c("chr", "position")


dat <- import.snpR.data(genos[,-c(1:2)], snp_meta, sample_meta)
pops <- unique(sample_meta[,3:4])
setwd("colony/")

for(i in i:nrow(pops)){
  tdat <- subset_snpR_data(dat, facets = "stream", subfacets =  pops[i,1])
  if(nrow(tdat@sample.meta) <= 1 | nrow(tdat@snp.meta) <= 10){
    out[[i]] <- data.frame(offspringID1 = character(), ofspringID2 = character(), Probability = numeric())
  }
  else{
    tdat <- subset_snpR_data(tdat, facets = "year", subfacets = pops[i,2])
    tdat <- filter_snps(tdat, maf = 0.01)
    if(nrow(tdat@sample.meta) <= 1 | nrow(tdat@snp.meta) <= 10){
      out[[i]] <- data.frame(OffspringID1 = character(), OffspringID2 = character(), Probability = numeric())
    }
    else{
      dir.create(paste0(pops[i,], collapse = "_"))
      setwd(paste0(pops[i,], collapse = "_"))
      write_colony_input(tdat, sampleIDs = "ID")
      call_colony("colony/colony_input.dat", colony_path = "D:/ZSL/Colony/colony2s.exe")
      out[[i]] <- read.csv("colony/colony_input.FullSibDyad", header = T, stringsAsFactors = F)
      setwd("../")
    }
  }
}
cout <- dplyr::bind_rows(out)

write.table(cout, "colony_sibs.txt", quote = F, sep = "\t", row.names = F)
