library(ggplot2); library(snpR)

sample_meta <- read.table("metadata.txt", header = T)
genos <- read.table("bamlist2_baits_geno.geno", stringsAsFactors = F)
sample_meta$stream <- paste0("str_", sample_meta$stream)
snp_meta <- genos[,1:2]
colnames(snp_meta) <- c("chr", "position")


dat <- import.snpR.data(genos[,-c(1:2)], snp_meta, sample_meta)
dat <- calc_hwe(dat, "stream.year")
hwe <- get.snpR.stats(dat, "stream.year", "hwe")
hwe$single$low.p <- ifelse(hwe$single$pHWE <= 0.000001, 1, 0)
bad.loci <- tapply(hwe$single$low.p, hwe$single[,c("chr", "position")], sum, na.rm = T)
bad.loci <- reshape2::melt(bad.loci)
bad.loci <- na.omit(bad.loci)
bad.loci <- bad.loci[which(bad.loci$value > 0),] # four snps to remove
good.loci <- which(!do.call(paste0, snp.meta(dat)[,1:2]) %in% do.call(paste0, bad.loci[,1:2]))
dat <- import.snpR.data(genos[good.loci,-c(1:2)], snp_meta[good.loci,], sample_meta)


pops <- unique(sample_meta[,3:4])
dir.create("colony")
setwd("colony/")


#=========run and source sibs==============
out <- vector("list", nrow(pops))

for(i in 1:nrow(pops)){
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
      call_colony("colony/colony_input.dat", colony_path = "C:/usr/bin/Colony/colony2s.exe")
      out[[i]] <- read.csv("colony/colony_input.FullSibDyad", header = T, stringsAsFactors = F)
      setwd("../")
    }
  }
}
empties <- which(lapply(out, nrow) == 0)
cout <- dplyr::bind_rows(out[-empties])

write.table(cout, "colony_sibs.txt", quote = F, sep = "\t", row.names = F)




#===============harvest and save ne===================
setwd("..")
nef <- list.files(path ="colony/", pattern = "Ne", recursive = T)
out <- data.frame(pop = gsub("/.+", "", nef),
                  ne.rand = numeric(length(nef)), 
                  lcl.rand = numeric(length(nef)),
                  ucl.rand = numeric(length(nef)),
                  ne.nonrand = numeric(length(nef)),
                  lcl.nonrand = numeric(length(nef)),
                  ucl.nonrand = numeric(length(nef)))

for(i in 1:length(nef)){
  dat <- readLines(paste0("colony/", nef[i]))
  dat <- gsub(" ", "", dat)
  dat <- gsub("^.+=", "", dat)
  out[i,-1] <- as.numeric(dat[c(5,6,7,12,13,14)])
}

out[out == 2147483647] <- Inf
write.table(out, "colony_ne.txt", quote = F, sep = "\t", row.names = F)
