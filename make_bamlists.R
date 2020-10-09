library(snpR);

# sibpurge
dat <- readRDS("sib_purged_snpRdata.RDS")
meta <- dat@sample.meta

pops <- unique(meta[,3:4])
for(i in 1:nrow(pops)){
  matches <- meta$ID[which(meta$stream == pops$stream[i] & meta$year == pops$year[i])]
  matches <- paste0("/home/hemstrow/LCT/from_mike/LCT/new_samples/", matches, ".bam")
  write.table(as.data.frame(matches), paste0("sibpurge_bamlists/sibpurge_bamlist_", paste0(pops$stream[i], "_", pops$year[i])), 
              row.names = F, col.names = F, quote = F)
}



# no sibpurge
dat <- readRDS("no_sibpurge_snpRdata.RDS")
meta <- dat@sample.meta

pops <- unique(meta[,3:4])
for(i in 1:nrow(pops)){
  matches <- meta$ID[which(meta$stream == pops$stream[i] & meta$year == pops$year[i])]
  matches <- paste0("/home/hemstrow/LCT/from_mike/LCT/new_samples/", matches, ".bam")
  write.table(as.data.frame(matches), paste0("no_sibpurge_bamlists/sibpurge_bamlist_", paste0(pops$stream[i], "_", pops$year[i])), 
              row.names = F, col.names = F, quote = F)
}
