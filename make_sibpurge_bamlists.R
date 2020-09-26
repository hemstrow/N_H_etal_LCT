library(snpR);
dat <- readRDS("sib_purged_snpRdata.RDS")
meta <- dat@sample.meta

pops <- unique(substr(meta$ID, 1, 12))
for(i in pops){
  matches <- meta$ID[grep(i, meta$ID)]
  matches <- paste0("/home/hemstrow/LCT/from_mike/LCT/new_samples/", matches, ".bam")
  write.table(as.data.frame(matches), paste0("sibpurge_bamlists/sibpurge_bamlist_", i), 
              row.names = F, col.names = F, quote = F)
}
