meta <- read.table("bamlist2")
meta[,1] <- gsub("new_samples\\/", "", meta[,1])
meta[,1] <- gsub("\\.bam", "", meta[,1])
meta$watershed <- substr(meta[,1], 1, 3)
meta$stream <- substr(meta[,1], 5, 7)
meta$year <- substr(meta[,1], 9, 12)
colnames(meta)[1] <- "ID"

write.table(meta, "metadata.txt", col.names = T, row.names = F, quote = F)
