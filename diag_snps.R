library(snpR); library(data.table)


check <- data.table::fread("combined_diagnostics.bed")
check$V2 <- check$V2 - 1
check$V3 <- check$V3 - 1
data.table::fwrite(check, "combined_diagnostics.bed", col.names = FALSE, row.names = FALSE, sep = "\t")



#=====================prepare data and calculate basic statistics/Table S2=============
# prepare data
sample_meta <- read.table("metadata.txt", header = T)
genos <- read.table("bamlist2_combined_diagnostics.rf_geno.geno", stringsAsFactors = F)
sample_meta$stream <- paste0("str_", sample_meta$stream)
sample_meta$stream[sample_meta$watershed == "LHU" & sample_meta$stream == "str_INC"] <- "str_INK"
## merge pops ran together during pva
pva_merge_key <- data.frame(old = c("SEC", "TRB", "SHE", "SNO"), new = "GRP")
for(i in 1:nrow(pva_merge_key)){
  sample_meta$stream <- gsub(pva_merge_key[i,1], pva_merge_key[i,2], sample_meta$stream)
}
snp_meta <- genos[,1:2]
colnames(snp_meta) <- c("chr", "position")


# prepare snp meta
RBT_diag <- data.table::fread("RBT_All_diag_Filt4-19-18_loci_list.rf", sep = ":")
RBT_diag$sp <- "RBT"
YCT_diag <- data.table::fread("YCT_Diagnostics_3-30-18.rf", sep = ":")
YCT_diag$sp <- "YCT"
ref <- data.table::fread("combined_diagnostics_ref.txt", header = F)
ref$chr <- stringr::str_extract(ref$V1, "^.+:")
ref$chr <- gsub("\\:", "", ref$chr)
ref$position <- stringr::str_extract(ref$V1, ":[0-9]+")
ref$position <- gsub("\\:", "", ref$position)
ord <- c(3,4,2)
ref <- ref[,..ord]
colnames(ref)[3] <- "ref"
diag <- rbind(RBT_diag, YCT_diag)
colnames(diag)[1:2] <- c("chr", "position")
diag$position <- as.numeric(diag$position)
ref$position <- as.numeric(ref$position)
ref$position <- ref$position + 1
diag <- merge(ref, diag, all.x = TRUE, by = c("chr", "position"))
snp_meta <- merge(snp_meta, diag, all.x = TRUE, by = c("chr", "position"), sort = FALSE)



dat <- import.snpR.data(genos[,-c(1:2)], snp_meta, sample_meta)
dat <- calc_maf(dat, c(".base", "stream", "watershed"))

res <- get.snpR.stats(dat, "watershed", "maf")$single
reso <- get.snpR.stats(dat, ".base", "maf")$single
keep.col <- c(6, 3:5, 7:10)
reso <- dplyr::arrange(reso[,keep.col], sp, chr, position)
data.table::fwrite(reso, "diagnostic_snp_results_overall.txt", sep = "\t")

rbt <- reso[reso$sp == "RBT",]
sum(rbt$ref == rbt$minor)/nrow(rbt)

yct <- reso[reso$sp == "YCT",]
sum(yct$ref == yct$major)/nrow(yct)
