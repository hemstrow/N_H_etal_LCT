# add n for both the single pop stats and ns for both pops for fst.
# remove the nomaf stuff from the output.
# email new sheets, notes on filters and name changes, and scripts.

library(ggplot2); library(snpR)

sample_meta <- read.table("metadata.txt", header = T)
genos <- read.table("bamlist2_baits_geno.geno", stringsAsFactors = F)
sample_meta$stream <- paste0("str_", sample_meta$stream)
sample_meta$stream[sample_meta$watershed == "LHU" & sample_meta$stream == "str_INC"] <- "str_INK"
snp_meta <- genos[,1:2]
colnames(snp_meta) <- c("chr", "position")

# ne
dat <- import.snpR.data(genos[,-c(1:2)], snp_meta, sample_meta)
neS <- calc_ne(dat, "stream", chr = "chr", NeEstimator_path = "D://Users/hemst/Downloads/Ne2-1.exe", methods = c("LD", "het", "Coan"))
neW <- calc_ne(dat, "watershed", chr = "chr", NeEstimator_path = "D://Users/hemst/Downloads/Ne2-1.exe", methods = c("LD", "het", "Coan"))
neA <- calc_ne(dat, chr = "chr", NeEstimator_path = "D://Users/hemst/Downloads/Ne2-1.exe", methods = c("LD", "het", "Coan"))
neSY <- calc_ne(dat, "stream.year", chr = "chr", NeEstimator_path = "D://Users/hemst/Downloads/Ne2-1.exe", methods = c("LD", "het", "Coan"))
neWY <- calc_ne(dat, "watershed.year", chr = "chr", NeEstimator_path = "D://Users/hemst/Downloads/Ne2-1.exe", methods = c("LD", "het", "Coan"))
neS$ne$facet <- "stream"
neW$ne$facet <- "watershed"
neA$ne$facet <- ".base"
neSY$ne$facet <- "stream.year"
neWY$ne$facet <- "watershed.year"
ne <- list(neS$ne, neW$ne, neSY$ne, neWY$ne)
ne <- dplyr::bind_rows(ne)

# other stats
dat <- filter_snps(dat, 0.05, hf_hets = 0.55, maf.facets = "stream") # 44 removed from low maf, 8 from high heterozygosity, leaving 845/897
dat <- calc_basic_snp_stats(dat, c("watershed", "stream", "year", "watershed.year", "stream.year"))
dat <- calc_het_hom_ratio(dat, c("watershed", "stream", "year", "watershed.year", "stream.year"))
fst <- calc_pairwise_fst(dat, c("watershed", "stream", "year", "watershed.year", "stream.year"), method = "genepop")
dat <- calc_pairwise_fst(dat, c("watershed", "stream", "year", "watershed.year", "stream.year"))

hh <- dat@sample.stats
ss <- dat@stats
ho <- tapply(ss$ho, ss$subfacet, mean, na.rm = T)
pi <- tapply(ss$pi, ss$subfacet, mean, na.rm = T)
comb_stats <- data.frame(pop = names(ho), ho = ho, pi = pi)
ss$n <- ss$maj.count + ss$min.count
comb_stats$n <- tapply(ss$n, ss$subfacet, function(x) max(x)/2)
fst <- fst[[2]]
fst$p1 <- unlist(strsplit(fst$comparison, "~"))[c(T,F)]
fst$p2 <- unlist(strsplit(fst$comparison, "~"))[c(F,T)]
fst$n1 <- comb_stats$n[match(fst$p1, comb_stats$pop)]
fst$n2 <- comb_stats$n[match(fst$p2, comb_stats$pop)]
fst <- fst[,-c(3:4)]

ggplot(neW$ne, aes(x = pop, y = LDNe, color = pop, shape = as.factor(pcrit))) +
  theme_bw() + geom_point(position = position_dodge(width=0.5)) + 
  geom_errorbar(aes(ymin = LDNe_lCIj, ymax = LDNe_uCIj), position = position_dodge(width=0.5)) + 
  scale_y_continuous(limits  = c(-5, 250))

p_ne_SY <- neSY$ne
p_ne_SY$stream <- substr(p_ne_SY$pop, 5, 7)
p_ne_SY$year <- substr(p_ne_SY$pop, 9, 12)

ggplot(p_ne_SY, aes(y = LDNe, x = as.factor(pcrit))) +
  facet_grid(stream~year) +
  theme_bw() + geom_point(position = position_dodge(width=0.5)) + 
  geom_errorbar(aes(ymin = LDNe_lCIj, ymax = LDNe_uCIj), position = position_dodge(width=0.5)) + 
  scale_y_continuous(limits  = c(-5, 250))

p_ne_S <- neS$ne
p_ne_S$stream <- substr(p_ne_S$pop, 5, 7)
ggplot(p_ne_S[-which(is.infinite(p_ne_S$LDNe_lCIj)),], aes(x = stream, y = LDNe)) +
  facet_wrap(~pcrit, ncol = 1) +
  theme_bw() + geom_point(position = position_dodge(width=0.5)) + 
  geom_errorbar(aes(ymin = LDNe_lCIj, ymax = LDNe_uCIj), position = position_dodge(width=0.5)) + 
  scale_y_continuous(limits  = c(-5, 250))

# save
write.csv(hh[,c(1:4, 7)], "het_hom_ratios.csv", row.names = F, quote = F)
write.csv(comb_stats, "pi_ho.csv", row.names = F, quote = F)
write.csv(fst, "fst.csv", row.names = F, quote = F)
write.csv(ne, "ne.csv", row.names = F, quote = F)
hhp <- ggplot(hh, aes(x = watershed, y = `Het/Hom`, color = watershed)) + 
  geom_point() + theme_bw() + scale_color_viridis_d()
hhpy <- ggplot(hh, aes(x = as.factor(year), y = `Het/Hom`, color = as.factor(year))) + 
  geom_point() + theme_bw() + scale_color_viridis_d()
pc <- plot_clusters(dat, "watershed", c("PCA", "umap"))
ggsave("het_hom_watershed.pdf", hhp, width = 11, height = 8.5)
ggsave("het_hom_year.pdf", hhpy, width = 11, height = 8.5)
ggsave("pca_watershed.pdf", pc$plots$pca, width = 11, height = 8.5)
ggsave("umap_watershed.pdf", pc$plots$umap, width = 11, height = 8.5)

