library(snpR); library(ggplot2)

dat <- readRDS("sib_purged_snpRdata.RDS")

#===============================calculate basic statistics==========================
dat <- filter_snps(dat, 0.05, HWE = 0.000001, maf.facets = "stream") # 44 removed from low maf, 8 from high heterozygosity, leaving 845/897
dat <- calc_basic_snp_stats(dat, c("watershed", "stream", "year", "watershed.year", "stream.year"))
dat <- calc_het_hom_ratio(dat, c("watershed", "stream", "year", "watershed.year", "stream.year"))
fst <- calc_pairwise_fst(dat, c("watershed", "stream", "year", "watershed.year", "stream.year"), method = "genepop")
dat <- fst[[1]]
fst <- fst[[2]]

# combine and reformat
hh <- dat@sample.stats
ss <- dat@stats
ho <- tapply(ss$ho, ss$subfacet, mean, na.rm = T)
pi <- tapply(ss$pi, ss$subfacet, mean, na.rm = T)
comb_stats <- data.frame(pop = names(ho), ho = ho, pi = pi)
ss$n <- ss$maj.count + ss$min.count
comb_stats$n <- tapply(ss$n, ss$subfacet, function(x) max(x)/2)
fst$p1 <- unlist(strsplit(fst$comparison, "~"))[c(T,F)]
fst$p2 <- unlist(strsplit(fst$comparison, "~"))[c(F,T)]
fst$n1 <- comb_stats$n[match(fst$p1, comb_stats$pop)]
fst$n2 <- comb_stats$n[match(fst$p2, comb_stats$pop)]
fst <- fst[,-c(3:4)]


#=============add in hybridization, ne===========
# get hybridization
hyb <- readr::read_csv("LCT_RBT_final hybridization 2-19-19.csv")
colnames(hyb)[1] <- "ID"
colnames(hyb)[5] <- "Relative.percent.RBT"
hyb$ID <- gsub("LHU_INC", "LHU_INK", hyb$ID)
hh$ID <- gsub("QUI_BAT_2013", "QUI_NFB_2013", hh$ID)
hh$stream <- gsub("str_BAT", "str_NFB", hh$stream)
comb <- merge(hh, hyb, by = "ID", all.x = T)

comb_stats$pop <- gsub("BAT", "NFB", comb_stats$pop)
comb$stream.year <- paste0(comb$stream, ".", comb$year)
heho_tab <- tapply(comb$'Het/Hom', comb$stream.year, mean)
heho_tab <- data.frame(locale = names(heho_tab), Het.Hom = as.numeric(heho_tab))
hyb_tab <- tapply(comb$Relative.percent.RBT, comb$stream.year, mean, na.omit = T)
hhyb_tab <- data.frame(locale = names(hyb_tab), Relative.percent.RBT = as.numeric(hyb_tab))
stats <- merge(heho_tab, comb_stats, by.x = "locale", by.y = "pop", all.x = T)
stats <- merge(stats, hyb_tab, by = "locale")

# get ne
ne <- read.csv("ne_sib_removed.csv")
ne$pop <- gsub("BAT", "NFB", ne$pop)
ne <- ne[ne$facet == "stream.year",]
ne <- ne[-which(is.infinite(ne$LDNe)),]
ne <- reshape2::dcast(ne[,c("pop", "pcrit", "LDNe")], pop ~ pcrit, value.var = "LDNe")
ne_colony <- read.table("colony_ne.txt", header = T)
ne_colony$pop <- paste0(substr(ne_colony$pop, 1, 7), ".", substr(ne_colony$pop, 9, 12))
ne_matched <- merge(ne_colony, ne, by = "pop")

stats$year <- gsub("^.+\\.", "", stats$locale)
stats <- merge(stats, ne_colony, by.x = "locale", by.y = "pop", all.x = T)




# save
write.csv(hh[,c(1:4, 7)], "het_hom_ratios.csv", row.names = F, quote = F)
write.csv(stats, "stream_year_combined_stats.csv", row.names = F, quote = F)
write.csv(fst, "fst.csv", row.names = F, quote = F)
write.csv(comb_stats, "pi_ho.csv")
# hhp <- ggplot(hh, aes(x = watershed, y = `Het/Hom`, color = watershed)) + 
#   geom_point() + theme_bw() + scale_color_viridis_d()
# hhpy <- ggplot(hh, aes(x = as.factor(year), y = `Het/Hom`, color = as.factor(year))) + 
#   geom_point() + theme_bw() + scale_color_viridis_d()
pc <- plot_clusters(dat, c("watershed.stream"), c("PCA", "umap"))
# ggsave("het_hom_watershed.pdf", hhp, width = 11, height = 8.5)
# ggsave("het_hom_year.pdf", hhpy, width = 11, height = 8.5)
# ggsave("pca_watershed.pdf", pc$plots$pca, width = 11, height = 8.5)
# ggsave("umap_watershed.pdf", pc$plots$umap, width = 11, height = 8.5)