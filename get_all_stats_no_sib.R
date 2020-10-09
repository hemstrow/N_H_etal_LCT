library(dplyr); library(ranger); library(lmerTest); library(lme4); library(ggplot2); library(openxlsx)
library(sjPlot); library(sjmisc); library(sjlabelled); library(plotly); library(snpR)
source("make_pred_and_make_AIC_comp.R")
wb <- openxlsx::createWorkbook(creator = "William Hemstrom")
openxlsx::addWorksheet(wb, "genetic_stats")
openxlsx::addWorksheet(wb, "AIC_tab")
openxlsx::addWorksheet(wb, "Fst")
openxlsx::addWorksheet(wb, "Tajima's D plot")
openxlsx::addWorksheet(wb, "Hyb Correction plot")
openxlsx::addWorksheet(wb, "rf Extinction plot")
openxlsx::addWorksheet(wb, "rf Abundance plot")






#=====================prepare data and calculate basic statistics=============
# prepare data
sample_meta <- read.table("metadata.txt", header = T)
genos <- read.table("bamlist2_baits_geno.geno", stringsAsFactors = F)
sample_meta$stream <- paste0("str_", sample_meta$stream)
sample_meta$stream[sample_meta$watershed == "LHU" & sample_meta$stream == "str_INC"] <- "str_INK"
snp_meta <- genos[,1:2]
colnames(snp_meta) <- c("chr", "position")
dat <- import.snpR.data(genos[,-c(1:2)], snp_meta, sample_meta)

# calculate basic stats
neSY <- calc_ne(dat, "stream.year", chr = "chr", NeEstimator_path = "D://ne_estimator/Ne2-1.exe", methods = c("LD"))
dat <- filter_snps(dat, HWE = 0.000001)
dat <- calc_basic_snp_stats(dat, c("stream.year"))
dat <- calc_het_hom_ratio(dat, c("stream.year"))
fst <- calc_pairwise_fst(dat, c("stream.year"), method = "genepop")
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

fstc <- reshape2::dcast(fst, p1 ~ p2, value.var = "overall_fst")
colnames(fstc)[1] <- "Overall Fst"
fstn1 <- reshape2::dcast(fst, p1 ~ p2, value.var = "n1")
colnames(fstn1)[1] <- "n pop 1"
fstn2 <- reshape2::dcast(fst, p1 ~ p2, value.var = "n2")
colnames(fstn2)[1] <- "n pop 2"
fst$nt <- fst$n1 + fst$n2
fstnt <- reshape2::dcast(fst, p1 ~ p2, value.var = "nt")
colnames(fstnt)[1] <- "n total"


openxlsx::writeData(wb, "Fst", fstc, keepNA = T)
openxlsx::writeData(wb, "Fst", fstn1, startRow = nrow(fstc) + 4, keepNA = T); at.row <- nrow(fstc) + 4
openxlsx::writeData(wb, "Fst", fstn2, startRow = at.row + nrow(fstc) + 4, keepNA = T)
openxlsx::writeData(wb, "Fst", fstnt, startRow = at.row + nrow(fstc) + 4, keepNA = T)

#=================add in hybridization info============
# get hybridization info
hyb <- readr::read_csv("LCT_RBT_final hybridization 2-19-19.csv")
colnames(hyb)[1] <- "ID"
colnames(hyb)[5] <- "Relative.percent.RBT"
# hh$ID <- gsub("QUI_BAT_2013", "QUI_NFB_2013", hh$ID)
# hh$stream <- gsub("str_BAT", "str_NFB", hh$stream)
comb <- merge(hh, hyb, by = "ID", all.x = T) # note, this will auto-subset by whatever samples are still in they analysis! (removing poorly sequenced individuals)

#================summerize by population==============
comb_stats$pop <- gsub("BAT", "NFB", comb_stats$pop)
comb$stream.year <- paste0(comb$stream, ".", comb$year)
heho_tab <- tapply(comb$'Het/Hom', comb$stream.year, mean)
heho_tab <- data.frame(locale = names(heho_tab), Het.Hom = as.numeric(heho_tab))
hyb_tab <- tapply(comb$Relative.percent.RBT, comb$stream.year, mean, na.omit = T)
hhyb_tab <- data.frame(locale = names(hyb_tab), Relative.percent.RBT = as.numeric(hyb_tab))
stats <- merge(heho_tab, comb_stats, by.x = "locale", by.y = "pop", all.x = T)
stats <- merge(stats, hhyb_tab, by = "locale")


#===============merge in LDNe info======================
colnames(neSY$ne)[1] <- "locale"
LDNe <- neSY$ne
LDNe <- LDNe[,c(1,2,3)]
LDNe$LDNe[is.infinite(LDNe$LDNe)] <- NA
LDNec <- reshape2::dcast(LDNe, locale ~ pcrit, value.var = "LDNe")
#LDNec$locale <- gsub("BAT", "NFB", LDNec$locale)
colnames(LDNec) <- c("locale", paste0("LDNe_", colnames(LDNec)[-1]))
stats <- merge(stats, LDNec, "locale", all.x = T)

#===============merge in colony Ne info================
ne_colony <- read.table("colony_ne.txt", header = T)
ne_colony$pop <- paste0(substr(ne_colony$pop, 1, 7), ".", substr(ne_colony$pop, 9, 12))
colnames(ne_colony)[1] <- "locale"
stats <- merge(stats, ne_colony, by = "locale", all.x = T)

#===============merge in PVA data======================
PVA <- readr::read_csv("Table_Thetas_hyb_PVA_03062020.csv")
PVA$`PVA Extinction` <- as.numeric(gsub("%", "", PVA$`PVA Extinction`))
PVA$`PVA Ext U95` <- as.numeric(gsub("%", "", PVA$`PVA Ext U95`))
PVA$`PVA Ext L95` <- as.numeric(gsub("%", "", PVA$`PVA Ext L95`))
PVA$locale <- paste0(PVA$Stream_code, ".", PVA$`Genetic Year`)

stats <- merge(stats, PVA, by = "locale", all = T)

#===============fix column names=======================
colnames(stats) <- gsub(" ", ".", colnames(stats))
colnames(stats) <- gsub("/", ".", colnames(stats))
colnames(stats) <- gsub("%", "percent", colnames(stats))
colnames(stats) <- gsub("’", "", colnames(stats))
colnames(stats) <- gsub(")", "", colnames(stats))
colnames(stats) <- gsub("\\(", "", colnames(stats))
colnames(stats) <- gsub("Genetic.Year", "year", colnames(stats))

#==============model AIC comparison============
hh.AIC_comp <- make_AIC_comp(stats, "Het.Hom")

ho.AIC_comp <- make_AIC_comp(stats, "ho")

pi.AIC_comp <- make_AIC_comp(stats, "pi")

# note: no effect
ne.AIC_comp <- make_AIC_comp(stats, "ne.rand")

ldne_.01.AIC_comp <- make_AIC_comp(stats, "LDNe_0.01")

ldne_.02.AIC_comp <- make_AIC_comp(stats, "LDNe_0.02")

thetaT.AIC_comp <- make_AIC_comp(stats, "Tajimas.theta.4Nu")

thetaW.AIC_comp <- make_AIC_comp(stats, "Wattersons.theta.4Nu")

AIC_tab <- dplyr::bind_rows(hh.AIC_comp$AIC, ho.AIC_comp$AIC, pi.AIC_comp$AIC, ne.AIC_comp$AIC,
                            ldne_.01.AIC_comp$AIC, ldne_.02.AIC_comp$AIC, thetaT.AIC_comp$AIC,
                            thetaW.AIC_comp$AIC)

#write.table(AIC_tab, "AIC_table_no_sibling_removal.txt", col.names = T, row.names = F, quote = F, sep = "\t")
openxlsx::writeData(wb, "AIC_tab", x = AIC_tab, keepNA = T)

#==============get hyb corrected values===============
hh.pred <- make_pred(hh.AIC_comp)
ggplot(hh.pred$pre, aes(x = Relative.percent.RBT, y = Het.Hom - pred_Het.Hom)) + geom_point() + theme_bw()

ho.pred <- make_pred(ho.AIC_comp)
ggplot(ho.pred$pred, aes(x = Relative.percent.RBT, y = ho - pred_ho)) + geom_point() + theme_bw()

pi.pred <- make_pred(pi.AIC_comp)
ggplot(pi.pred$pred, aes(x = Relative.percent.RBT, y = pi - pred_pi)) + geom_point() + theme_bw()

# note: no effect
ne.pred <- make_pred(ne.AIC_comp)
ggplot(ne.pred$pred, aes(x = Relative.percent.RBT, y = ne.rand - pred_ne.rand)) + geom_point() + theme_bw()

ldne_.01.pred <- make_pred(ldne_.01.AIC_comp)
ggplot(ldne_.01.pred$pred, aes(x = Relative.percent.RBT, y = LDNe_0.01 - pred_LDNe_0.01)) + geom_point() + theme_bw()

ldne_.02.pred <- make_pred(ldne_.02.AIC_comp)
ggplot(ldne_.02.pred$pred, aes(x = Relative.percent.RBT, y = LDNe_0.02 - pred_LDNe_0.02)) + geom_point() + theme_bw()

# ldne_.05.pred <- make_pred(stats, "LDNe_0.05")
# ggplot(ldne_.05.pred$pred, aes(x = Relative.percent.RBT, y = LDNe_0.05 - pred_LDNe_0.05)) + geom_point() + theme_bw()


thetaT.pred <- make_pred(thetaT.AIC_comp)
ggplot(thetaT.pred$pred, aes(x = Relative.percent.RBT, y = Tajimas.theta.4Nu - pred_Tajimas.theta.4Nu)) + geom_point() + theme_bw()

thetaW.pred <- make_pred(thetaW.AIC_comp)
ggplot(thetaW.pred$pred, aes(x = Relative.percent.RBT, y = Wattersons.theta.4Nu - pred_Wattersons.theta.4Nu)) + geom_point() + theme_bw()

#=================create a table of model results============================
mlist <- list(hh.pred$mod, ho.pred$mod, pi.pred$mod, ne.pred$mod, ldne_.01.pred$mod, ldne_.02.pred$mod,
              thetaT.pred$mod, thetaW.pred$mod)
tab_model(mlist,
          show.aic = T, 
          digits.re = 7,
          digits = 7,
          collapse.ci = T, 
          pred.labels = c("(Intercept)", "%RBT", "year"), 
          show.icc = F, file = "Model_summary_no_sib_removal.html", CSS = css_theme("cells"))



#=================create exportable table of predicted/obs stats=============
merge.list <- c("Creek", "Basin", "year", "Relative.percent.RBT")
stats_table <- merge(pi.pred$pred, ho.pred$pred, by = merge.list, all = T)
stats_table <- merge(stats_table, hh.pred$pred, by = merge.list, all = T)             
stats_table <- merge(stats_table, thetaT.pred$pred, by = merge.list, all = T)             
stats_table <- merge(stats_table, thetaW.pred$pred, by = merge.list, all = T)
stats_table$pred_theta_diff <- (stats_table$pred_Tajimas.theta.4Nu - stats_table$pred_Wattersons.theta.4Nu)/rowMeans(stats_table[,c(11,13)])
stats_table$theta_diff <- (stats_table$Tajimas.theta.4Nu - stats_table$Wattersons.theta.4Nu)/rowMeans(stats_table[,c(12,14)])
stats_table <- merge(stats_table, ne.pred$pred, by = merge.list, all = T)
stats_table <- merge(stats_table, ldne_.01.pred$pred, by = merge.list, all = T)
stats_table <- merge(stats_table, ldne_.02.pred$pred, by = merge.list, all = T)
stats_table <- merge(stats_table, stats[,c("Creek", "Basin", "year", "LDNe_0.05")], by = c("Creek", "Basin", "year"), all.x = T)
stats_table <- merge(stats_table, stats[,c("Creek", "Basin", "year", "n")], by = c("Creek", "Basin", "year"))

#write.table(stats_table, "hyb_corrected_diversity_estimates_no_sib_removal.txt", sep = "\t", col.names = T, row.names = F, quote = F)
openxlsx::writeData(wb, "genetic_stats", x = stats_table, keepNA = T)

#===============hyb correction plot==========
mst <- reshape2::melt(stats_table[,-which(colnames(stats_table) == "n")], id.vars = c("Creek", "Basin", "year", "Relative.percent.RBT"))
mst$type <- ifelse(grepl("pred_", mst$variable), "corrected", "observed")
mst$variable <- gsub("pred_", "", mst$variable)
mst <- reshape2::dcast(mst, Creek + year + Basin  + variable + Relative.percent.RBT ~ type, value.var = "value")

colnames(mst)[4] <- "stat"
p <- ggplot(mst[-which(mst$stat == "LDNe_0.05"),], aes(x = observed, y = corrected, color = Relative.percent.RBT)) + 
  facet_wrap(~stat, scales = "free") + geom_point() +
  geom_abline(slope = 1, intercept = 0) + theme_bw() + 
  scale_color_viridis_c(direction = -1, end = .9) + labs(color = "%RBT") +
  theme(strip.background = element_blank(), strip.text = element_text(size = 20))

ggsave("hyb_correction_plot_no_sib_removal.pdf", p, "pdf", width = 11, height = 8.5)
p
openxlsx::insertPlot(wb, "Hyb Correction plot", width = 11, height = 8.5)

#===============theta plot==================
theta_pred <- stats_table[, c("Creek", "Basin", "year",
                              "pred_Tajimas.theta.4Nu", "pred_Wattersons.theta.4Nu", "pred_theta_diff")]
theta_pred$type <- "corrected"
colnames(theta_pred) <- gsub("pred_", "", colnames(theta_pred))


theta_obs <- stats_table[, c("Creek", "Basin", "year",
                             "Tajimas.theta.4Nu", "Wattersons.theta.4Nu", "theta_diff")]
theta_obs$type <- "observed"
theta_plot_dat <- rbind(theta_pred, theta_obs)
theta_plot_dat <- merge(theta_plot_dat, stats_table[,c("Basin", "Creek", "year", "Relative.percent.RBT")],
                        by = c("Basin", "Creek", "year"))

## merge and plot
p2 <- ggplot(theta_plot_dat, 
             aes(x = Tajimas.theta.4Nu, y = Wattersons.theta.4Nu, size = theta_diff, color = Relative.percent.RBT,
                 fill = Relative.percent.RBT)) +
  geom_point(alpha = 0.5) + geom_abline(slope = 1, intercept = 0) + facet_wrap(~type, ncol = 1) + theme_bw() +
  theme(strip.background = element_blank()) + ylab("Watterson's Theta") + xlab("Tajima's Theta") +
  scale_color_viridis_c() + scale_fill_viridis_c()

ggsave("theta_plot_no_sib_removal.pdf", p2, "pdf", width = 11, height = 8.5)
p2
openxlsx::insertPlot(wb, "Tajima's D plot", width = 11, height = 8.5)


#================rf for PVA================
rfstats <- stats


# remove anything where the genetic year varies too strongly from demographic year
missmatch <- which(abs(as.numeric(rfstats$year) - rfstats$PVA.Abund.Year) > 5)
rfstats <- rfstats[-missmatch,]

# clean
rfstats <- rfstats[-which(is.na(rfstats$PVA.Extinction)),]
rfstats$Basin <- as.factor(rfstats$Basin)
rfstats$Creek <- as.factor(rfstats$Creek)
rfstats$PVA.Extinction <- scale(rfstats$PVA.Extinction)
rfstats$PVA.Abundance.2.5 <- scale(rfstats$PVA.Abundance.2.5)
rfstats$pi <- scale(rfstats$pi)
rfstats$ho <- scale(rfstats$ho)
rfstats$Wattersons.theta.4Nu <- scale(rfstats$Wattersons.theta.4Nu)
rfstats$Tajimas.theta.4Nu <- scale(rfstats$Tajimas.theta.4Nu)
rfstats$ne.rand[is.infinite(rfstats$ne.rand)] <- NA
rfstats$ne.rand <- scale(rfstats$ne.rand)




# run a random forest against extinction risk
rfstats_ext <- rfstats[,c("PVA.Extinction", "Het.Hom", "ho", "pi", "Wattersons.theta.4Nu", "Tajimas.theta.4Nu", 
                    "Relative.percent.RBT", "ne.rand","Basin", "Creek")]

# full model
mod <- ranger::ranger(dependent.variable.name = "PVA.Extinction", 
                      data = rfstats_ext[-which(is.na(rfstats_ext$Relative.percent.RBT) | is.na(rfstats_ext$ne.rand)),], 
                      num.trees = 100000, 
                      mtry = ncol(rfstats_ext) - 1,
                      num.threads = 4, 
                      verbose = T, importance = "permutation")
mod$variable.importance
pd <- data.frame(obs = rfstats_ext[-which(is.na(rfstats_ext$Relative.percent.RBT) | is.na(rfstats_ext$ne.rand)),]$PVA.Extinction,
                 pred = mod$predictions)
prf1 <- ggplot(pd, aes(obs, pred)) + geom_point() + theme_bw()
ggsave("rf_pred_plot_no_sib_removal.pdf", prf1, "pdf", width = 11, height = 8.5)

prf1
openxlsx::insertPlot(wb, "rf Extinction plot", width = 11, height = 8.5)



# run a random forest against PVA.Abundance
rfstats_ab <- rfstats[,c("PVA.Abundance.2.5", "Het.Hom", "ho", "pi", "Wattersons.theta.4Nu", "Tajimas.theta.4Nu", 
                      "Relative.percent.RBT", "ne.rand", "Basin", "Creek")]
rfstats_ab <- na.omit(rfstats_ab)

# full model
mod <- ranger::ranger(dependent.variable.name = "PVA.Abundance.2.5", 
                      data = rfstats_ab, 
                      num.trees = 100000, 
                      mtry = ncol(rfstats_ab) - 1,
                      num.threads = 4,
                      verbose = T, importance = "permutation")
mod$variable.importance

pd <- data.frame(obs = rfstats_ab$PVA.Abundance.2.5,
                 pred = mod$predictions)
prf2 <- ggplot(pd, aes(obs, pred)) + geom_point() + theme_bw()
prf2
openxlsx::insertPlot(wb, "rf Abundance plot", width = 11, height = 8.5)


#================save============
openxlsx::saveWorkbook(wb, "no_sib_removal_stats.xlsx", overwrite = T)
