library(dplyr); library(ranger); library(lmerTest); library(lme4); library(ggplot2); library(openxlsx)
library(sjPlot); library(sjmisc); library(sjlabelled); library(plotly); library(gridExtra);library(snpR)
wb <- openxlsx::createWorkbook(creator = "William Hemstrom")


openxlsx::addWorksheet(wb, "Table S2")
openxlsx::addWorksheet(wb, "Table S3")
openxlsx::addWorksheet(wb, "Figure 2")
openxlsx::addWorksheet(wb, "Figure 3")
openxlsx::addWorksheet(wb, "Figure 4")
openxlsx::addWorksheet(wb, "Figure 5")
openxlsx::addWorksheet(wb, "thetas")
openxlsx::addWorksheet(wb, "per-stat comparison")
openxlsx::addWorksheet(wb, "Correction model AICs")




nsites <- 60589






#=====================prepare data and calculate basic statistics/Table S2=============
# prepare data
sample_meta <- read.table("metadata.txt", header = T)
genos <- read.table("bamlist2_baits_geno.geno", stringsAsFactors = F)
sample_meta$stream <- paste0("str_", sample_meta$stream)
sample_meta$stream[sample_meta$watershed == "LHU" & sample_meta$stream == "str_INC"] <- "str_INK"
snp_meta <- genos[,1:2]
colnames(snp_meta) <- c("chr", "position")
dat <- import.snpR.data(genos[,-c(1:2)], snp_meta, sample_meta)

# calculate basic stats
# filter via HWE in every pop (removed if out of HWE in any pop)
dat <- calc_hwe(dat, "stream.year")
hwe <- get.snpR.stats(dat, "stream.year", "hwe")
hwe$single$low.p <- ifelse(hwe$single$pHWE <= 0.000001, 1, 0)
bad.loci <- tapply(hwe$single$low.p, hwe$single[,c("chr", "position")], sum, na.rm = T)
bad.loci <- reshape2::melt(bad.loci)
bad.loci <- na.omit(bad.loci)
bad.loci <- bad.loci[which(bad.loci$value > 0),] # four snps to remove
good.loci <- which(!do.call(paste0, snp.meta(dat)[,1:2]) %in% do.call(paste0, bad.loci[,1:2]))
dat <- import.snpR.data(genos[good.loci,-c(1:2)], snp_meta[good.loci,], sample_meta)

# datf <- filter_snps(dat, HWE = 0.000001)
# dat <- filter_snps(dat, HWE = 0.000001, min_ind = .6) # alt filtering
neSY <- calc_ne(dat, "stream.year", chr = "chr", NeEstimator_path = "C://usr/bin/Ne2-1.exe", methods = c("LD"))
neSY <- get.snpR.stats(neSY, "stream.year", "pop")
dat <- calc_basic_snp_stats(dat, c("stream.year"))
dat <- calc_het_hom_ratio(dat, c("stream.year"))
dat <- calc_tajimas_d(dat, "stream.year.chr", par = 6)
tajimas_D <- get.snpR.stats(dat, "stream.year.chr", "weighted.means")
tajimas_D <- dplyr::filter(tajimas_D, snp.subfacet == ".base")
tajimas_D <- tajimas_D[,-c(1, 3, 4, 7)]
colnames(tajimas_D)[2:3] <- gsub("weighted_mean_", "", colnames(tajimas_D)[2:3])
tajimas_D$ws.theta <- tajimas_D$ws.theta/nsites # devide by total number of sequenced sites to change the denominator
tajimas_D$ts.theta <- tajimas_D$ts.theta/nsites # devide by total number of sequenced sites to change the denominator

# dat <- calc_weighted_stats(dat, "stream.year")
# check <- get.snpR.stats(dat, "stream.year", "pop")
# check$weighted_mean_pi <- check$weighted_mean_pi/nsites
# check$weighted_mean_ho <- check$weighted_mean_ho/nsites

# combine and reformat
hh <- dat@sample.stats
ss <- dat@stats
ss$n <- ss$maj.count + ss$min.count
comb_stats <- tapply(ss$n, ss$subfacet, function(x) max(x)/2)
comb_stats <- data.frame(subfacet = names(comb_stats), n = as.numeric(comb_stats))
stats <- get.snpR.stats(dat, "stream.year", c("ho", "pi"))
stats <- na.omit(stats$weighted.means)
colnames(stats) <- gsub("weighted_mean_", "", colnames(stats))
stats <- stats[,-c(3:4)]
comb_stats <- merge(comb_stats, stats, by = "subfacet", all = T)
comb_stats <- dplyr::filter(comb_stats, subfacet != ".base")

comb_stats$ho <- comb_stats$ho/nsites
comb_stats$pi <- comb_stats$pi/nsites
comb_stats <- comb_stats[,-3]
colnames(comb_stats)[1] <- "pop"


# add in hybridization info
## get hybridization info
hyb <- openxlsx::read.xlsx("Table S1 20210111.xlsx")
hyb_key <- openxlsx::read.xlsx("Table S1 20210111.xlsx", sheet = 2)
hyb <- merge(hyb, hyb_key, by.x = "Previous.ID", by.y = "New.ID")
colnames(hyb)[7] <- "ID"
colnames(hyb)[4] <- "Relative.percent.RBT"
comb <- merge(hh, hyb, by = "ID", all.x = T) # note, this will auto-subset by whatever samples are still in they analysis! (removing poorly sequenced individuals)

# summerize by population
comb_stats$pop <- gsub("BAT", "NFB", comb_stats$pop)
comb$stream.year <- paste0(comb$stream, ".", comb$year)
heho_tab <- get.snpR.stats(dat, "stream.year", "het_hom_ratio")$weighted.means
heho_tab <- na.omit(heho_tab)

heho_tab <- data.frame(locale = heho_tab$subfacet, Het.Hom = heho_tab$`weighted_mean_Het/Hom`)
hyb_tab <- tapply(comb$Relative.percent.RBT, comb$stream.year, mean, na.omit = T)
hhyb_tab <- data.frame(locale = names(hyb_tab), Relative.percent.RBT = as.numeric(hyb_tab))
stats <- merge(heho_tab, comb_stats, by.x = "locale", by.y = "pop", all.x = T)
stats <- merge(stats, hhyb_tab, by = "locale")

# merge in tajima's D/thetas
colnames(tajimas_D)[1] <- "locale"
stats <- merge(stats, tajimas_D)
stats$theta_diff <- (stats$ts.theta - stats$ws.theta)/rowMeans(stats[,c("ts.theta", "ws.theta")])

# merge in LDNe info
neSY <- neSY[,-1]
colnames(neSY)[1] <- "locale"
LDNe <- neSY
rm(neSY)
LDNe[LDNe == Inf] <- NA
stats <- merge(stats, LDNe, "locale", all.x = T)

# merge in colony Ne info
ne_colony <- read.table("colony_ne.txt", header = T)
ne_colony$pop <- paste0(substr(ne_colony$pop, 1, 7), ".", substr(ne_colony$pop, 9, 12))
colnames(ne_colony)[1] <- "locale"
stats <- merge(stats, ne_colony, by = "locale", all.x = T)

# merge in PVA data
PVA <- readr::read_csv("Table_Thetas_hyb_PVA_03062020.csv")
PVA$`PVA Extinction` <- as.numeric(gsub("%", "", PVA$`PVA Extinction`))
PVA$`PVA Ext U95` <- as.numeric(gsub("%", "", PVA$`PVA Ext U95`))
PVA$`PVA Ext L95` <- as.numeric(gsub("%", "", PVA$`PVA Ext L95`))
PVA$locale <- paste0(PVA$Stream_code, ".", PVA$`Genetic Year`)

stats <- merge(stats, PVA, by = "locale", all = T)


# fix column names
colnames(stats) <- gsub(" ", ".", colnames(stats))
colnames(stats) <- gsub("/", ".", colnames(stats))
colnames(stats) <- gsub("%", "percent", colnames(stats))
colnames(stats) <- gsub("’", "", colnames(stats))
colnames(stats) <- gsub(")", "", colnames(stats))
colnames(stats) <- gsub("\\(", "", colnames(stats))
colnames(stats) <- gsub("Genetic.Year", "year", colnames(stats))


# stats2 <- read.xlsx("Tables_and_figures_corrected.xlsx", sheet = 1)
# 
# 
# s2m <- merge(stats[,1:4], stats2, by = "locale")
# plot(s2m$pi.x, s2m$pi.y)
#==============get the corrected stats using lms================
# run models
source("make_pred_and_make_AIC_comp.R")
# hh.AIC_comp <- make_AIC_comp(stats, "Het.Hom")

ho.AIC_comp <- make_AIC_comp(stats, "ho")

pi.AIC_comp <- make_AIC_comp(stats, "pi")

# note: no effect
# ne.AIC_comp <- make_AIC_comp(stats, "ne.rand")

ldne_.01.AIC_comp <- make_AIC_comp(stats, "LDNe_0.01")

thetaT.AIC_comp <- make_AIC_comp(stats, "ts.theta")

thetaW.AIC_comp <- make_AIC_comp(stats, "ws.theta")

# ldne_.05.AIC_comp <- make_AIC_comp(stats, "LDNe_0.05") # no models work and aren't singular...

## write the AIC table
AIC_tab <- dplyr::bind_rows(ho.AIC_comp$AIC, pi.AIC_comp$AIC, thetaT.AIC_comp$AIC,
                            thetaW.AIC_comp$AIC)

openxlsx::writeData(wb, "Correction model AICs", x = AIC_tab, keepNA = T)

# predict to get values
# hh.pred <- make_pred(hh.AIC_comp)
# ggplot(hh.pred$pre, aes(x = Relative.percent.RBT, y = Het.Hom - pred_Het.Hom)) + geom_point() + theme_bw()

ho.pred <- make_pred(ho.AIC_comp)
ggplot(ho.pred$pred, aes(x = Relative.percent.RBT, y = ho - pred_ho)) + geom_point() + theme_bw()

pi.pred <- make_pred(pi.AIC_comp)
ggplot(pi.pred$pred, aes(x = Relative.percent.RBT, y = pi - pred_pi)) + geom_point() + theme_bw()

# note: no effect
# ne.pred <- make_pred(ne.AIC_comp)
# ggplot(ne.pred$pred, aes(x = Relative.percent.RBT, y = ne.rand - pred_ne.rand)) + geom_point() + theme_bw()

# ldne_.01.pred <- make_pred(ldne_.01.AIC_comp)
# ggplot(ldne_.01.pred$pred, aes(x = Relative.percent.RBT, y = LDNe_0.01 - pred_LDNe_0.01)) + geom_point() + theme_bw()

# ldne_.05.pred <- make_pred(ldne_.05.AIC_comp)
# ggplot(ldne_.05.pred$pred, aes(x = Relative.percent.RBT, y = LDNe_0.05 - pred_LDNe_0.05)) + geom_point() + theme_bw()


thetaT.pred <- make_pred(thetaT.AIC_comp)
ggplot(thetaT.pred$pred, aes(x = Relative.percent.RBT, y = ts.theta - pred_ts.theta)) + geom_point() + theme_bw()

thetaW.pred <- make_pred(thetaW.AIC_comp)
ggplot(thetaW.pred$pred, aes(x = Relative.percent.RBT, y = ws.theta - pred_ws.theta)) + geom_point() + theme_bw()

# merge them together
merge.list <- c("Creek", "Basin", "year", "Relative.percent.RBT")
stats_table <- merge(pi.pred$pred, ho.pred$pred, by = merge.list, all = T)
# stats_table <- merge(stats_table, hh.pred$pred, by = merge.list, all = T)             
stats_table <- merge(stats_table, thetaT.pred$pred, by = merge.list, all = T)             
stats_table <- merge(stats_table, thetaW.pred$pred, by = merge.list, all = T)
stats_table$pred_theta_diff <- (stats_table$pred_ts.theta - stats_table$pred_ws.theta)/rowMeans(stats_table[,c("pred_ts.theta", "pred_ws.theta")])
stats_table$theta_diff <- (stats_table$ts.theta - stats_table$ws.theta)/rowMeans(stats_table[,c("ts.theta", "ws.theta")])
# stats_table <- merge(stats_table, ne.pred$pred, by = merge.list, all = T)
# stats_table <- merge(stats_table, ldne_.01.pred$pred, by = merge.list, all = T)
stats_table <- merge(stats_table, stats[,c("locale", "Basin", "Creek", "year",
                                           "ne.rand",
                                           "LDNe_0.01",
                                           "LDNe_0.05",
                                           "LDNe_lCIj_0.05", "LDNe_uCIj_0.05",
                                           "LDNe_lCIj_0.01", "LDNe_uCIj_0.01",
                                           "PVA.Extinction", "PVA.Ext.L95", "PVA.Ext.U95",
                                           "PVA.Abundance.50", "PVA.Abundance.2.5", "PVA.Abund.Year", "n")], by = c("Creek", "Basin", "year"), all.x = T)
colnames(stats_table) <- gsub("pred_", "corrected_", colnames(stats_table))

Table.S2 <- stats_table

stats_table <- Table.S2[,c("locale", "Basin", "Creek", "year",
                           "corrected_pi", "corrected_ho",
                           "corrected_ts.theta", "corrected_ws.theta", "corrected_theta_diff", "ne.rand",
                           "LDNe_0.01",
                           "PVA.Extinction", "PVA.Ext.L95", "PVA.Ext.U95",
                           "PVA.Abundance.50", "PVA.Abundance.2.5", "PVA.Abund.Year", "n")]
colnames(stats_table) <- gsub("corrected_", "", colnames(stats_table))

openxlsx::writeData(wb, "Table S2", x = Table.S2[-which(colnames(Table.S2) == "locale"),], keepNA = T)

#================Figure 2==================
# read in PCs
hyb$watershed <- gsub("_.+_.+_.+", "", hyb$ID)
hyb$samp <- substr(hyb$ID, 1, 7)

# PC1 31.8%
# PC2 0.93%

# plot
col.scale <- RColorBrewer::brewer.pal(5, "Dark2")
myColors <- c("black", col.scale) # LTR, MAG, MCD, RBT, SFG
hyb$use_watershed <- hyb$watershed
hyb$use_watershed[which(!hyb$use_watershed %in% c("LTR", "MAG", "MCD", "RBT", "SFH"))] <- "Other"
hyb$use_watershed <- factor(hyb$use_watershed, levels = c("Other", "LTR", "MAG", "MCD", "RBT", "SFH"))

col.scale <- scale_colour_manual(name = "use_watershed",values = myColors)

Fig2a <- ggplot(hyb, aes(x = PC1, y = PC2, color = use_watershed)) + 
  geom_point(size = 6) + theme_bw() +
  col.scale + ylab("PC1 (31.8%)") + xlab("PC2 (0.93%)") +
  ggtitle("A") + guides(color = guide_legend(title = "Watershed"))
# ggsave(filename = "updated_PCA_plot.pdf", device = "pdf", plot = pd, width = 11, height = 8.5)
Fig2a



# import and prep the NGSadmix plot
setwd("NGSadmix/")
K2 <- read.table("bamlist1_omy_adm.beagle.k2.qopt")
K3 <- read.table("bamlist1_omy_adm.beagle.k3.qopt")
K2 <- cbind(K2, K3, hyb)
K2 <- K2[which(K2$watershed %in% c("QUI", "MAG", "SFH", "LTR", "MCD", "RBT")),]
write.table(K2[,1:2], "K2.sub.qopt", col.names = F, row.names = F, quote = F)
write.table(K2[,3:5], "K3.sub.qopt", col.names = F, row.names = F, quote = F)
facet.levels <- c(K2$watershed, "MAG", "SFH", "LTR", "MCD", "RBT")
facet.levels <- facet.levels[-which(duplicated(facet.levels, fromLast = T))]
NGSadmix <- plot_structure(".sub.qopt", k = 3, facet = K2$watershed, 
                           facet.order = facet.levels, qsort_K = 2, qsort = 2, 
                           separator_thickness = .5, separator_color = "white")
NGSadmix$plot <- NGSadmix$plot + xlab("Watershed") + ggtitle("B")

grid.arrange(Fig2a, NGSadmix$plot)
openxlsx::insertPlot(wb, "Figure 2", width = 11, height = 8.5)
setwd("..")

#================Figure 3==================
# basic correlations:
## melt data
year.matches <- which(abs(stats_table$year  - stats_table$PVA.Abund.Year) <= 5)
cor_dat <- stats_table[year.matches,]
cor_dat <- reshape2::melt(cor_dat[,c("locale", "Basin", "Creek", "year", "pi", "ho", "LDNe_0.01",
                                     "ts.theta", "ws.theta", "theta_diff", "ne.rand", "PVA.Extinction", "PVA.Abundance.50")], 
                          id.vars = c("locale", "Basin", "Creek", "year", "PVA.Extinction", "PVA.Abundance.50"))

## change variable names to be pretty
pretty_names_tab <- c(ho = "H[o]",
                      pi = "pi",
                      ws.theta = "theta[W]",
                      ts.theta = "theta[T]",
                      theta_diff = "theta[diff]",
                      LDNe_0.01 = "LDNe[0.01]",
                      LDNe_0.05 = "LDNe[0.05]",
                      ne.rand = "Ne[Colony]")


vlabeller <- as_labeller(pretty_names_tab, default = label_parsed)

## plots
cor_plot_ab <- ggplot(cor_dat, aes(y = log10(PVA.Abundance.50), x = value)) + 
  facet_wrap(~variable, scales = "free_x", labeller = vlabeller, strip.position = "bottom") + 
  geom_point() + theme_bw() + 
  theme(strip.background = element_blank(), axis.title.x = element_blank(), 
        strip.placement = "outside", axis.text.x = element_text(size = 10, angle = 90), 
        strip.text = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        panel.spacing.x = unit(1, "lines")) + ggtitle("A") + ylab("log10(Estimated PVA Abundance)")



cor_plot_ext <- ggplot(cor_dat, aes(y = log10(PVA.Extinction), x = value)) + 
  facet_wrap(~variable, scales = "free_x", labeller = vlabeller, strip.position = "bottom") + geom_point() +
  theme_bw() + 
  theme(strip.background = element_blank(), axis.title.x = element_blank(), 
        strip.placement = "outside", axis.text.x = element_text(size = 10, angle = 90), 
        strip.text = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        panel.spacing.x = unit(1, "lines")) + ggtitle("B") + ylab("log10(% Estimated PVA Extinction Risk)")




# rf

rfstats <- stats_table[year.matches,]

## clean
rfstats$Basin <- as.factor(rfstats$Basin)
rfstats$Creek <- as.factor(rfstats$Creek)
rfstats$PVA.Extinction <- scale(rfstats$PVA.Extinction)
rfstats$PVA.Abundance.50 <- scale(rfstats$PVA.Abundance.50)
rfstats$pi <- scale(rfstats$pi)
rfstats$ho <- scale(rfstats$ho)
rfstats$ws.theta <- scale(rfstats$ws.theta)
rfstats$ts.theta <- scale(rfstats$ts.theta)




## run a random forest against extinction risk
rfstats_ext <- rfstats[,c("PVA.Extinction","ne.rand", "ho", "pi", "ws.theta", "ts.theta")]

## full model
rf_ext <- ranger::ranger(dependent.variable.name = "PVA.Extinction", 
                         data = rfstats_ext, 
                         num.trees = 100000, 
                         mtry = ncol(rfstats_ext) - 1,
                         num.threads = 4, 
                         verbose = T, importance = "permutation")

pred_rf <- data.frame(obs = rfstats_ext$PVA.Extinction,
                      pred = rf_ext$predictions, statistic = "Extinction")

# rsq 0.05

## run a random forest against PVA.Abundance
rfstats_ab <- rfstats[,c("PVA.Abundance.50","ne.rand", "ho", "pi", "ws.theta", "ts.theta")]
rfstats_ab <- na.omit(rfstats_ab)



## full model
rf_ab <- ranger::ranger(dependent.variable.name = "PVA.Abundance.50", 
                        data = rfstats_ab, 
                        num.trees = 100000, 
                        mtry = ncol(rfstats_ab) - 1,
                        num.threads = 4,
                        verbose = T, importance = "permutation")

# rsq -.46

pred_rf <- rbind(pred_rf,
                 data.frame(obs = rfstats_ab$PVA.Abundance.50,
                            pred = rf_ab$predictions,
                            statistic = "Abundance"))


## plot rf data
pred_rf$new_statistic <- ifelse(pred_rf$statistic == "Abundance", "log10(Estimated PVA Abundance)", "log10(% Estimated PVA Extinction Risk)")
pred_rf$new_statistic <- factor(pred_rf$new_statistic, levels = c("log10(Estimated PVA Abundance)", "log10(% Estimated PVA Extinction Risk)"))
rf_p <- ggplot(pred_rf, aes(x = obs, y = pred)) + geom_point() +
  theme_bw() + facet_wrap(~new_statistic, scales = "free") + theme(strip.background = element_blank()) +
  ylab("Predicted") + xlab("Observed") + ggtitle("C")


# combine plots
grid.arrange(cor_plot_ab + theme(plot.margin = unit(c(.4, .4, .4, .4), "cm")), 
             cor_plot_ext + theme(plot.margin = unit(c(.4, .4, .4, .4), "cm")), 
             rf_p + theme(plot.margin = unit(c(.4, .4, .4, .4), "cm")),
             layout_matrix = matrix(c(1,2, 3, 3), byrow = T, nrow = 2),
             heights = c(1, .5))
openxlsx::insertPlot(wb, "Figure 3", width = 11, height = 8.5)


#===============Table S3, part 1================
mod_stats <- stats_table[year.matches,]
mod_stats <- mod_stats[,c("Basin", "Creek", "year", "pi", "ho", "Het.Hom",
                          "ts.theta", "ws.theta", "theta_diff", "ne.rand", "PVA.Extinction", "PVA.Abundance.50")]
run_stats <- c("pi", "ho", "ts.theta", "ws.theta", "theta_diff", "ne.rand")

# clean
mod_stats$Basin <- as.factor(mod_stats$Basin)
mod_stats$Creek <- as.factor(mod_stats$Creek)
mod_stats$PVA.Abundance.50 <- scale(mod_stats$PVA.Abundance.50)
mod_stats$PVA.Extinction <- scale(mod_stats$PVA.Extinction)

models <- data.frame(stat = numeric(length(run_stats)),
                     Ext = numeric(length(run_stats)),
                     AB = numeric(length(run_stats)),
                     Ext.model = character(length(run_stats)),
                     stringsAsFactors = F)

for(i in 1:length(run_stats)){
  # prep
  tdat <- mod_stats[,c("Creek", "Basin", "PVA.Abundance.50", "PVA.Extinction", run_stats[i])]
  colnames(tdat)[5] <- "stat"
  bads <- which(is.na(tdat$stat) | is.infinite(tdat$stat))
  if(length(bads) > 0){
    tdat <- tdat[-bads,]
  }
  tdat$stat <- scale(tdat$stat)
  
  # run models, checking if the lmer worked for ext and running a basic lm if not
  AB <- lmerTest::lmer(PVA.Abundance.50 ~ stat + (1|Basin) + (1|Basin:Creek), tdat)
  Ext <- try(lmerTest::lmer(PVA.Extinction ~ stat + (1|Basin) + (1|Basin:Creek), tdat), silent = T)
  out <- numeric(2)
  out[2] <- summary(AB)$coefficients[2,5]
  
  ## check if the lmer for Ext is good:
  goodlme <- T
  if(class(Ext) == "try-error"){goodlme <- F}
  else if(!is.null(summary(Ext)$optinfo$conv$lme4$code)){goodlme <- F}
  
  
  if(!goodlme){
    cat("Statistic: ", run_stats[i], " failed to converge with lmer.\n")
    Ext <- lm(PVA.Extinction ~ stat, tdat)
    out[1] <- summary(Ext)$coefficients[2,4]
    models[i, 4] <- "lm"
  }
  else if(goodlme){
    cat("Statistic: ", run_stats[i], " converged with lmer.\n")
    out[1] <- summary(Ext)$coefficients[2,5]
    models[i, 4] <- "glmm"
  }
  
  models[i,1] <- run_stats[i]
  models[i,2:3] <- out
}

Table.S3 <- models



#==============Figure 4==========
change_dat <- reshape2::melt(stats_table[,c("locale", "Basin", "Creek", "year", "pi", "ho",
                                            "ts.theta", "ws.theta", "theta_diff", "ne.rand", "PVA.Extinction", "PVA.Abundance.50")], 
                             id.vars = c("locale", "Basin", "Creek", "year", "PVA.Extinction", "PVA.Abundance.50"))

change_dat$locale <- paste0(change_dat$Basin, "_", change_dat$Creek)
min_year <- tapply(change_dat$year, change_dat$locale, min, na.rm = T)
change_dat$min_year <- min_year[match(change_dat$locale, names(min_year))]
change_dat$year_delta <- change_dat$year - change_dat$min_year

change_plot <- ggplot(change_dat, aes(x = year_delta, y = value, color = locale)) +
  theme_bw() + facet_wrap(~variable, scales = "free_y", strip.position = "left", labeller = vlabeller) +
  geom_line(aes(group = locale), alpha = 0.5) +
  geom_point() + 
  theme(legend.position = "none", strip.background = element_blank(), 
        axis.title.y = element_blank(), strip.placement = "outside", 
        strip.text = element_text(size = 12), axis.text = element_text(size = 10)) +
  scale_color_viridis_d() + xlab("Years Since First Sample")
change_plot
openxlsx::insertPlot(wb, "Figure 4", width = 11, height = 8.5)


#==============Figure 5==================
change_dat <- stats_table[year.matches,]
change_dat <- reshape2::melt(change_dat[,c("locale", "Basin", "Creek", "year", "pi", "ho",
                                           "ts.theta", "theta_diff", "ne.rand", "PVA.Extinction", "PVA.Abundance.50")], 
                             id.vars = c("locale", "Basin", "Creek", "year", "PVA.Extinction", "PVA.Abundance.50"))
# note: dropping ws.theta because it had a fixed effect for year

change_dat$locale <- paste0(change_dat$Basin, "_", change_dat$Creek)
min_year <- tapply(change_dat$year, change_dat$locale, min, na.rm = T)
change_dat$min_year <- min_year[match(change_dat$locale, names(min_year))]
change_dat$year_delta <- change_dat$year - change_dat$min_year

stat_deltas <- data.frame(Creek = character(0), Basin = character(0), locale = character(0), variable = character(0), difference = numeric(0), year_delta = numeric(0))
locales <- unique(change_dat$locale)
for(i in 1:length(locales)){
  matches <- which(change_dat$locale == locales[i])
  years <- unique(change_dat[matches,]$year)
  if(length(years) == 1){
    next
  }
  
  ymin <- min(years)
  ymax <- max(years)
  min_stats <- change_dat[intersect(matches, which(change_dat$year == ymin)),]
  max_stats <- change_dat[intersect(matches, which(change_dat$year == ymax)),]
  comp_stats <- merge(max_stats[,c("variable", "value", "PVA.Abundance.50")], min_stats[,c("variable", "value", "PVA.Abundance.50")],
                      by = "variable")
  delta <- comp_stats$value.y - comp_stats$value.x
  
  ab.delta <- comp_stats$PVA.Abundance.50.y - comp_stats$PVA.Abundance.50.x
  
  stat_deltas <- rbind(stat_deltas, cbind(change_dat[matches[1],2:3], locale = locales[i], variable = comp_stats[,1], difference = delta, year_delta = ymax - ymin, PVA.Abundance.Delta = ab.delta))
}

stat_deltas <- merge(stat_deltas, stats_table[,c("Creek", "Basin", "PVA.Extinction")], by = c("Creek", "Basin"))
colnames(stat_deltas)[4] <- c("statistic")

stat_deltas_m <- reshape2::melt(stat_deltas, id.vars = colnames(stat_deltas)[-c(7:8)])



# pretty_names_tab2 <- list("Het.Hom" = "Het/Hom",
#                           "ho" = bquote(H[o]),
#                           "pi" = bquote(pi),
#                           "ws.theta" = bquote(theta[W]),
#                           "ts.theta" = bquote(theta[T]),
#                           "theta_diff" = bquote(theta[diff]),
#                           "LDNe_0.01" = bquote(LDNe[0.01]),
#                           "ne.rand" = bquote(Ne[Colony]))

pretty_names_tab2 <- c(pretty_names_tab, PVA.Extinction = "Estimated~PVA~Extinction~Risk", PVA.Abundance.Delta = "Delta*~Estimated~PVA~Abundance")
vlabeller <- as_labeller(pretty_names_tab2, default = label_parsed)



# stat_deltas_m$statistic <- factor(stat_deltas_m$statistic, levels = unique(stat_deltas_m$statistic),
#                                   labels = unlist(pretty_names_tab2[match(unique(stat_deltas_m$statistic), names(pretty_names_tab2))]))
# stat_deltas_m$variable <- as.character(stat_deltas_m$variable)
# stat_deltas_m$variable[which(stat_deltas_m$variable == "PVA.Extinction")] <- "Extinction"
# stat_deltas_m$variable[which(stat_deltas_m$variable == "PVA.Abundance.Delta")] <- "Delta* Abundance"
stat_deltas_m$rate_of_change <- stat_deltas_m$difference/stat_deltas_m$year_delta

ext_plot <- ggplot(stat_deltas_m, aes(y = value, x = rate_of_change))+
  geom_point() +
  facet_grid(variable~statistic, scales = "free", labeller = vlabeller, switch = "both") + theme_bw() + 
  theme(strip.background = element_blank(), axis.title.y = element_blank(), strip.placement = "outside", 
        axis.text.x = element_text(size = 10, angle = 90), strip.text = element_text(size = 12)) +
  xlab(bquote(Delta*Statistic/Year))

ext_plot
openxlsx::insertPlot(wb, "Figure 5", width = 11, height = 8.5)


# stat_deltasc <- reshape2::dcast(stat_deltas, Creek + Basin + locale + year_delta + PVA.Extinction ~ variable, value.var = "difference", mean)

#===============Table S3 part 2================
stat_deltasc <- reshape2::dcast(stat_deltas, Creek + Basin + locale + year_delta + PVA.Extinction + PVA.Abundance.Delta~ statistic, value.var = "difference", mean)
mod_stats <- stat_deltasc

# clean
mod_stats$Basin <- as.factor(mod_stats$Basin)
mod_stats$Creek <- as.factor(mod_stats$Creek)
mod_stats$PVA.Abundance.Delta <- scale(mod_stats$PVA.Abundance.Delta)
mod_stats$PVA.Extinction <- scale(mod_stats$PVA.Extinction)

run_stats <- c("pi", "ho", "ts.theta", "theta_diff")

models <- data.frame(stat = numeric(length(run_stats)),
                     Ext = numeric(length(run_stats)),
                     AB = numeric(length(run_stats)),
                     Ext.model = "lm",
                     stringsAsFactors = F)

for(i in 1:length(run_stats)){
  # prep
  tdat <- mod_stats[,c("Creek", "Basin", "PVA.Abundance.Delta", "PVA.Extinction", run_stats[i])]
  colnames(tdat)[5] <- "stat"
  bads <- which(is.na(tdat$stat) | is.infinite(tdat$stat))
  if(length(bads) > 0){
    tdat <- tdat[-bads,]
  }
  tdat$stat <- scale(tdat$stat)
  
  # run models, checking if the lmer worked for ext and running a basic lm if not
  AB <- lm(PVA.Abundance.Delta ~ stat, tdat)
  Ext <- lm(PVA.Extinction ~ stat, tdat)
  out <- c(summary(Ext)$coefficients[2,4], summary(AB)$coefficients[2,4])
  
  models[i,1] <- paste0("Delta_", run_stats[i])
  models[i,2:3] <- out
}

# combine
colnames(models) <- c("Statistic", "Extinction", "Abundance", "Model Type")
colnames(Table.S3) <- c("Statistic", "Extinction", "Abundance", "Model Type")
Table.S3 <- rbind(Table.S3, models)

# save
openxlsx::writeData(wb, "Table S3", x = Table.S3, keepNA = T)


#================correction plot===============
mst <- Table.S2[,c("Creek", "Basin", "year", "Relative.percent.RBT", "corrected_pi", "pi", "corrected_ho", "ho",
                   "corrected_ts.theta", "ts.theta", "corrected_ws.theta", "ws.theta", "corrected_theta_diff", "theta_diff")]
mst <- reshape2::melt(mst, id.vars = c("Creek", "Basin", "year", "Relative.percent.RBT"))
mst$type <- ifelse(grepl("corrected_", mst$variable), "corrected", "observed")
mst$variable <- gsub("corrected_", "", mst$variable)
mst <- reshape2::dcast(mst, Creek + year + Basin  + variable + Relative.percent.RBT ~ type, value.var = "value")

colnames(mst)[4] <- "stat"
p <- ggplot(mst, aes(x = observed, y = corrected, color = Relative.percent.RBT)) + 
  facet_wrap(~stat, scales = "free", labeller = vlabeller, strip.position = "bottom") + 
  geom_point(size = 4) +
  geom_abline(slope = 1, intercept = 0) + theme_bw() + 
  scale_color_viridis_c(direction = 1, end = .85) + labs(color = "%RBT") +
  theme(strip.background = element_blank(), strip.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90), axis.text = element_text(size = 10), 
        strip.placement = "outside") +
  xlab("Corrected Value") + ylab("Observed Value")


p
openxlsx::insertPlot(wb, "per-stat comparison", width = 11, height = 8.5)
saveRDS(mst, "mst.RDS")

#================thetas=====================
theta_pred <- Table.S2[,c("Creek", "Basin", "year",
                          "corrected_ts.theta", "corrected_ws.theta", "corrected_theta_diff")]
theta_pred$type <- "corrected"
colnames(theta_pred) <- gsub("corrected_", "", colnames(theta_pred))


theta_obs <- Table.S2[,c("Creek", "Basin", "year",
                         "ts.theta", "ws.theta", "theta_diff")]
theta_obs$type <- "observed"
theta_plot_dat <- rbind(theta_pred, theta_obs)
theta_plot_dat <- merge(theta_plot_dat, Table.S2[,c("Basin", "Creek", "year", "Relative.percent.RBT")],
                        by = c("Basin", "Creek", "year"))

theta_plot_dat$samp <- paste0(theta_plot_dat$Basin, "_", theta_plot_dat$Creek, "_", theta_plot_dat$year)


## merge and plot
p2 <- ggplot(theta_plot_dat, 
             aes(x = ts.theta, y = ws.theta, size = theta_diff, color = Relative.percent.RBT,
                 fill = Relative.percent.RBT, label = samp)) +
  geom_point(alpha = 0.5) + geom_abline(slope = 1, intercept = 0) + facet_wrap(~type, ncol = 1) + theme_bw() +
  theme(strip.background = element_blank()) + ylab(expression(theta[W])) + xlab(expression(theta[T])) +
  scale_color_viridis_c() + scale_fill_viridis_c()

p2

saveRDS(theta_plot_dat, "theta_plot_dat.RDS")
openxlsx::insertPlot(wb, "thetas", width = 11, height = 8.5)


#================save============
openxlsx::saveWorkbook(wb, "Tables_and_figures_corrected.xlsx", overwrite = T)

