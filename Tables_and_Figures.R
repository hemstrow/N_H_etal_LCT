library(dplyr); library(ranger); library(lmerTest); library(lme4); library(ggplot2); library(openxlsx)
library(sjPlot); library(sjmisc); library(sjlabelled); library(plotly); library(snpR); library(gridExtra)
wb <- openxlsx::createWorkbook(creator = "William Hemstrom")


openxlsx::addWorksheet(wb, "Table S2")
openxlsx::addWorksheet(wb, "Table S3")
openxlsx::addWorksheet(wb, "Figure 2")
openxlsx::addWorksheet(wb, "Figure 3")
openxlsx::addWorksheet(wb, "Figure 4")
openxlsx::addWorksheet(wb, "Figure 5")



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
dat <- filter_snps(dat, HWE = 0.000001)
neSY <- calc_ne(dat, "stream.year", chr = "chr", NeEstimator_path = "D://ne_estimator/Ne2-1.exe", methods = c("LD"))
neSY <- get.snpR.stats(neSY, "stream.year", "pop")
dat <- calc_basic_snp_stats(dat, c("stream.year"))
dat <- calc_het_hom_ratio(dat, c("stream.year"))
dat <- calc_tajimas_d(dat, "stream.year")
tajimas_D <- get.snpR.stats(dat, "stream.year", "single.window")
tajimas_D$ws.theta <- tajimas_D$ws.theta/nsites # devide by total number of sequenced sites to change the denominator
tajimas_D$ts.theta <- tajimas_D$ts.theta/nsites # devide by total number of sequenced sites to change the denominator

# combine and reformat
hh <- dat@sample.stats
ss <- dat@stats
ho <- tapply(ss$ho, ss$subfacet, mean, na.rm = T)
pi <- tapply(ss$pi, ss$subfacet, mean, na.rm = T)
comb_stats <- data.frame(pop = names(ho), ho = ho, pi = pi)
ss$n <- ss$maj.count + ss$min.count
comb_stats$n <- tapply(ss$n, ss$subfacet, function(x) max(x)/2)
comb_stats$ho <- comb_stats$ho/nsites
comb_stats$pi <- comb_stats$pi/nsites


# add in hybridization info
## get hybridization info
hyb <- readr::read_csv("LCT_RBT_final hybridization 2-19-19.csv")
colnames(hyb)[1] <- "ID"
colnames(hyb)[5] <- "Relative.percent.RBT"
comb <- merge(hh, hyb, by = "ID", all.x = T) # note, this will auto-subset by whatever samples are still in they analysis! (removing poorly sequenced individuals)

# summerize by population
comb_stats$pop <- gsub("BAT", "NFB", comb_stats$pop)
comb$stream.year <- paste0(comb$stream, ".", comb$year)
heho_tab <- tapply(comb$'Het/Hom', comb$stream.year, mean)
heho_tab <- data.frame(locale = names(heho_tab), Het.Hom = as.numeric(heho_tab))
hyb_tab <- tapply(comb$Relative.percent.RBT, comb$stream.year, mean, na.omit = T)
hhyb_tab <- data.frame(locale = names(hyb_tab), Relative.percent.RBT = as.numeric(hyb_tab))
stats <- merge(heho_tab, comb_stats, by.x = "locale", by.y = "pop", all.x = T)
stats <- merge(stats, hhyb_tab, by = "locale")

# merge in tajima's D/thetas
tajimas_D <- tajimas_D[,c(2,7,8)]
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

# get into a tidy table for paper
Table.S2 <- stats[,c("locale", "Basin", "Creek", "year", "Het.Hom", "ho", "pi", "ws.theta", "ts.theta", "theta_diff",
                     "LDNe_0.05", "LDNe_lCIj_0.05", "LDNe_uCIj_0.05",
                     "LDNe_0.01", "LDNe_lCIj_0.01", "LDNe_uCIj_0.01",
                     "ne.rand", 
                     "PVA.Extinction", "PVA.Ext.L95", "PVA.Ext.U95",
                     "PVA.Abundance.50", "PVA.Abundance.2.5", "PVA.Abund.Year")]

openxlsx::writeData(wb, "Table S2", x = cbind(Table.S2, n = stats$n), keepNA = T)

#================Figure 2==================
# read in PCs
pd <- read.table("bamlist1_omy_pca.pc_scores.txt", stringsAsFactors = F)

# add hyb percent
PCA.hyb <- readr::read_csv("LCT_RBT_final hybridization 2-19-19.csv")
colnames(PCA.hyb)[1] <- "ID"
colnames(PCA.hyb)[5] <- "Relative.percent.RBT"

hyb.fix <- PCA.hyb$ID[which(!PCA.hyb$ID %in% pd$Lab)]
pd.fix <- pd$Lab[which(!pd$Lab %in% PCA.hyb$ID)]

hyb.fix <- data.frame(hyb.lab = hyb.fix, index = substr(hyb.fix, 9, 16), stringsAsFactors = F)
pd.fix <- data.frame(pd.lab = pd.fix, index = substr(pd.fix, 9, 16), stringsAsFactors = F)
fix.lab <- merge(hyb.fix, pd.fix, "index", all = T)
fix.lab <- fix.lab[-1,]

matches <- match(fix.lab$pd.lab, pd$Lab)

pd$Lab[matches] <- fix.lab$hyb.lab
pd <- merge(hyb, pd, by.x = "ID", by.y = "Lab")
pd$watershed <- gsub("_.+_.+_.+", "", pd$ID)
pd$samp <- substr(pd$ID, 1, 7)

# PC1 31.8%
# PC2 0.93%

# plot
Fig2 <- ggplot(pd, aes(x = PC1, y = PC2, color = watershed, size = Relative.percent.RBT, label = samp)) + 
  geom_label() + theme_bw() +
  scale_color_viridis_d() + ylab("PC1 (31.8%)") + xlab("PC2 (0.93%)")
# ggsave(filename = "updated_PCA_plot.pdf", device = "pdf", plot = pd, width = 11, height = 8.5)
Fig2
openxlsx::insertPlot(wb, "Figure 2", width = 11, height = 8.5)


#================Figure 3==================
# basic correlations:
## melt data
year.matches <- which(abs(Table.S2$year  - Table.S2$PVA.Abund.Year) <= 5)
cor_dat <- Table.S2[year.matches,]
cor_dat <- reshape2::melt(cor_dat[,-c(1, 12, 13, 15, 16, 19, 20, 22, 23)], id.vars = colnames(cor_dat)[c(2:4, 18, 21)])

## change variable names to be pretty
pretty_names_tab <- list("Het.Hom" = "Het/Hom",
                         "ho" = bquote(H[o]),
                         "pi" = bquote(pi),
                         "ws.theta" = bquote(theta[W]),
                         "ts.theta" = bquote(theta[T]),
                         "theta_diff" = bquote(theta[diff]),
                         "LDNe_0.05" = bquote(LDNe[0.05]),
                         "LDNe_0.01" = bquote(LDNe[0.01]),
                         "ne.rand" = bquote(Ne[Colony]))

vlabeller <- function(variable,value){
  return(pretty_names_tab[value])
}

## plots
cor_plot_ab <- ggplot(cor_dat, aes(y = log10(PVA.Abundance.50), x = value)) + 
  facet_wrap(~variable, scales = "free_x", labeller = vlabeller, strip.position = "bottom") + 
  geom_point() + theme_bw() + 
  theme(strip.background = element_blank(), axis.title.x = element_blank(), 
        strip.placement = "outside", axis.text = element_text(size = 6), panel.spacing.x = unit(1, "lines")) 



cor_plot_ext <- ggplot(cor_dat, aes(y = log10(PVA.Extinction), x = value)) + 
  facet_wrap(~variable, scales = "free_x", labeller = vlabeller, strip.position = "bottom") + geom_point() +
  theme_bw() + 
  theme(strip.background = element_blank(), axis.title.x = element_blank(), 
        strip.placement = "outside", axis.text = element_text(size = 6), panel.spacing.x = unit(1, "lines"))




# rf

rfstats <- Table.S2[year.matches,]

## clean
rfstats$Basin <- as.factor(rfstats$Basin)
rfstats$Creek <- as.factor(rfstats$Creek)
rfstats$PVA.Extinction <- scale(rfstats$PVA.Extinction)
rfstats$PVA.Abundance.50 <- scale(rfstats$PVA.Abundance.50)
rfstats$pi <- scale(rfstats$pi)
rfstats$ho <- scale(rfstats$ho)
rfstats$ws.theta <- scale(rfstats$ws.theta)
rfstats$ts.theta <- scale(rfstats$ts.theta)
rfstats$ne.rand[is.infinite(rfstats$ne.rand)] <- NA
rfstats$ne.rand <- scale(rfstats$ne.rand)




## run a random forest against extinction risk
rfstats_ext <- rfstats[,c("PVA.Extinction", "Het.Hom", "ho", "pi", "ws.theta", "ts.theta",
                          "ne.rand","Basin", "Creek")]

## full model
rf_ext <- ranger::ranger(dependent.variable.name = "PVA.Extinction", 
                        data = rfstats_ext[-which(is.na(rfstats_ext$ne.rand)),], 
                        num.trees = 100000, 
                        mtry = ncol(rfstats_ext) - 1,
                        num.threads = 4, 
                        verbose = T, importance = "permutation")

pred_rf <- data.frame(obs = rfstats_ext[-which(is.na(rfstats_ext$ne.rand)),]$PVA.Extinction,
                      pred = rf_ext$predictions, statistic = "Extinction")


## run a random forest against PVA.Abundance
rfstats_ab <- rfstats[,c("PVA.Abundance.50", "Het.Hom", "ho", "pi", "ws.theta", "ts.theta",
                         "ne.rand", "Basin", "Creek")]
rfstats_ab <- na.omit(rfstats_ab)

## full model
rf_ab <- ranger::ranger(dependent.variable.name = "PVA.Abundance.50", 
                      data = rfstats_ab, 
                      num.trees = 100000, 
                      mtry = ncol(rfstats_ab) - 1,
                      num.threads = 4,
                      verbose = T, importance = "permutation")

pred_rf <- rbind(pred_rf,
                 data.frame(obs = rfstats_ab$PVA.Abundance.50,
                            pred = rf_ab$predictions,
                            statistic = "Abundance"))


## plot rf data
rf_p <- ggplot(pred_rf, aes(x = obs, y = pred)) + geom_point() +
  theme_bw() + facet_wrap(~statistic, scales = "free") + theme(strip.background = element_blank()) +
  ylab("Predicted") + xlab("Observed")


# combine plots
grid.arrange(cor_plot_ext, cor_plot_ab, rf_p, layout_matrix = matrix(c(1,2, 3, 3), byrow = T, nrow = 2),
             heights = c(1, .5))
openxlsx::insertPlot(wb, "Figure 3", width = 11, height = 8.5)


#===============Table S3, part 1================
mod_stats <- Table.S2[year.matches,]
mod_stats <- mod_stats[,-c(1, 12, 13, 15, 16, 19, 20, 22, 23)]
run_stats <- colnames(mod_stats)[-c(1:3, 13, 14)]

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
change_dat <- reshape2::melt(change_dat[,-c(1, 12, 13, 15, 16, 19, 20, 22, 23)], id.vars = colnames(change_dat)[c(2:4, 18, 21)])

change_dat$locale <- paste0(change_dat$Basin, "_", change_dat$Creek)
min_year <- tapply(change_dat$year, change_dat$locale, min, na.rm = T)
change_dat$min_year <- min_year[match(change_dat$locale, names(min_year))]
change_dat$year_delta <- change_dat$year - change_dat$min_year

change_plot <- ggplot(change_dat, aes(x = year_delta, y = value, color = locale)) +
  theme_bw() + facet_wrap(~variable, scales = "free_y", strip.position = "left", labeller = vlabeller) +
  geom_line(aes(group = locale), alpha = 0.5) +
  geom_point() + 
  theme(legend.position = "none", strip.background = element_blank(), axis.title.y = element_blank(), strip.placement = "outside") +
  scale_color_viridis_d() + xlab("Years Since First Sample")
change_plot
openxlsx::insertPlot(wb, "Figure 4", width = 11, height = 8.5)


#==============Figure 5==================
change_dat <- Table.S2[year.matches,]
change_dat <- reshape2::melt(change_dat[,-c(1, 12, 13, 15, 16, 19, 20, 22, 23)], id.vars = colnames(change_dat)[c(2:4, 18, 21)])

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
  
  stat_deltas <- rbind(stat_deltas, cbind(change_dat[matches[1],1:2], locale = locales[i], variable = comp_stats[,1], difference = delta, year_delta = ymax - ymin, PVA.Abundance.Delta = ab.delta))
}

stat_deltas <- merge(stat_deltas, Table.S2[,c("Creek", "Basin", "PVA.Extinction")], by = c("Creek", "Basin"))
colnames(stat_deltas)[4] <- c("statistic")

stat_deltas_m <- reshape2::melt(stat_deltas, id.vars = colnames(stat_deltas)[-c(7:8)])



pretty_names_tab2 <- list("Het.Hom" = "Het/Hom",
                         "ho" = bquote(H[o]),
                         "pi" = bquote(pi),
                         "ws.theta" = bquote(theta[W]),
                         "ts.theta" = bquote(theta[T]),
                         "theta_diff" = bquote(theta[diff]),
                         "LDNe_0.05" = bquote(LDNe[0.05]),
                         "LDNe_0.01" = bquote(LDNe[0.01]),
                         "ne.rand" = bquote(Ne[Colony]))

stat_deltas_m$statistic <- factor(stat_deltas_m$statistic, levels = unique(stat_deltas_m$statistic),
                                  labels = unlist(pretty_names_tab2[match(unique(stat_deltas_m$statistic), names(pretty_names_tab2))]))
stat_deltas_m$variable <- as.character(stat_deltas_m$variable)
stat_deltas_m$variable[which(stat_deltas_m$variable == "PVA.Extinction")] <- "Extinction"
stat_deltas_m$variable[which(stat_deltas_m$variable == "PVA.Abundance.Delta")] <- "Delta* Abundance"


ext_plot <- ggplot(stat_deltas_m, aes(y = value, x = difference/year_delta))+
  geom_point() +
  facet_grid(variable~statistic, scales = "free", labeller = label_parsed, switch = "y") + theme_bw() + 
  theme(strip.background = element_blank(), axis.title.y = element_blank(), strip.placement = "outside", axis.text.x = element_text(size = 4)) +
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

run_stats <- run_stats[-c(7:8)]
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
Table.S3 <- rbind(Table.S3, models)
colnames(Table.S3) <- c("Statistic", "Extinction", "Abundance", "Extinction Model")

# save
openxlsx::writeData(wb, "Table S3", x = Table.S3, keepNA = T)

#================save============
openxlsx::saveWorkbook(wb, "Tables_and_figures.xlsx", overwrite = T)

