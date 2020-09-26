library(dplyr); library(ranger); library(lmerTest); library(lme4); library(ggplot2)
library(sjPlot); library(sjmisc); library(sjlabelled); library(plotly)

#=================prep data=====================

stats <- read.csv("stream_year_combined_stats.csv")

# prep and merge the PVA data
PVA <- readr::read_csv("Table_Thetas_hyb_PVA_03062020.csv")
PVA$`PVA Extinction` <- as.numeric(gsub("%", "", PVA$`PVA Extinction`))
PVA$`PVA Ext U95` <- as.numeric(gsub("%", "", PVA$`PVA Ext U95`))
PVA$`PVA Ext L95` <- as.numeric(gsub("%", "", PVA$`PVA Ext L95`))

PVA$Group <- gsub("LHU_INC", "LHU_INK", PVA$Group)
PVA$locale <- gsub("(.+)_(.+)_(.+)_(.+)", "str_\\3\\.\\4", PVA$Group)
PVA$locale <- gsub("BAT", "NFB", PVA$locale)
stats <- merge(stats, PVA, by = "locale", all = T)
colnames(stats) <- gsub(" ", ".", colnames(stats))
colnames(stats) <- gsub("/", ".", colnames(stats))
colnames(stats) <- gsub("%", "percent", colnames(stats))
colnames(stats) <- gsub("’", "", colnames(stats))
colnames(stats) <- gsub(")", "", colnames(stats))
colnames(stats) <- gsub("\\(", "", colnames(stats))

## add in the updated thetas from sibpurge
spt <- read.table("sibpurge_bamlists/cat_indF_thetas.txt")
colnames(spt) <- c("locale", "Wattersons.theta.4Nu", "Tajimas.theta.4Nu")
spt$locale <- gsub("sibpurge_bamlist_", "", spt$locale)
spt$locale <- gsub("_baits_thetas-indF2.thetas.global", "", spt$locale)
spt$locale <- gsub("LHU_INC", "LHU_INK", spt$locale)
spt$locale <- gsub("(.+)_(.+)_(.+)", "str_\\2\\.\\3", spt$locale)
spt$locale <- gsub("BAT", "NFB", spt$locale)
stats$Wattersons.theta.4Nu <- NULL
stats$Tajimas.theta.4Nu <- NULL
stats <- merge(stats, spt, by = "locale", all.x = T)

## Add in the Ne values from LDNe
LDNe <- read.table("ne_sib_removed.csv", header = T, sep = ",")
LDNe <- LDNe[which(LDNe$facet == "stream.year"),]
LDNe$pop <- gsub("BAT", "NFB", LDNe$pop)
colnames(LDNe)[1] <- "locale"
LDNe <- LDNe[,c(1,2,9)]
LDNe$LDNe[is.infinite(LDNe$LDNe)] <- NA
LDNec <- reshape2::dcast(LDNe, locale ~ pcrit, value.var = "LDNe")
colnames(LDNec) <- c("locale", paste0("LDNe_", colnames(LDNec)[-1]))
stats <- merge(stats, LDNec, by = "locale", all.x = T)

ostats <- stats


## zero out cases where the genetics and PVA years aren't within 5 years
missmatch <- which(abs(as.numeric(stats$year) - stats$PVA.Abund.Year) > 5)
stats <- stats[-missmatch,]



#================glms========================

mod1 <- lmerTest::lmer(Het.Hom ~ (1|Basin/Creek) + year + Relative.percent.RBT, ostats)
mod2 <- lmerTest::lmer(Het.Hom ~ (1|Basin/Creek) + (1|year) + Relative.percent.RBT, ostats) # best
mod3 <- lmerTest::lmer(Het.Hom ~ (1|Basin/Creek)  + Relative.percent.RBT, ostats)


mod1.pi <- lmerTest::lmer(pi ~ (1|Basin/Creek) + year + Relative.percent.RBT, ostats) # use for consistancy
mod2.pi <- lmerTest::lmer(pi ~ (1|Basin/Creek) + (1|year) + Relative.percent.RBT, ostats) # singular...
mod3.pi <- lmerTest::lmer(pi ~ (1|Basin/Creek) + as.factor(year) + Relative.percent.RBT, ostats)
mod4.pi <- lmerTest::lmer(pi ~ (1|Basin/Creek) + Relative.percent.RBT, ostats) # best other

mod1.ho <- lmerTest::lmer(ho ~ (1|Basin/Creek) + year + Relative.percent.RBT, ostats)
mod2.ho <- lmerTest::lmer(ho ~ (1|Basin/Creek) + (1|year) + Relative.percent.RBT, ostats) # singular
mod3.ho <- lmerTest::lmer(ho ~ (1|Basin/Creek) + Relative.percent.RBT, ostats) # doesn't converge


mod1.wt <- lmerTest::lmer(scale(Wattersons.theta.4Nu) ~ (1|Basin/Creek) + year + Relative.percent.RBT, ostats)
mod2.wt <- lmerTest::lmer(scale(Wattersons.theta.4Nu) ~ (1|Basin/Creek) + (1|year) + Relative.percent.RBT, ostats) # best
mod3.wt <- lmerTest::lmer(scale(Wattersons.theta.4Nu) ~ (1|Basin/Creek) + Relative.percent.RBT, ostats)


mod1.tt <- lmerTest::lmer(scale(Tajimas.theta.4Nu) ~ (1|Basin/Creek) + year + Relative.percent.RBT, ostats)
mod2.tt <- lmerTest::lmer(scale(Tajimas.theta.4Nu) ~ (1|Basin/Creek) + (1|year) + Relative.percent.RBT, ostats) # best
mod3.tt <- lmerTest::lmer(scale(Tajimas.theta.4Nu) ~ (1|Basin/Creek)  + Relative.percent.RBT, ostats)

mod1.td <- lmerTest::lmer(Theta.deviation ~ (1|Basin/Creek) + year + Relative.percent.RBT, ostats)
mod2.td <- lmerTest::lmer(Theta.deviation ~ (1|Basin/Creek) + (1|year) + Relative.percent.RBT, ostats) # best
mod3.td <- lmerTest::lmer(Theta.deviation ~ (1|Basin/Creek) + Relative.percent.RBT, ostats)

mod1.ne <- lmerTest::lmer(ne.rand ~ (1|Basin) + Creek  + year + Relative.percent.RBT, 
                          ostats[-which(is.infinite(ostats$ne.rand)),]) # best, only one that isn't singular and converges

mod2.ne <- lmerTest::lmer(ne.rand ~ (1|Basin/Creek) + (1|year) + Relative.percent.RBT, ostats[-which(is.infinite(ostats$ne.rand)),])
mod3.ne <- lmerTest::lmer(ne.rand ~ (1|Basin/Creek) + Relative.percent.RBT, ostats[-which(is.infinite(ostats$ne.rand)),])
mod4.ne <- lmerTest::lmer(ne.rand ~ (1|Basin) + Creek + Relative.percent.RBT, ostats[-which(is.infinite(ostats$ne.rand)),])
mod5.ne <- lm(ne.rand ~ Basin + Creek + Relative.percent.RBT, ostats[-which(is.infinite(ostats$ne.rand)),])


mod1.ne <- lmerTest::lmer(ne.rand ~ (1|Basin/Creek) + year + Relative.percent.RBT, ostats[-which(is.infinite(ostats$ne.rand)),])


# pretty output
tab_model(mod2,
          mod1.pi,
          mod1.ho,
          mod2.wt,
          mod2.tt,
          show.aic = T, 
          digits.re = 7, 
          collapse.ci = T,
          dv.labels = c("Het:Hom", "pi", "Ho", "scale(Watterson's Theta)",
                        "scale(Tajima's Theta)"), 
          pred.labels = c("(Intercept)", "%RBT"), show.icc = F, file = "Model_summary.html", CSS = css_theme("cells"))

#===========================random forest on extinction risk=================
# define function
# cvloo <- function(rfstats, response, trees){
#   cv.out <- matrix(NA, nrow(rfstats), 2)
#   colnames(cv.out) <- c("obs", "predicted")
#   for(i in 1:nrow(rfstats)){
#     cat("CV run: ", i, "\n")
#     
#     holdout <- i
#     hmod <- ranger::ranger(dependent.variable.name = response, 
#                            data = rfstats[-holdout,], 
#                            num.trees = trees, 
#                            mtry = ncol(rfstats) - 1,
#                            num.threads = 4,
#                            verbose = T)
#     
#     hold <- rfstats[holdout,,drop = F]
#     cv.out[i,] <- c(hold[,response], predict(hmod, data = hold)$predictions) 
#   }
#   
#   return(as.data.frame(cv.out))
# }
# 
colnames(stats)[23:24] <- c("Wt", "Tt")
stats <- stats[-which(is.na(stats$PVA.Extinction)),]
stats$Basin <- as.factor(stats$Basin)
stats$Creek <- as.factor(stats$Creek)
stats$PVA.Extinction <- scale(stats$PVA.Extinction)
stats$PVA.Abundance.2.5 <- scale(stats$PVA.Abundance.2.5)
stats$pi <- scale(stats$pi)
stats$ho <- scale(stats$ho)
stats$Wt <- scale(stats$Wt)
stats$Tt <- scale(stats$Tt)
stats$ne.rand[is.infinite(stats$ne.rand)] <- NA
stats$ne.rand <- scale(stats$ne.rand)
# 
# m1 <- lmerTest::lmer(PVA.Extinction ~ (1|year) + Het.Hom + ho + pi + Wt + Tt + (1|Creek/Basin), stats)
# m2 <- lmerTest::lmer(PVA.Extinction ~ (1|year) + (1:Creek/Basin) + ne.rand + Relative.percent.RBT + ho + pi + Wt + Tt + Het.Hom, stats[-which(is.infinite(stats$ne.rand)),])
# m3 <- lmerTest::lmer(PVA.Abundance.2.5 ~ (1|year) + (1:Creek/Basin) + ne.rand + Relative.percent.RBT + ho + pi + Wt + Tt + Het.Hom, stats[-which(is.infinite(stats$ne.rand)),])
# 

# run a random forest against extinction risk
rfstats <- stats[,c("PVA.Extinction", "Het.Hom", "ho", "pi", "Wt", "Tt", "Relative.percent.RBT", "ne.rand", "Basin", "Creek")]
## cross validate, leave-one-out
# cv.er <- cvloo(rfstats, "PVA.Extinction", 100000)
# plot(cv.er$obs, cv.er$predicted)
# cor(cv.er$obs, cv.er$predicted)^2

# full model
mod <- ranger::ranger(dependent.variable.name = "PVA.Extinction", 
                      data = rfstats[-which(is.na(rfstats$Relative.percent.RBT) | is.na(rfstats$ne.rand)),], 
                      num.trees = 100000, 
                      mtry = ncol(rfstats) - 1,
                      num.threads = 4, 
                      verbose = T, importance = "permutation")
mod$variable.importance
pd <- data.frame(obs = rfstats[-which(is.na(rfstats$Relative.percent.RBT) | is.na(rfstats$ne.rand)),]$PVA.Extinction,
                 pred = mod$predictions)
ggsave("rf_pred_plot.pdf", ggplot(pd, aes(obs, pred)) + geom_point() + theme_bw(), "pdf", width = 11, height = 8.5)

# run a random forest against PVA.Abundance
rfstats <- stats[,c("PVA.Abundance.2.5", "Het.Hom", "ho", "pi", "Wt", "Tt", "Relative.percent.RBT", "ne.rand", "Basin", "Creek")]
rfstats <- na.omit(rfstats)
## cross validate, leave-one-out
#cv.er <- cvloo(rfstats, "PVA.Abundance.2.5", 100000)
plot(cv.er$obs, cv.er$predicted)
cor(cv.er$obs, cv.er$predicted)^2

# full model
mod <- ranger::ranger(dependent.variable.name = "PVA.Abundance.2.5", 
                      data = rfstats, 
                      num.trees = 100000, 
                      mtry = ncol(rfstats) - 1,
                      num.threads = 4,
                      verbose = T, importance = "permutation")
mod$variable.importance

#=========================corrections for hybrids===============
make_pred <- function(dat, vars){
  mstats <- dat[,c("Basin", "Creek", "Relative.percent.RBT", "year", vars)]
  mstats <- na.omit(mstats)
  infs <- which(is.infinite(mstats[,vars]))
  if(length(infs) > 0){
    mstats <- mstats[-infs,]
  }
  formula <- paste0(vars, "~ (1|Basin/Creek) + Relative.percent.RBT + (1|year)")
  mod <- lme4::lmer(formula, mstats)
  if(isSingular(mod)){
    formula <- paste0(vars, "~ (1|Basin/Creek) + Relative.percent.RBT + year")
    mod <- lme4::lmer(formula, mstats)
  }
  if(vars == "ne.rand"){
    formula <- paste0(vars, "~ (1|Basin) + Creek  + year + Relative.percent.RBT")
    mod <- lme4::lmer(formula, mstats)
  }
  
  unique_streams <- unique(mstats[,c("Creek", "Basin", "year")])
  unique_streams <- na.omit(cbind.data.frame(unique_streams, Relative.percent.RBT = 0))
  unique_streams$Creek <- as.character(unique_streams$Creek)
  unique_streams$Basin <- as.character(unique_streams$Basin)
  
  pred <- predict(mod, newdata = unique_streams)
  pred <- cbind(unique_streams, pred)
  colnames(pred)[ncol(pred)] <- paste0("pred_", vars)
  rownames(pred) <- 1:nrow(pred)
  pred$Relative.percent.RBT <- NULL
  pred <- merge(pred, dat[,c("Basin", "Creek", "year", "Relative.percent.RBT", vars)])

  return(list(pred = pred, mod = mod))
}

hh.pred <- make_pred(ostats, "Het.Hom")
ggplot(hh.pred$pre, aes(x = Relative.percent.RBT, y = Het.Hom - pred_Het.Hom)) + geom_point() + theme_bw()

ho.pred <- make_pred(ostats, "ho")
ggplot(ho.pred$pred, aes(x = Relative.percent.RBT, y = ho - pred_ho)) + geom_point() + theme_bw()

pi.pred <- make_pred(ostats, "pi")
ggplot(pi.pred$pred, aes(x = Relative.percent.RBT, y = pi - pred_pi)) + geom_point() + theme_bw()

# note: no effect
ne.pred <- make_pred(ostats, "ne.rand")
ggplot(ne.pred$pred, aes(x = Relative.percent.RBT, y = ne.rand - pred_ne.rand)) + geom_point() + theme_bw()

thetaT.pred <- make_pred(ostats, "Tajimas.theta.4Nu")
ggplot(thetaT.pred$pred, aes(x = Relative.percent.RBT, y = Tajimas.theta.4Nu - pred_Tajimas.theta.4Nu)) + geom_point() + theme_bw()

thetaW.pred <- make_pred(ostats, "Wattersons.theta.4Nu")
ggplot(thetaW.pred$pred, aes(x = Relative.percent.RBT, y = Wattersons.theta.4Nu - pred_Wattersons.theta.4Nu)) + geom_point() + theme_bw()


# combined table
merge.list <- c("Creek", "Basin", "year", "Relative.percent.RBT")
stats_table <- merge(pi.pred$pred, ho.pred$pred, by = merge.list, all = T)
stats_table <- merge(stats_table, hh.pred$pred, by = merge.list, all = T)             
stats_table <- merge(stats_table, thetaT.pred$pred, by = merge.list, all = T)             
stats_table <- merge(stats_table, thetaW.pred$pred, by = merge.list, all = T)
stats_table <- merge(stats_table, ne.pred$pred, by = merge.list, all = T)
stats_table$theta_diff <- (stats_table$Tajimas.theta.4Nu - stats_table$Wattersons.theta.4Nu)/rowMeans(stats_table[,c(12,14)])
stats_table$pred_theta_diff <- (stats_table$pred_Tajimas.theta.4Nu - stats_table$pred_Wattersons.theta.4Nu)/rowMeans(stats_table[,c(11,13)])
stats_table <- merge(stats_table, ostats[,c("Creek", "Basin", "year", "n")], by = c("Creek", "Basin", "year"))
stats_table <- merge(stats_table, ostats[,c("Creek", "Basin", "year", "LDNe_0.01", "LDNe_0.02", "LDNe_0.05")], by = c("Creek", "Basin", "year"), all.x = T)

write.table(stats_table, "hyb_corrected_diversity_estimates.txt", sep = "\t", col.names = T, row.names = F, quote = F)

# plot
mst <- reshape2::melt(stats_table[,-which(colnames(stats_table) == "n")], id.vars = c("Creek", "Basin", "year", "Relative.percent.RBT"))
mst$type <- ifelse(grepl("pred_", mst$variable), "corrected", "observed")
mst$variable <- gsub("pred_", "", mst$variable)
mst <- reshape2::dcast(mst, Creek + year + Basin  + variable + Relative.percent.RBT ~ type, value.var = "value")

colnames(mst)[4] <- "stat"
p <- ggplot(mst, aes(x = observed, y = corrected, color = Relative.percent.RBT)) + 
  facet_wrap(~stat, scales = "free") + geom_point() +
  geom_abline(slope = 1, intercept = 0) + theme_bw() + 
  scale_color_viridis_c(direction = -1, end = .9) + labs(color = "%RBT") +
  theme(strip.background = element_blank(), strip.text = element_text(size = 20))


ggsave(filename = "corrected_stats_plot.pdf", p,
       width = 11, height = 8.5)


# plot tajima's trajectory
## thetas
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

ggsave(filename = "corrected_thetas_plot.pdf", p2,
       width = 11, height = 8.5)

ggplotly(p2)


# Ne
nep <- ggplot(ostats[-which(is.na(ostats$Basin)),], aes(x = Creek, y = log10(ne.rand), color = as.factor(year), fill = as.factor(year))) + 
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = log10(lcl.rand), ymax = log10(ucl.rand)), position=position_dodge(width=0.5)) + 
  scale_y_continuous(limits = c(0, 3)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90), strip.background = element_blank()) +
  facet_wrap(~Basin, scales = "free_x")

nepnl <- ggplot(ostats[-which(is.na(ostats$Basin)),], aes(x = Creek, y = ne.rand, color = as.factor(year), fill = as.factor(year))) + 
  geom_point(position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = lcl.rand, ymax = ucl.rand), position=position_dodge(width=0.5)) + 
  scale_y_continuous(limits = c(0, 750)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90), strip.background = element_blank()) +
  facet_wrap(~Basin, scales = "free_x")

ggsave(filename = "ne_plot.pdf", nep,
       width = 11, height = 8.5)


ggsave(filename = "ne_plot_natural_scale.pdf", nepnl,
       width = 11, height = 8.5)


