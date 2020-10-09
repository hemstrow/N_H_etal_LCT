library(ggplot2); library(plotly)

# read in data
# sample_meta <- read.table("metadata.txt", header = T)
# sample_meta$stream <- paste0("str_", sample_meta$stream)
# sample_meta$stream[sample_meta$watershed == "LHU" & sample_meta$stream == "str_INC"] <- "str_INK"
pd <- read.table("bamlist1_omy_pca.pc_scores.txt", stringsAsFactors = F)
# pd <- merge(pd, sample_meta, by.x = "Lab", by.y = "ID")

# add hyb percent
hyb <- readr::read_csv("LCT_RBT_final hybridization 2-19-19.csv")
colnames(hyb)[1] <- "ID"
colnames(hyb)[5] <- "Relative.percent.RBT"

hyb.fix <- hyb$ID[which(!hyb$ID %in% pd$Lab)]
pd.fix <- pd$Lab[which(!pd$Lab %in% hyb$ID)]

hyb.fix <- data.frame(hyb.lab = hyb.fix, index = substr(hyb.fix, 9, 16), stringsAsFactors = F)
pd.fix <- data.frame(pd.lab = pd.fix, index = substr(pd.fix, 9, 16), stringsAsFactors = F)
fix.lab <- merge(hyb.fix, pd.fix, "index", all = T)
fix.lab <- fix.lab[-1,]

matches <- match(fix.lab$pd.lab, pd$Lab)

pd$Lab[matches] <- fix.lab$hyb.lab
pd <- merge(hyb, pd, by.x = "ID", by.y = "Lab")
pd$watershed <- gsub("_.+_.+_.+", "", pd$ID)
pd$samp <- substr(pd$ID, 1, 7)

# plot
ggsave(filename = "updated_PCA_plot.pdf", device = "pdf", plot = ggplot(pd, aes(x = PC1, y = PC2, color = watershed, size = Relative.percent.RBT, label = samp)) + 
  geom_label() + theme_bw() +
  scale_color_viridis_d(), width = 11, height = 8.5)

ggplotly(ggplot(pd, aes(x = PC1, y = PC2, color = watershed, size = Relative.percent.RBT)) + geom_point() + theme_bw() +
  scale_color_viridis_d())
