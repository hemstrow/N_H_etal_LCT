library(ggplot2); library(plotly)

# read in data
sample_meta <- read.table("metadata.txt", header = T)
sample_meta$stream <- paste0("str_", sample_meta$stream)
sample_meta$stream[sample_meta$watershed == "LHU" & sample_meta$stream == "str_INC"] <- "str_INK"
pd <- read.table("bamlist1_omy_pca.pc_scores.txt")
pd <- merge(pd, sample_meta, by.x = "Lab", by.y = "ID")

# add hyb percent
hyb <- readr::read_csv("LCT_RBT_final hybridization 2-19-19.csv")
colnames(hyb)[1] <- "ID"
colnames(hyb)[5] <- "Relative.percent.RBT"
pd <- merge(pd, hyb, by.x = "Lab", by.y = "ID")

# plot
ggsave(filename = "updated_PCA_plot.pdf", device = "pdf", plot = ggplot(pd, aes(x = PC1, y = PC2, color = watershed, size = Relative.percent.RBT, label = watershed)) + 
  geom_label() + theme_bw() +
  scale_color_viridis_d(), width = 11, height = 8.5)

ggplotly(ggplot(pd, aes(x = PC1, y = PC2, color = watershed, size = Relative.percent.RBT)) + geom_point() + theme_bw() +
  scale_color_viridis_d())
