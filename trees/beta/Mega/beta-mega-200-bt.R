# dependencies

library("ape")
library("Biostrings")
library("ggplot2")
library("ggtree")

# this nwk tree has been produced by MEGAX with:
# input = msa/beta/TDP1beta_msa.fasta
# Model = ML
# No bootstrap = 400
# Sub model = JTT + G + I + F -> found w/ best model selection function by MEGA

# read tree
treebeta <- read.tree(file = "trees/beta/Mega/mega-beta-200bt.nwk")
treebeta$tip.label <- gsub("_", " ",treebeta$tip.label)

# highlight truncatule node
b <- groupOTU(treebeta, .node = c("M. truncatula"))

# tree visualization
plotbeta <- ggtree(b, layout = "circular", branch.length = "none",aes(color=group)) +
    ggtitle("ML 400 bootstrap TDP1beta tree") +
    #geom_text(aes(label=node)) +
    theme(plot.title = element_text(hjust = 0.5, size = 30))  +
    scale_color_manual(values = c("black", "firebrick")) +
    geom_tiplab2(size=8, align = T, hjust = -0.04)

# flip and add bootstrap values
plotbeta <- flip(plotbeta, 51,43)
bt_list <-  round(as.numeric(plotbeta$data$label)*100)

# add taxonomic strips
plotbeta <- plotbeta +
  geom_label2(aes(label = bt_list),color="black", size = 3, fill="lightblue", fontface="bold") +
  geom_strip(
    35,36, label = "Ancestral\nPlants",align = T, angle = 30, color = "#880085", 
    offset = 12, offset.text = 0.7, fontsize = 10, barsize = 1.5) + 
  geom_strip(
    29,34, label = "Commelinids\nmonocots" ,align = T, angle = 30, color = "deepskyblue4",
    offset = 12, offset.text = 1.2, fontsize = 10, barsize = 1.5) + 
  geom_strip(
    28,2, label = "Rosids", align = T, angle = 30, color = "yellowgreen", 
    offset = 12, offset.text = 6.3, fontsize = 10, barsize = 1.5) +
  geom_strip(
    10,12, label = "Asterids", align = T, angle = 30, color = "#4a6200",
    offset = 12, offset.text = 0.7, fontsize = 10, barsize = 1.5) 

# read metadata to assign correct tip shape and color
metadata_tip_beta <- readxl::read_xlsx(path = "data/TDP_Dataset.xlsx", sheet = 6)
metadata_200bt_beta <- metadata_tip_beta[match(plotbeta$data$label, metadata_tip_beta$Name),]

dfxl_beta <- as.data.frame(metadata_200bt_beta)

def_plotbeta <- plotbeta + 
  geom_tippoint(aes(shape = as.factor(dfxl_beta$PCH_lab), fill = as.character(dfxl_beta$Color_lab)), col = "black", size = 7.2) + 
  scale_fill_manual(values = c("yellowgreen","#4a6200","deepskyblue4","#880085")) + 
  scale_shape_manual(values = c(21,21,24)) +
  theme(legend.position = "none")

#save
ggsave(filename = "trees/beta/Mega/ML-200bt-beta.pdf", plot = def_plotbeta, device = "pdf",
  height = 15, width = 15)