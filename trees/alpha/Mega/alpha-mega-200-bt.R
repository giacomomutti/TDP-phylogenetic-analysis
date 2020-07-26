# dependencies

library("ape")
library("Biostrings")
library("ggplot2")
library("ggtree")


# this nwk tree has been produced by MEGAX with:
# input = msa/alpha/TDP1alpha_msa.fasta
# Model = ML
# No bootstrap = 200
# Sub model = JTT + G -> found w/ best model selection function by MEGA

# read tree
treealpha <- read.tree(file = "trees/alpha/Mega/mega-alpha-200bt.nwk")
treealpha$tip.label <- gsub("_", " ", treealpha$tip.label)

# highlight truncatula node
a <- groupOTU(treealpha, .node = c("M. truncatula"))

# tree visualization
plotalpha <- ggtree(a, layout = "circular", branch.length = "none", aes(color = group)) +
  ggtitle("ML 200 bootstrap TDP1alpha tree ") +
  geom_text(aes(label = node)) +
  theme(plot.title = element_text(hjust = 0.5, size = 30)) +
  scale_color_manual(values = c("black", "firebrick")) +
  geom_tiplab2(size = 8, align = T)

# rotate, flip and add bootstrap values
bt_list <- round(as.numeric(plotalpha$data$label) * 100)
plotalpha <- rotate(plotalpha, node = 117)
plotalpha <- flip(plotalpha, 116, 102)
plotalpha <- flip(plotalpha, 97, 140)
plotalpha <- flip(plotalpha, 166, 165)
plotalpha <- flip(plotalpha, 152, 62)

# add toaxonomic strips
plotalpha <- plotalpha +
  geom_label2(aes(label = bt_list), color = "black", size = 3, fill = "lightblue", fontface = "bold") +
  geom_cladelabel(
    node = 88, label = "Fungi", align = T, angle = 30, color = "brown",
    offset = 9.83, offset.text = 0.4, fontsize = 10, barsize = 1.5, extend = 0.5
  ) +
  geom_strip(
      82, 87, label = "Animals", align = T, angle = 30, color = "slateblue4", 
      offset = 9.83, offset.text = 0.4, fontsize = 10, barsize = 1.5) +
  geom_strip(
      78, 81, label = "Algae", align = T, angle = 30, color = "#ff6700", 
      offset = 9.83, offset.text = 0.4, fontsize = 10, barsize = 1.5) +
  geom_strip(
      73, 74, label = "Ancestral\nPlants", align = T, angle = 30, color = "#880085", 
      offset = 9.83, offset.text = 0.78, fontsize = 10, barsize = 1.5) +
  geom_strip(
      77, 76, label = "Gymnosperms", align = T, angle = 30, color = "violet", 
      offset = 9.83, offset.text = 0.4, fontsize = 10, barsize = 1.5) +
  geom_strip(
      61, 59, label = "Liliod\nmonocots", align = T, angle = 30, color = "aquamarine3", 
      offset = 9.83, offset.text = 1.2, fontsize = 10, barsize = 1.5) +
  geom_strip(
      72, 62, label = "Commelinids\nmonocots", align = T, angle = 30, color = "deepskyblue4", 
      offset = 9.83, offset.text = 1.7, fontsize = 10, barsize = 1.5) +
  geom_cladelabel(
      node = 57, label = "Basal\neudicots", align = T, angle = 38, color = "darkgreen", 
      offset = 9.83, offset.text = 2, fontsize = 10, barsize = 1.5, extend = 0.5) +
  geom_strip(
      56, 18, label = "Rosids", align = T, angle = 30, color = "yellowgreen", 
      offset = 9.83, offset.text = 4.4, fontsize = 10, barsize = 1.5) +
  geom_strip(
      17, 55, label = "Superasterids", align = T, angle = 30, color = "#4a6200", 
      offset = 9.83, offset.text = 0.4, fontsize = 10, barsize = 1.5)

# read metadata to assign correct tip shape and color
metadata_tip <- readxl::read_xlsx(path = "data/TDP_Dataset.xlsx", sheet = 3)
metadata_200bt <- metadata_tip[match(plotalpha$data$label, metadata_tip$Name), ]
dfxl <- as.data.frame(metadata_200bt)

def_plotalpha <- plotalpha +
  geom_tippoint(aes(shape = as.factor(dfxl$PCH_lab), fill = as.character(dfxl$Color_lab)), col = "black", size = 7.2) +
  scale_fill_manual(values = c(
    "#ff6700", "aquamarine3", "#880085", "violet",
    "slateblue4", "#4a6200", "darkgreen", "yellowgreen", "brown", "deepskyblue4"
  )) +
  scale_shape_manual(values = c(22, 21, 21, 23, 25, 24)) +
  theme(legend.position = "none")

# save
ggsave(
  filename = "trees/alpha/Mega/ML-200bt-alpha.pdf", plot = def_plotalpha, device = "pdf",
  height = 25, width = 25)