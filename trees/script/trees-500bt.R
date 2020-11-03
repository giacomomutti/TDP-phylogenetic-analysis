# dependencies

library("ape")
library("Biostrings")
library("ggplot2")
library("ggtree")
library("glue")

#####TDP1 ALPHA######

# this nwk tree has been produced by MEGAX with:
# input = msa/alpha/TDP1alpha_msa.fasta
# Model = ML
# No bootstrap = 500
# Sub model = JTT + G -> found w/ best model selection function by MEGA

# read tree

treealpha <- read.tree(file = "trees/data/alpha_500bt.nwk")

treealpha$tip.label <- gsub("_", " ", treealpha$tip.label)
treealpha$tip.label <- gsub(" subsp. vesca", " ", treealpha$tip.label)
treealpha$tip.label <- gsub(" Japonica Group", " ", treealpha$tip.label)
treealpha$tip.label <- gsub(" var. chinensis", " ", treealpha$tip.label)
treealpha$tip.label <- gsub(" subsp. lyrata", " ", treealpha$tip.label)
treealpha$tip.label <- gsub(" subsp. pepo", " ", treealpha$tip.label)
treealpha$tip.label <- gsub(" var. radiata", " ", treealpha$tip.label)
treealpha$tip.label <- gsub(" var. scolymus", " ", treealpha$tip.label)
treealpha$tip.label <- gsub(" subsp. lyrata", " ", treealpha$tip.label)
treealpha$tip.label <- gsub(" subsp. tauschii", " ", treealpha$tip.label)
treealpha$tip.label <- gsub(" C-169", " ", treealpha$tip.label)
treealpha$tip.label <- trimws(treealpha$tip.label, "right")

# highlight truncatula node
a <- groupOTU(treealpha, .node = c("M. truncatula"))

# tree visualization
plotalpha <- ggtree(a, layout = "circular", branch.length = "none", aes(color = group)) +
  # geom_text(aes(label = node)) +
  geom_label2(aes(subset = (as.numeric(label) * 100 < 99.5), label = round(as.numeric(label) * 100, 0)),
              color = "black",
              size = 2, label.size = 0.1,
              fill = "lightblue", fontface = "bold"
  ) +
  scale_color_manual(values = c("black", "firebrick")) +
  geom_tiplab2(size = 3.5, align = T, hjust = -0.09, fontface = 'italic') +
  theme(legend.position = "none")
# rotate, flip and add bootstrap values
plotalpha <- rotate(plotalpha, node = 164)
plotalpha <- flip(plotalpha, 97, 140)
plotalpha <- flip(plotalpha, 116, 102)
plotalpha <- flip(plotalpha, 118, 17)
# plotalpha <- flip(plotalpha, 152, 62)


off <- 5.7
off_text <- 0.05
barsize <- 1
fontsize <- 5

# add taxonomic strips
plotalpha <- plotalpha +
  geom_cladelabel(
    node = 88, "Fungi", align = T, angle = 30,
    color = "brown", offset = off, offset.text = off_text + 0.01,
    fontsize = fontsize, barsize = barsize, extend = 0.2
  ) +
  geom_strip(82, 87, "Animals",
             align = T, angle = 30,
             color = "slateblue4", offset = off, offset.text = off_text,
             fontsize = fontsize, barsize = barsize
  ) +
  geom_strip(78, 81, "Algae",
             align = T, angle = 30,
             color = "#ff6700", offset = off, offset.text = off_text + 0.09,
             fontsize = fontsize, barsize = barsize
  ) +
  geom_strip(73, 75, "Ancestral\nPlants",
             align = T, angle = 30,
             color = "#880085", offset = off, offset.text = off_text + 0.6,
             fontsize = fontsize, barsize = barsize
  ) +
  geom_strip(77, 76, "Gymnosperms",
             align = T, angle = 30,
             color = "violet", offset = off, offset.text = off_text + 0.06,
             fontsize = fontsize, barsize = barsize
  ) +
  geom_strip(62, 72, "Monocots",
             align = T, angle = 30,
             color = "deepskyblue4", offset = off, offset.text = off_text + 0.34,
             fontsize = fontsize, barsize = barsize
  ) +
  geom_cladelabel(57, "Basal\neudicots",
                  align = T, angle = 38,
                  color = "darkgreen", offset = off, offset.text = off_text + 1,
                  fontsize = fontsize, barsize = barsize, extend = 0.4
  ) +
  geom_strip(56, 25, "Rosids",
             align = T, angle = 30,
             color = "yellowgreen", offset = off, offset.text = off_text + 3,
             fontsize = fontsize, barsize = barsize
  ) +
  geom_strip(17, 55, "Superasterids",
             align = T, angle = 30,
             color = "#4a6200", offset = off, offset.text = off_text + 0.4,
             fontsize = fontsize, barsize = barsize
  )

# read metadata to assign correct tip shape and color
metadata_tip <- readxl::read_xlsx(path = "data/TDP_Dataset.xlsx", sheet = 3)
metadata_500bt <- metadata_tip[match(plotalpha$data$label, metadata_tip$Name), ]
dfxl <- as.data.frame(metadata_500bt)

plotalpha <- plotalpha +
  geom_tippoint(aes(
    shape = as.factor(dfxl$PCH_lab),
    fill = as.character(dfxl$Color_lab)
  ), col = "black", size = 3.5) +
  scale_fill_manual(values = c(
    "#ff6700", "#880085", "violet",
    "slateblue4", "#4a6200", "darkgreen",
    "yellowgreen", "brown", "deepskyblue4"
  )) +
  scale_shape_manual(values = c(22, 21, 21, 23, 25, 24)) +
  theme(legend.position = "none")

ggsave(
  filename = "trees/plots/ML-500bt-alpha.pdf",
  plot = plotalpha, device = cairo_pdf, height = 15, width = 15
)
ggsave(
  filename = "trees/plots/ML-500bt-alpha.png",
  plot = plotalpha, device = "png",
  height = 15, width = 15, dpi = 350
)

### Tree with MSA plot

treealpha <- read.tree("trees/data/alpha_500bt.nwk")
treealpha$tip.label <- gsub("_", " ", treealpha$tip.label)

p1 <- ggtree(treealpha) + geom_tiplab(size = 3, fontface = 'italic')
p1 <- msaplot(p1, "msa/alpha/TDP1alpha_msa_rev.fasta", offset = 5.8, width = 1.8)

off <- 2.5
off_text <- 0.03
barsize <- 0.5
fontsize <- 3

# highlight truncatule node
p1 <- p1 + theme(legend.position = "none") +
  geom_cladelabel(
    node = 88, "Fungi", align = T, color = "brown",
    offset = off, offset.text = off_text,
    fontsize = fontsize, barsize = barsize, extend = 0.6
  ) +
  geom_strip(82, 87, "Animals",
             align = T, color = "slateblue4",
             offset = off, offset.text = off_text,
             fontsize = fontsize, barsize = barsize
  ) +
  geom_strip(78, 81, "Algae",
             align = T, color = "#ff6700",
             offset = off, offset.text = off_text,
             fontsize = fontsize, barsize = barsize
  ) +
  geom_strip(74, 75, "Ancestral Plants",
             align = T, color = "#880085",
             offset = off, offset.text = off_text,
             fontsize = fontsize, barsize = barsize
  ) +
  geom_cladelabel(
    node = 73, label = "", align = T, color = "#880085",
    offset = off, offset.text = off_text,
    fontsize = fontsize, barsize = barsize, extend = 0.6
  ) +
  geom_strip(77, 76, "Gymnosperms",
             align = T, color = "violet",
             offset = off, offset.text = off_text,
             fontsize = fontsize, barsize = barsize
  ) +
  geom_strip(62, 72, "Monocots",
             align = T, color = "deepskyblue4",
             offset = off, offset.text = off_text,
             fontsize = fontsize, barsize = barsize
  ) +
  geom_cladelabel(57, "Basal Eudicots",
                  align = T, color = "darkgreen",
                  offset = off, offset.text = off_text,
                  fontsize = fontsize, barsize = barsize, extend = 0.6
  ) +
  geom_strip(40, 2, "Rosids",
             align = T, color = "yellowgreen",
             offset = off, offset.text = off_text,
             fontsize = fontsize, barsize = barsize
  ) +
  geom_strip(45, 55, "Superasterids",
             align = T, color = "#4a6200",
             offset = off, offset.text = off_text,
             fontsize = fontsize, barsize = barsize
  ) +
  geom_cladelabel(17,
                  label = "", align = T, color = "#4a6200",
                  offset = off, offset.text = off_text,
                  fontsize = fontsize, barsize = barsize, extend = 0.6
  ) +
  geom_cladelabel(56,
                  label = "", align = T, color = "yellowgreen",
                  offset = off, offset.text = off_text,
                  fontsize = fontsize, barsize = barsize, extend = 0.6
  )

# save
ggsave(
  filename = "trees/plots/ML-500bt-alpha_msa.pdf",
  plot = p1, device = cairo_pdf, height = 8, width = 10
)

ggsave(
  filename = "trees/plots/ML-500bt-alpha_msa.png",
  plot = p1, device = "png", height = 8, width = 10, dpi = 350
)

##########TDP1 BETA#############

# this nwk tree has been produced by MEGAX with:
# input = msa/beta/TDP1beta_msa.fasta
# Model = ML
# No bootstrap = 500
# Sub model = JTT + G + I + F -> found w/ best model selection function by MEGA

# read beta tree
treebeta <- read.tree(file = "trees/data/beta_500bt.nwk")

treebeta$tip.label <- gsub("_", " ", treebeta$tip.label)
treebeta$tip.label <- gsub(" subsp. vesca", " ", treebeta$tip.label)
treebeta$tip.label <- gsub(" Japonica Group", " ", treebeta$tip.label)
treebeta$tip.label <- gsub(" var. chinensis", " ", treebeta$tip.label)
treebeta$tip.label <- gsub(" subsp. lyrata", " ", treebeta$tip.label)
treebeta$tip.label <- trimws(treebeta$tip.label, "right")

# highlight truncatula node
mt <- groupOTU(treebeta, .node = c("M. truncatula"))

# tree visualization
plotbeta <- ggtree(mt,
  layout = "circular", branch.length = "none",
  aes(color = group)
) +
  # geom_text(aes(label=node)) +
  geom_label2(aes(subset = (as.numeric(label) * 100 < 99.5), label = round(as.numeric(label) * 100, 0)),
    color = "black",
    size = 2, label.size = 0.1,
    fill = "lightblue",
    fontface = "bold"
  ) +
  scale_color_manual(values = c("black", "firebrick")) +
  geom_tiplab2(size = 3, align = T, hjust = -0.09, fontface = 'italic') +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.position = "none"
  )

# flip and add bootstrap values
plotbeta <- flip(plotbeta, 41, 64)

off <- 5.8
off_text <- 0.2
fontsize <- 4
barsize <- 1

# add taxonomic strips
plotbeta <- plotbeta +
  geom_strip(35, 36,
    label = "Ancestral\nPlants", align = T, angle = 30,
    color = "#880085", offset = off, offset.text = off_text,
    fontsize = fontsize, barsize = barsize
  ) +
  geom_strip(29, 34,
    label = "Monocots", align = T, angle = 30,
    color = "deepskyblue4", offset = off, offset.text = off_text,
    fontsize = fontsize, barsize = barsize
  ) +
  geom_strip(28, 2,
    label = "Rosids", align = T, angle = 30,
    color = "yellowgreen", offset = off, offset.text = off_text,
    fontsize = fontsize, barsize = barsize, hjust = 1
  ) +
  geom_strip(25, 27,
    label = "Asterids", align = T, angle = 30,
    color = "#4a6200", offset = off, offset.text = off_text,
    fontsize = fontsize, barsize = barsize
  )

# read metadata to assign correct tip shape and color
metadata_tip_beta <- readxl::read_xlsx(path = "data/TDP_Dataset.xlsx", sheet = 6)
metadata_500bt_beta <- metadata_tip_beta[match(plotbeta$data$label, metadata_tip_beta$Name), ]
dfxl_beta <- as.data.frame(metadata_500bt_beta)

plotbeta <- plotbeta +
  geom_tippoint(aes(
    shape = as.factor(dfxl_beta$PCH_lab),
    fill = as.character(dfxl_beta$Color_lab)
  ),
  col = "black", size = 3
  ) +
  scale_fill_manual(values = c("yellowgreen", "#4a6200", "deepskyblue4", "#880085")) +
  scale_shape_manual(values = c(21, 21, 24)) +
  theme(legend.position = "none")

# save
ggsave(
  filename = "trees/plots/ML-500bt-beta.pdf", plot = plotbeta,
  device = cairo_pdf,
  height = 8, width = 8
)

ggsave(
  filename = "trees/plots/ML-500bt-beta.png", plot = plotbeta,
  device = "png", dpi = 350,
  height = 8, width = 8
)

### Tree with MSA plot

treebeta <- read.tree(file = "trees/data/beta_500bt.nwk")
treebeta$tip.label <- gsub("_", " ", treebeta$tip.label)

p <- ggtree(treebeta) + geom_tiplab(size = 3, fontface = 'italic') # geom_treescale()
p <- msaplot(p, "msa/beta/TDP1beta_msa_rev.fasta", offset = 1.7, width = 2)

# add taxonomic bars

p <- p + theme(legend.position = "none") +
  geom_strip(35, 36,
             label = "Ancestral\nplants", align = T, color = "#880085",
             offset = 0.9, offset.text = 0.03, fontsize = 3, barsize = 0.5
  ) +
  geom_strip(29, 34,
             label = "Monocots", align = T, color = "deepskyblue4",
             offset = 0.9, offset.text = 0.03, fontsize = 3, barsize = 0.5
  ) +
  geom_strip(28, 2,
             label = "Rosids", align = T, color = "yellowgreen",
             offset = 0.9, offset.text = 0.03, fontsize = 3, barsize = 0.5
  ) +
  geom_strip(25, 27,
             label = "Asterids", align = T, color = "#4a6200",
             offset = 0.9, offset.text = 0.03, fontsize = 3, barsize = 0.5
  )

# save
ggsave(filename = "trees/plots/ML-500bt-beta_msa.pdf", plot = p, device = cairo_pdf, height = 8, width = 10)
ggsave(filename = "trees/plots/ML-500bt-beta_msa.png", plot = p, device = "png", height = 8, width = 10, dpi = 350)