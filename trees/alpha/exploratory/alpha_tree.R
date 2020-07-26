######### ALPHA NJ/UPGMA TREE #########

# dependencies
library(readxl)
library(ggplot2)
library(ggtree)
library(phangorn)
library(cowplot)
library(ape)

# load alpha msa produced from "msa/alpha/alpha_msa.R"  

TDPalpha_aln <- readRDS("msa/alpha/TDP1alpha_msa.rds")

# compute similarity and identity matrix

# this point is crucial! some other packages that may be useful: "protr", "phanghorn"

matrix_id <- as.matrix(seqinr::dist.alignment(TDPalpha_aln, "identity", gap = 10))
matrix_sim <- as.matrix(seqinr::dist.alignment(TDPalpha_aln, "similarity", gap = 10)) # fitch margulis model is the only option

# however when you give this to mega it won't matter because
# it will use the best sub. model after optimal model search -> see excel

##### NJ ALPHA #####

# extract metadate from data/TDP_metadata.xlsx
tip_metadataalpha <- readxl::read_xlsx(path = "data/TDP_Dataset.xlsx", sheet = 2)
#build NJ tree from similarity Fitch matrix
nj_tree <- ape::nj(matrix_sim)

treenj <- ggtree::groupOTU(nj_tree, .node = c("M. truncatula"))

njtree_alpha <- ggtree(treenj, aes(color = group)) +
  scale_color_manual(values = c("black", "firebrick")) +
  geom_tiplab(hjust = -0.07, align = T, linesize = .5) +
  ggtitle("NJ Tree of TDP1alpha") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  ggplot2::xlim(NA, 0.77) +
  geom_strip(2, 7, label = "Animals", align = T, angle = 30, color = "slateblue4",
    offset = 0.193, offset.text = 0.007, fontsize = 4, barsize = 1) +
  geom_cladelabel(node = 1, label = "Fungi", align = T, angle = 30, color = "brown", 
    offset = 0.193, offset.text = 0.007, fontsize = 4, barsize = 1, extend = 0.5) +
  geom_strip(88, 85, label = "Algae", align = T, angle = 30, color = "#ff6700", 
    offset = 0.193, offset.text = 0.007, fontsize = 4, barsize = 1) +
  geom_strip(64, 80, label = "Ancestral Plants\nand Gymnosperms", align = T, angle = 30, 
    color = "#880085", offset = 0.193, offset.text = 0.01, fontsize = 4, barsize = 1) +
  geom_strip(72, 83, label = "Monocots", align = T, angle = 30, color = "aquamarine3", 
    offset = 0.193, offset.text = 0.007, fontsize = 4, barsize = 1) +
  geom_strip(37, 55, label = "Eudicots", align = T, angle = 30, color = "yellowgreen", 
    offset = 0.193, offset.text = 0.007, fontsize = 4, barsize = 1)

njtree_alpha <- flip(njtree_alpha, 166, 1)

tip_metadatanj <- tip_metadataalpha[match(njtree_alpha$data$label, tip_metadataalpha$Name), ]

# add colored tips
njtree_alpha <- njtree_alpha + 
  geom_tippoint(aes(fill = as.character(tip_metadatanj$Color_lab), 
    shape = as.factor(tip_metadatanj$PCH_lab)), col = "black", size = 3.5) + 
  scale_fill_manual(values = c("slateblue4", "yellowgreen", "#880085", "brown", "#ff6700", "aquamarine3")) + 
  scale_shape_manual(values = c(22, 21, 21, 23, 25, 24))

#gridnjalpha <- plot_grid(njtree_alpha, labels = "AUTO")

# TODO fix fungi label outside plot

ggsave(filename = "trees/alpha/exploratory/NJ_alpha.pdf", plot = njtree_alpha, 
       height = 13, width = 10, device = "pdf")

# write nwk file of nj tree
write.tree(as.phylo(njtree_alpha), file = "trees/alpha/exploratory/nj_alpha.nwk", append = FALSE, digits = 10, tree.names = FALSE)

# save tree plot
#saveRDS(njtree_alpha, file = "trees/both/nj_alpha.rds")

##### UPGMA ALPHA ######

upgma_tree <- phangorn::upgma(matrix_sim)

treeupgma <- groupOTU(upgma_tree, .node = c("M. truncatula"))

upgma_alpha <- ggtree(treeupgma, aes(color = group)) +
  scale_color_manual(values = c("black", "firebrick")) +
  geom_tiplab(hjust = -0.09) +
  ggtitle("UPGMA Tree of TDP1alpha") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  ggplot2::xlim(NA, 0.5)

tip_metadataupgma <- tip_metadataalpha[match(upgma_alpha$data$label, tip_metadataalpha$Name), ]

upgma_alpha <- upgma_alpha + geom_tippoint(aes(fill = as.character(tip_metadataupgma$Color_lab), 
  shape = as.factor(tip_metadataupgma$PCH_lab)), col = "black", size = 4) + 
  scale_fill_manual(values = c("slateblue4","yellowgreen","#880085","brown","#ff6700","aquamarine3")) +
  scale_shape_manual(values = c(22, 21, 21, 23, 25, 24)) 

# uglier and without annotations because we decided to onlyh keep nj

ggsave(filename = "trees/alpha/exploratory/UPGMA_alpha.pdf", plot = upgma_alpha, 
  height = 13, width = 10, device = "pdf")

# write nwk file of upgma
write.tree(as.phylo(upgma_alpha), file = "trees/alpha/exploratory/upgma_alpha.nwk", append = FALSE, digits = 10, tree.names = FALSE)

# save tree plot
#saveRDS(upgma_alpha, file = "trees/both/upgma_alpha.rds")

grid <- cowplot::plot_grid(upgma_alpha, njtree_alpha,  labels = "AUTO" )
ggsave(filename = "trees/alpha/exploratory/alpha_trees.pdf", plot = grid, 
  width = 20, height = 13, device = "pdf")