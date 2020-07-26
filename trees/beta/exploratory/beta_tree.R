######### BETA NJ/UPGMA TREE #########

# dependencies

library(readxl)
library(ggplot2)
library(ggtree)
library(phangorn)
library(ape)
library(cowplot)

# load beta msa  produced from "msa/beta/beta_msa.R"  
TDPbeta_aln <- readRDS("msa/beta/TDP1beta_msa.rds")

# compute similarity and identity matrix

# this point is crucial! some other packages that may be useful: "protr", "phanghorn"

matrix_id <- as.matrix(seqinr::dist.alignment(TDPbeta_aln, "identity", gap = 10))
matrix_sim <- as.matrix(seqinr::dist.alignment(TDPbeta_aln, "similarity", gap = 10)) # fitch margulis model is the only option

# however when you give this to mega it won't matter because
# it will use the best sub. model -> see excel

# extract metadate from data/TDP_metadata.xlsx

#build NJ tree from similarity Fitch matrix
nj_tree <- ape::nj(matrix_sim)

treenj <- ggtree::groupOTU(nj_tree, .node = c("M. truncatula"))

##### NJ TREE #####

# extract metadate from data/TDP_metadata.xlsx
tip_metadatabeta <- readxl::read_xlsx(path = "data/TDP_Dataset.xlsx", sheet = 5)

#build NJ tree from similarity Fitch matrix
nj_tree <- ape::nj(matrix_sim)

treenj <- groupOTU(nj_tree, .node = c("M. truncatula"))

nj_beta <- ggtree(treenj, aes(color = group)) +
  scale_color_manual(values = c("black", "firebrick")) +
  geom_tiplab(hjust = -0.06, align = T, linesize = .5) +
  geom_tippoint() +
  ggtitle("NJ Tree of TDP1beta") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  ggplot2::xlim(NA, 0.42)

nj_beta <- ggtree::rotate(nj_beta, 41)

tip_metadatanjbeta <- tip_metadatabeta[match(nj_beta$data$label, tip_metadatabeta$Name), ]

nj_beta <- nj_beta + 
  geom_tippoint(aes(fill = as.character(tip_metadatanjbeta$Color_lab), shape = as.factor(tip_metadatanjbeta$PCH_lab)), 
    col = "black", size = 4) + 
  scale_fill_manual(values = c("yellowgreen", "aquamarine3", "#880085")) + 
  scale_shape_manual(values = c(21, 21, 24))

ggsave(filename = "trees/beta/exploratory/NJ_beta.pdf", plot = nj_beta, 
  height = 13, width = 10, device = "pdf")

# write nwk file of nj tree
write.tree(as.phylo(nj_beta), file = "trees/beta/exploratory/nj_beta.nwk", append = FALSE, digits = 10, tree.names = FALSE)

# save tree plot
#saveRDS(nj_beta, file = "trees/both/nj_beta.rds")

###### UPGMA beta #####

upgma_treebeta <- phangorn::upgma(matrix_sim)

treeupgmabeta <- groupOTU(upgma_treebeta, .node = c("M. truncatula"))

upgma_beta <- ggtree(treeupgmabeta, aes(color = group)) +
  geom_tiplab(hjust = +1.2) +
  scale_color_manual(values = c("black", "firebrick")) +
  ggtitle("UPGMA Tree of TDP1beta") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  scale_x_reverse() +
  #geom_text(aes(label=node)) + 
  scale_y_reverse() +
  xlim(0.6, 0)
  
upgma_beta <- ggtree::rotate(upgma_beta, 38)
upgma_beta <- ggtree::rotate(upgma_beta, 40)
upgma_beta <- upgma_beta + geom_strip(3, 6, label = "Monocots", align = T, angle = 30, color = "aquamarine3",
    offset = 0.29, barsize = 1, offset.text =  -0.002, fontsize = 4) +
  geom_strip(1, 13, label = "Ancestral Plants", align = T, angle = 30, color = "#880085",
    offset = 0.29, barsize = 1, extend = 0.15, offset.text = -0.002, fontsize = 4) +
  geom_strip(29, 33, label = "Eudicots", align = T, angle = 30, color = "yellowgreen", 
    offset = 0.29, barsize = 1, offset.text = -0.002, fontsize = 4)

tip_metadataupgmabeta <- tip_metadatabeta[match(upgma_beta$data$label, tip_metadatabeta$Name), ]

upgma_beta <- upgma_beta + 
  geom_tippoint(aes(fill = as.factor(tip_metadataupgmabeta$Color_lab), shape = as.factor(tip_metadataupgmabeta$PCH_lab)), col = "black", size = 4) +
  scale_fill_manual(values = c("yellowgreen", "aquamarine3", "#880085")) +
  scale_shape_manual(values = c(21, 21, 24))   

ggsave(filename = "trees/beta/exploratory/upgma_beta.pdf", plot = upgma_beta, 
  height = 13, width = 10, device = "pdf")

# write nwk file of nj tree
write.tree(as.phylo(upgma_beta), file = "trees/beta/exploratory/upgma_beta.nwk", append = FALSE, digits = 10, tree.names = FALSE)

# save tree plot
#saveRDS(upgma_beta, file = "trees/both/upgma_beta.rds")

grid <- cowplot::plot_grid(upgma_beta, nj_beta,  labels = "AUTO" )

ggsave(filename = "trees/beta/exploratory/beta_trees.pdf", plot = grid, 
  width = 20, height = 13, device = "pdf")