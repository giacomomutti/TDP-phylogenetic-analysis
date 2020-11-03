######### ALPHA NJ/UPGMA TREE #########

# dependencies
library(readxl)
library(ggplot2)
library(ggtree)
library(phangorn)
library(cowplot)
library(ape)

# load alpha msa produced from "msa/alpha/alpha_msa.R"  

TDPboth_aln <- readRDS("msa/both/TDP1both_msa.rds")

# compute similarity and identity matrix

# this point is crucial! some other packages that may be useful: "protr", "phanghorn"

matrix_id <- as.matrix(seqinr::dist.alignment(TDPboth_aln, "identity", gap = 10))
matrix_sim <- as.matrix(seqinr::dist.alignment(TDPboth_aln, "similarity", gap = 10)) # fitch margulis model is the only option

# however when you give this to mega it won't matter because
# it will use the best sub. model after optimal model search -> see excel

##### NJ ALPHA #####

# extract metadate from data/TDP_metadata.xlsx
#build NJ tree from similarity Fitch matrix
nj_tree <- ape::nj(matrix_sim)

treenj <- ggtree::groupOTU(nj_tree, .node = c("M. truncatula"))

njtree_both <- ggtree(treenj, aes(color = group)) +
  scale_color_manual(values = c("black", "firebrick")) +
  geom_tiplab(hjust = -0.07, align = T, linesize = .5) +
  ggtitle("NJ Tree of TDP1alpha and beta") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  ggplot2::xlim(NA, 0.77) 

ggsave(filename = "trees/both/exploratory/NJ_both.pdf", plot = njtree_both, 
       height = 13, width = 10, device = "pdf")


# write nwk file of nj tree
write.tree(as.phylo(njtree_both), file = "trees/both/exploratory/nj_both.nwk", append = FALSE, digits = 10, tree.names = FALSE)

# save tree plot
#saveRDS(njtree_alpha, file = "trees/both/nj_alpha.rds")

##### UPGMA ALPHA ######

upgma_tree <- phangorn::upgma(matrix_sim)

treeupgma <- groupOTU(upgma_tree, .node = c("M. truncatula"))

upgma_both <- ggtree(treeupgma, aes(color = group)) +
  scale_color_manual(values = c("black", "firebrick")) +
  geom_tiplab(hjust = -0.09) +
  ggtitle("UPGMA Tree of TDP1alpha and beta") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  ggplot2::xlim(NA, 0.5)

# uglier and without annotations because we decided to onlyh keep nj

ggsave(filename = "trees/both/exploratory/UPGMA_both.pdf", plot = upgma_both, 
       height = 13, width = 10, device = "pdf")

# write nwk file of upgma
write.tree(as.phylo(upgma_alpha), file = "trees/both/exploratory/upgma_both.nwk", append = FALSE, digits = 10, tree.names = FALSE)

# save tree plot
#saveRDS(upgma_alpha, file = "trees/both/upgma_alpha.rds")

grid <- cowplot::plot_grid(upgma_alpha, njtree_alpha,  labels = "AUTO" )
ggsave(filename = "trees/both/exploratory/both_trees.pdf", plot = grid, 
       width = 20, height = 13, device = "pdf")
