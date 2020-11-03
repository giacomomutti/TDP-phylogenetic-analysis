######### ALPHA NJ/UPGMA TREE #########

# dependencies
library(readxl)
library(ggplot2)
library(ggtree)
library(phangorn)
library(cowplot)
library(ape)

# load alpha msa produced from "msa/alpha/alpha_msa.R"

TDPalpha_aln <- readRDS("msa/alpha/TDP1alpha_msa_rev.rds")

# compute similarity and identity matrix

# this point is crucial! some other packages that may be useful: "protr", "phanghorn"

matrix_id <- as.matrix(seqinr::dist.alignment(TDPalpha_aln, "identity", gap = 10))
matrix_sim <- as.matrix(seqinr::dist.alignment(TDPalpha_aln, "similarity", gap = 10)) # fitch margulis model is the only option

# however when you give this to mega it won't matter because
# it will use the best sub. model after optimal model search -> see excel

##### NJ ALPHA #####

# extract metadate from data/TDP_metadata.xlsx
tip_metadataalpha <- readxl::read_xlsx(path = "data/TDP_Dataset.xlsx", sheet = 2)
# build NJ tree from similarity Fitch matrix
nj_tree <- ape::nj(matrix_sim)

treenj <- ggtree::groupOTU(nj_tree, .node = c("M. truncatula"))

njtree_alpha <- ggtree(treenj, aes(color = group)) +
  scale_color_manual(values = c("black", "firebrick")) +
  geom_tiplab(hjust = -0.07, align = T, linesize = .5, fontface = "italic") +
  ggtitle("NJ Tree of TDP1alpha") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  ggplot2::xlim(NA, 0.77) +
  geom_strip(2, 7,
    label = "Animals", align = T, angle = 30, color = "slateblue4",
    offset = 0.193, offset.text = 0.007, fontsize = 4, barsize = 1
  ) +
  geom_cladelabel(
    node = 1, label = "Fungi", align = T, angle = 30, color = "brown",
    offset = 0.193, offset.text = 0.007, fontsize = 4, barsize = 1, extend = 0.5
  ) +
  geom_strip(88, 85,
    label = "Algae", align = T, angle = 30, color = "#ff6700",
    offset = 0.193, offset.text = 0.007, fontsize = 4, barsize = 1
  ) +
  geom_strip(64, 80,
    label = "Ancestral Plants\nand Gymnosperms", align = T, angle = 30,
    color = "#880085", offset = 0.193, offset.text = 0.01, fontsize = 4, barsize = 1
  ) +
  geom_strip(72, 83,
    label = "Monocots", align = T, angle = 30, color = "aquamarine3",
    offset = 0.193, offset.text = 0.007, fontsize = 4, barsize = 1
  ) +
  geom_strip(37, 55,
    label = "Eudicots", align = T, angle = 30, color = "yellowgreen",
    offset = 0.193, offset.text = 0.007, fontsize = 4, barsize = 1
  ) +
  ggplot2::ylim(0, 88.5)

njtree_alpha <- flip(njtree_alpha, 166, 1)

tip_metadatanj <- tip_metadataalpha[match(njtree_alpha$data$label, tip_metadataalpha$Name), ]

# add colored tips
njtree_alpha <- njtree_alpha +
  geom_tippoint(aes(
    fill = as.character(tip_metadatanj$Color_lab),
    shape = as.factor(tip_metadatanj$PCH_lab)
  ), col = "black", size = 3.5) +
  scale_fill_manual(values = c("slateblue4", "yellowgreen", "#880085", "brown", "#ff6700", "aquamarine3")) +
  scale_shape_manual(values = c(22, 21, 21, 23, 25, 24))

# gridnjalpha <- plot_grid(njtree_alpha, labels = "AUTO")

# TODO fix fungi label outside plot

# ggsave(
#   filename = "trees/plots/NJ_alpha.pdf", plot = njtree_alpha,
#   height = 13, width = 10, device = cairo_pdf
# )
# ggsave(
#   filename = "trees/plots/NJ_alpha.png", plot = njtree_alpha,
#   height = 13, width = 10, device = "png", dpi = 350
# )

# write nwk file of nj tree
write.tree(as.phylo(njtree_alpha), file = "trees/exploratory/nj_alpha.nwk", append = FALSE, digits = 10, tree.names = FALSE)

##### UPGMA ALPHA ######

upgma_tree <- phangorn::upgma(matrix_sim)

treeupgma <- groupOTU(upgma_tree, .node = c("M. truncatula"))

upgma_alpha <- ggtree(treeupgma, aes(color = group)) +
  scale_color_manual(values = c("black", "firebrick")) +
  geom_tiplab(hjust = -0.09, fontface = "italic") +
  ggtitle("UPGMA Tree of TDP1alpha") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  ggplot2::xlim(NA, 0.5)

tip_metadataupgma <- tip_metadataalpha[match(upgma_alpha$data$label, tip_metadataalpha$Name), ]

upgma_alpha <- upgma_alpha + geom_tippoint(aes(
  fill = as.character(tip_metadataupgma$Color_lab),
  shape = as.factor(tip_metadataupgma$PCH_lab)
), col = "black", size = 4) +
  scale_fill_manual(values = c("slateblue4", "yellowgreen", "#880085", "brown", "#ff6700", "aquamarine3")) +
  scale_shape_manual(values = c(22, 21, 21, 23, 25, 24))

# uglier and without annotations because we decided to onlyh keep nj

# ggsave(
#   filename = "trees/plots/UPGMA_alpha.pdf", plot = upgma_alpha,
#   height = 13, width = 10, device = cairo_pdf
# )
# ggsave(
#   filename = "trees/plots/UPGMA_alpha.ong", plot = upgma_alpha,
#   height = 13, width = 10, device = 'png', dpi = 350
# )

# write nwk file of upgma
write.tree(as.phylo(upgma_alpha), file = "trees/exploratory/upgma_alpha.nwk", append = FALSE, digits = 10, tree.names = FALSE)

# plot grid
grid_alpha <- cowplot::plot_grid(upgma_alpha, njtree_alpha, labels = "AUTO")

ggsave(
  filename = "trees/exploratory/alpha_trees.pdf", plot = grid_alpha,
  width = 20, height = 13, device = cairo_pdf
)

ggsave(
  filename = "trees/exploratory/alpha_trees.png", plot = grid_alpha,
  width = 20, height = 13, device = 'png', dpi = 350
)

######### BETA NJ/UPGMA TREE #########

# load beta msa  produced from "msa/beta/beta_msa.R"
TDPbeta_aln <- readRDS("msa/beta/TDP1beta_msa_rev.rds")

# compute similarity and identity matrix

# this point is crucial! some other packages that may be useful: "protr", "phanghorn"

matrix_id <- as.matrix(seqinr::dist.alignment(TDPbeta_aln, "identity", gap = 10))
matrix_sim <- as.matrix(seqinr::dist.alignment(TDPbeta_aln, "similarity", gap = 10)) # fitch margulis model is the only option

# however when you give this to mega it won't matter because
# it will use the best sub. model -> see excel

# extract metadate from data/TDP_metadata.xlsx

# build NJ tree from similarity Fitch matrix
nj_tree <- ape::nj(matrix_sim)

treenj <- ggtree::groupOTU(nj_tree, .node = c("M. truncatula"))

##### NJ TREE #####

# extract metadate from data/TDP_metadata.xlsx
tip_metadatabeta <- readxl::read_xlsx(path = "data/TDP_Dataset.xlsx", sheet = 5)

# build NJ tree from similarity Fitch matrix
nj_tree <- ape::nj(matrix_sim)

treenj <- groupOTU(nj_tree, .node = c("M. truncatula"))

nj_beta <- ggtree(treenj, aes(color = group)) +
  scale_color_manual(values = c("black", "firebrick")) +
  geom_tiplab(hjust = -0.06, align = T, linesize = .5, fontface = "italic") +
  geom_tippoint() +
  ggtitle("NJ Tree of TDP1beta") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  ggplot2::xlim(NA, 0.42)

nj_beta <- ggtree::rotate(nj_beta, 41)

tip_metadatanjbeta <- tip_metadatabeta[match(nj_beta$data$label, tip_metadatabeta$Name), ]

nj_beta <- nj_beta +
  geom_tippoint(aes(fill = as.character(tip_metadatanjbeta$Color_lab), shape = as.factor(tip_metadatanjbeta$PCH_lab)),
    col = "black", size = 4
  ) +
  scale_fill_manual(values = c("yellowgreen", "aquamarine3", "#880085")) +
  scale_shape_manual(values = c(21, 21, 24))

# ggsave(
#   filename = "trees/beta/exploratory/NJ_beta.pdf", plot = nj_beta,
#   height = 13, width = 10, device = "pdf"
# )

# write nwk file of nj tree
write.tree(as.phylo(nj_beta), file = "trees/exploratory/nj_beta.nwk", append = FALSE, digits = 10, tree.names = FALSE)

# save tree plot
# saveRDS(nj_beta, file = "trees/both/nj_beta.rds")

###### UPGMA beta #####

upgma_treebeta <- phangorn::upgma(matrix_sim)

treeupgmabeta <- groupOTU(upgma_treebeta, .node = c("M. truncatula"))

upgma_beta <- ggtree(treeupgmabeta, aes(color = group)) +
  geom_tiplab(hjust = +1.2, fontface = "italic") +
  scale_color_manual(values = c("black", "firebrick")) +
  ggtitle("UPGMA Tree of TDP1beta") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  scale_x_reverse() +
  # geom_text(aes(label=node)) +
  scale_y_reverse()

upgma_beta <- ggtree::rotate(upgma_beta, 38)
upgma_beta <- ggtree::rotate(upgma_beta, 40)
upgma_beta <- upgma_beta + geom_strip(3, 6,
  label = "Monocots", align = T, angle = 30, color = "aquamarine3",
  offset = 0.16, barsize = 1, offset.text = 0.043, fontsize = 4
) +
  geom_strip(1, 13,
    label = "Ancestral Plants", align = T, angle = 30, color = "#880085",
    offset = 0.16, barsize = 1, extend = 0.15, offset.text = .067, fontsize = 4
  ) +
  geom_strip(29, 33,
    label = "Eudicots", align = T, angle = 30, color = "yellowgreen",
    offset = 0.16, barsize = 1, offset.text = 0.041, fontsize = 4
  )

tip_metadataupgmabeta <- tip_metadatabeta[match(upgma_beta$data$label, tip_metadatabeta$Name), ]

upgma_beta <- upgma_beta +
  geom_tippoint(aes(fill = as.factor(tip_metadataupgmabeta$Color_lab), shape = as.factor(tip_metadataupgmabeta$PCH_lab)), col = "black", size = 4) +
  scale_fill_manual(values = c("yellowgreen", "aquamarine3", "#880085")) +
  scale_shape_manual(values = c(21, 21, 24))

# ggsave(
#   filename = "trees/beta/exploratory/upgma_beta.pdf", plot = upgma_beta,
#   height = 13, width = 10, device = "pdf"
# )
# ggsave(
#   filename = "trees/beta/exploratory/upgma_beta.png", plot = upgma_beta,
#   height = 13, width = 10, device = "png", dpi = 350
# )

# write nwk file of nj tree
write.tree(as.phylo(upgma_beta), file = "trees/exploratory/upgma_beta.nwk", append = FALSE, digits = 10, tree.names = FALSE)

# save tree plot
# saveRDS(upgma_beta, file = "trees/both/upgma_beta.rds")

grid_beta <- cowplot::plot_grid(upgma_beta, nj_beta, labels = "AUTO")

ggsave(
  filename = "trees/exploratory/beta_trees.pdf", plot = grid_beta,
  width = 20, height = 13, device = cairo_pdf
)
ggsave(
  filename = "trees/exploratory/beta_trees.png", plot = grid_beta,
  width = 20, height = 13, device = 'png', dpi = 350
)