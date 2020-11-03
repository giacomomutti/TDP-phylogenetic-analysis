#dependencies 

library(treeio)
library(ggtree)
library(ggplot2)
library(castor)

# extract subset tree from subset of species
# to visualize it beside GSDS and MEME plots

alpha <- read.tree("trees/data/alpha_500bt.nwk")
alpha$tip.label <- gsub("_"," ", alpha$tip.label)

# extract subset species of interest
subset_a <- alpha$tip.label[c(33, 66, 67, 72, 69,
                           16, 61, 11, 12, 14,
                           38, 34, 35, 36, 43,
                           43, 46, 50, 51, 54,
                           55, 28, 29, 31, 25,
                           20, 18, 19, 5, 3,
                           1, 15, 26, 27,
                           75, 74)]

subset_tree_alpha <- get_subtree_with_tips(alpha, only_tips = subset_a)
a <- ggtree(subset_tree_alpha$subtree, branch.length = "none") + ggtree::geom_tiplab() + ggplot2::xlim(0,18)

ggsave(plot = a, filename = "gsds-meme/gsds_subset_alpha.pdf", 
                device = cairo_pdf, width = 3.5, height = 3*3.25)

beta <- read.tree("trees/data/beta_500bt.nwk")
beta$tip.label <- gsub("_"," ", beta$tip.label)
# should be 26 seqs
subset_b <- beta$tip.label[c(36, 32, 33, 34, 29,
                             30, 35, 28, 26, 27, 
                             21, 24, 22, 18, 17, 
                             9, 14, 15, 12, 13, 
                             6, 7, 8, 4, 3, 2)]

subset_tree_beta <- get_subtree_with_tips(beta, only_tips = subset_b)

b <- ggtree(subset_tree_beta$subtree, branch.length = "none") + ggtree::geom_tiplab() + ggplot2::xlim(0,15)

ggsave(plot = b, filename = "gsds-meme/gsds_subset_beta.pdf", 
                device = cairo_pdf, width = 3.5, height = 3*3.25)

# then modify figure with inkscape