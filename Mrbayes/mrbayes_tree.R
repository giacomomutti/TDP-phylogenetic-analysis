library(treeio)
library(ggtree)

# tree produced by MrBayes v 3.2 ] using the fixed Jones invgamma amino acid model, 
# calculated for 1 million generations, sampling every 500 generations and a
# burnin of 25% for the calculation of the consensus tree


# nexus created from figtree
# tree rerooted using animals as outgroup and then
# exported through figtree
# This way I avoided problem of rerooting the tree in R and
# keeping the correct support values

contre <- read.mrbayes(file = "Mrbayes/data/mrbayes_def.nexus")
contre@phylo$tip.label <- gsub("_", " ", contre@phylo$tip.label)
contre@phylo$tip.label <- sub(" ", ". ", contre@phylo$tip.label)


p1 <- ggtree(contre, size = 0.6) +
  # geom_label(aes(label = node)) +
  geom_tiplab(align = F, linesize = .1, size = 2.5, fontface = "italic") +
  geom_text2(aes(
    label = round(as.numeric(prob), 3),
    subset = as.numeric(prob) < 0.99
  ),
  hjust = 1.1, vjust = -0.4, size = 1.5, color = "firebrick"
  ) +
  ggplot2::xlim(0, 5.5)

p1 <- flip(p1, 163, 127)

offset1 <- 0.35
offsettext1 <- -0.06
fontsize1 <- 4
offset2 <- 0.4
offsettext2 <- 0.02

p1 <- p1 +
  geom_strip(24, 1,
    label = "TDP1\u03B2", align = T, hjust = "center", offset = offset1, fontsize = fontsize1,
    color = "darkgreen", barsize = 2.5, angle = 90, offset.text = offsettext1
  ) +
  geom_strip(37, 118, "TDP1\u03b1",
    align = T, hjust = "center", offset = offset1, fontsize = fontsize1,
    color = "firebrick", barsize = 2.5, angle = 90, offset.text = offsettext1
  ) +
  geom_strip(36, 1, "Ancestral",
    align = T, hjust = "left", barsize = 1.5,
    offset = offset2, offset.text = offsettext2, color = "#880085"
  ) +
  geom_strip(6, 5, "Monocots",
    align = T, hjust = "left", barsize = 1.5,
    offset = offset2, offset.text = offsettext2, color = "deepskyblue4"
  ) +
  geom_strip(32, 24, "Rosids",
    align = T, hjust = "left", barsize = 1.5,
    offset = offset2, offset.text = offsettext2, color = "yellowgreen"
  ) +
  geom_strip(34, 35, "Asterids",
    align = T, hjust = "left", barsize = 1.5,
    offset = offset2, offset.text = offsettext2, color = "#4a6200"
  ) +
  geom_strip(123, 119, "Animals",
    align = T, hjust = "left", barsize = 1.5,
    offset = offset2, offset.text = offsettext2, color = "slateblue4"
  ) +
  geom_cladelabel(118, "Fungi",
    align = T, offset.text = -0.025, color = "brown",
    offset = offset2, barsize = 1.5, extend = 0.5
  ) +
  geom_strip(37, 38, "Algae",
    align = T, hjust = "left", barsize = 1.5,
    offset = offset2, offset.text = offsettext2, color = "#ff6700"
  ) +
  geom_strip(44, 43, "Gymnosperms",
    align = T, hjust = "left", barsize = 1.5,
    offset = offset2, offset.text = offsettext2, color = "violet"
  ) +
  geom_strip(42, 45, "Ancestral Plants",
    align = T, hjust = "left", barsize = 1.5,
    offset = offset2, offset.text = offsettext2, color = "#880085"
  ) +
  geom_strip(53, 60, "Monocots",
    align = T, hjust = "left", barsize = 1.5,
    offset = offset2, offset.text = offsettext2, color = "deepskyblue4"
  ) +
  geom_strip(64, 73, "Superasterids",
    align = T, hjust = "left", barsize = 1.5,
    offset = offset2, offset.text = offsettext2, color = "#4a6200"
  ) +
  geom_strip(111, 92, "Rosids",
    align = T, hjust = "left", barsize = 1.5,
    offset = offset2, offset.text = offsettext2, color = "yellowgreen"
  ) +
  geom_cladelabel(61, "Basal Eudicots",
    align = T, offset.text = -0.025, color = "darkgreen",
    offset = offset2, barsize = 1.5, extend = 0.5
  )

# ggsave(plot = p1, filename = "Mrbayes/plots/mrbayes_both.svg",
#        height = 10, width = 15)

ggsave(
  plot = p1, filename = "Mrbayes/plots/mrbayes_both.pdf", device = cairo_pdf,
  height = 10, width = 15
)
ggsave(
  plot = p1, filename = "Mrbayes/plots/mrbayes_both.png",
  height = 10, width = 15, dpi = 400
)


contre <- read.mrbayes(file = "Mrbayes/data/mrbayes_def.nexus")
# contre@phylo$tip.label <- gsub("_"," ", contre@phylo$tip.label)

p2 <- ggtree(contre) + geom_tiplab(fontface = "italic")
p2 <- msaplot(p2, "Mrbayes/data/both.fasta", offset = 2.2, width = 1.7) +
  theme(legend.position = "none")

ggsave(
  plot = p2, filename = "Mrbayes/plots/mrbayes_msa.png",
  height = 18, width = 15, dpi = 350
)
ggsave(
  plot = p2, filename = "Mrbayes/plots/mrbayes_msa.pdf",
  height = 18, width = 15, device = cairo_pdf
)
