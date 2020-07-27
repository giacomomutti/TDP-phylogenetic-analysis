########## TDP1 ALPHA MSA ##########

# this script creates the MSA tdp1 alpha w/ clustal Omega
# it produces:
# MSA.fasta
# picture of the MSA
# MSA.rds

# dependencies
# Import library

library(seqinr)
library(msa)

# Homologous proteins found with BLASTP w/ query sequence NCBI accession: XP_003622687
# -> see data/TDP_metadata.xlsx

# list of fasta files
raw_files <- list.files(path = "data/TDP1_alpha_fasta/", pattern = "fasta")

# read all fasta
fasta_list <- NULL
for (i in 1:length(raw_files)) {
  fasta_list[gsub(".fasta", "", raw_files[i])] <- seqinr::read.fasta(
    file = paste0("data/TDP1_alpha_fasta/", raw_files[i], sep = ""),
    seqtype = "AA", as.string = TRUE, set.attributes = T
  )
}

# extract names of species
alpha_sp <- lapply(fasta_list, attr, "Annot")
species <- gsub(".*\\[", "", gsub("\\]", "", alpha_sp))
family <- sub("([A-Z]).*", "\\1.", species)
specie <- sub(".*? ", "", species)

sp <- paste(family, specie, sep = " ")
length(sp)

# unlist fasta
tdp1_unlisted <- unlist(fasta_list, use.names = T)
# set species names as id
names(tdp1_unlisted) <- sp

# msa w/ CLustal Omega

TDPalpha_aln <- msa::msa(tdp1_unlisted, "ClustalOmega", type = "protein", verbose = T)

# convert MSA to a seqinr msa
TDPalpha_aln_seqinr <- msaConvert(TDPalpha_aln, type = "seqinr::alignment")

# write fasta of msa -> this will be used by MEGA
write.fasta(as.list(TDPalpha_aln_seqinr$seq), names = TDPalpha_aln_seqinr$nam, file.out = "msa/alpha/TDP1alpha_msa.fasta", as.string = F)

# rphast removed from CRAN so rphast::write.msa does not work
# write.msa(x = TDPalpha_aln_seqinr, format = "FASTA", file = "TDP1alpha_msa_v2.fasta")
saveRDS(TDPalpha_aln_seqinr, file = "msa/alpha/TDP1alpha_msa.rds")

# create image of alignment of representative species, this needs latex!
# i have some problems with output=pdf which wuld be better

alpha_prettyprint <- msa::msaPrettyPrint(TDPalpha_aln, file = "msa/alpha/TDP_alpha_aln.tex",
  subset = c(2, 85, 63, 77, 37, 71), # sig. species
  output = "tex", showNames = "left", showLogo = "top", logoColors = "rasmol",
  shadingMode = "similar", showConsensus = "bottom", showLegend = T, askForOverwrite = FALSE,
  furtherCode = c("\\defconsensus{.}{lower}{upper}", "\\showruler{1}{top}"))

# Found it at the end of a code, is it useful???

# Tdp1_msa <- write.msa(TDPalpha_alignment2, "Tdp1_msa.fasta", format = "FASTA")
# fasta <- read.fasta("Tdp1.msa")
# tdp1_alignment.msa <- import.fasta("Tdp1.msa")
# tdp1_alignment.fa <- export.fasta(tdp1_alignment.fa, outfile = "Tdp1_alignment.fasta")