########## TDP1 BETA MSA ##########

# this script creates the MSA tdp1 beta w/ clustal Omega
# it produces:
# MSA.fasta
# picture of the MSA
# MSA.rds

# dependencies

library(seqinr)
library(msa)

# Homologous proteins found with BLASTP w/ query sequence NCBI accession: AET04910

# -> see data/TDP_metadata.xlsx

# list of fasta files
raw_files_beta <- list.files(path = "data/TDP1_beta_fasta/", pattern = "fasta")

# read all fasta
fasta_list_beta <- NULL
for (i in 1:length(raw_files_beta)) {
  fasta_list_beta[gsub(".fasta", "", raw_files_beta[i])] <- seqinr::read.fasta(
    file = paste0("data/TDP1_beta_fasta/", raw_files_beta[i], sep = ""), 
    seqtype = "AA", as.string = TRUE, set.attributes = T)
}

# extract names of species
beta_sp <- lapply(fasta_list_beta, attr, "Annot")
species <- gsub(".*\\[", "", gsub("\\]", "", beta_sp))
family <- sub("([A-Z]).*", "\\1.", species)
specie <- sub(".*? ", "", species)

sp <- paste(family, specie, sep = " ")
length(sp)

# unlist fasta
tdp1beta_unlisted <- unlist(fasta_list_beta, use.names = T)
# J is a not standard AA that correspond to XLE (leucine or isoleucine)
# here it is substitued with Leu
tdp1beta_unlisted <- gsub("J", "L", tdp1beta_unlisted)

# set species names as id
names(tdp1beta_unlisted) <- sp

# msa w/ CLustal Omega
TDPbeta_aln <- msa::msa(tdp1beta_unlisted, "ClustalOmega", type = "protein", verbose = T)

# convert MSA ta a seqinr msa -> this will be used by MEGA
TDPbeta_aln_seqinr <- msaConvert(TDPbeta_aln, type = "seqinr::alignment")
write.fasta(as.list(TDPbeta_aln_seqinr$seq), names = TDPbeta_aln_seqinr$nam, file.out = "msa/beta/TDP1beta_msa.fasta", as.string = F)

saveRDS(TDPbeta_aln_seqinr, file = "msa/beta/TDP1beta_msa.rds")

# create image of alignment of representative species, this needs latex!
# i have some problems with output=pdf which wuld be better

beta_prettyprint <- msaPrettyPrint(TDPbeta_aln, file = "msa/beta/TDP_beta_aln.tex",
  subset = c(1, 4, 14, 11),# sig. species
  output = "tex", showNames = "left", showLogo = "top", logoColors = "rasmol", 
  shadingMode = "similar", showConsensus = "bottom", showLegend = T, askForOverwrite = FALSE, 
  furtherCode = c("\\defconsensus{.}{lower}{upper}", "\\showruler{1}{top}"))