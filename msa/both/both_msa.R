library(seqinr)
library(msa)

raw_files_alpha <- list.files(path = "data/TDP1_alpha_fasta/", pattern = "fasta")
raw_files_beta <- list.files(path = "data/TDP1_beta_fasta/", pattern = "fasta")

raw_files <- c(raw_files_alpha,raw_files_beta)

fasta_list_alpha <- NULL

for (i in 1:length(raw_files_alpha)) {
  fasta_list_alpha[gsub(".fasta", "", raw_files_alpha[i])] <- seqinr::read.fasta(
    file = paste0("data/TDP1_alpha_fasta/", raw_files_alpha[i], sep = ""),
    seqtype = "AA", as.string = TRUE, set.attributes = T
  )
}

fasta_list_beta <- NULL


for (i in 1:length(raw_files_beta)) {
  fasta_list_beta[gsub(".fasta", "", raw_files_beta[i])] <- seqinr::read.fasta(
    file = paste0("data/TDP1_beta_fasta/", raw_files_beta[i], sep = ""),
    seqtype = "AA", as.string = TRUE, set.attributes = T
  )
}


# extract names of species
header <- lapply(fasta_list_alpha, attr, "Annot")
species <- gsub(".*\\[", "", gsub("\\]", "", header))
family <- sub("([A-Z]).*", "\\1.", species)
specie <- sub(".*? ", "", species)

sp_alpha <- paste(family, specie, "alpha",sep = " ")


header_b <- lapply(fasta_list_beta, attr, "Annot")
species_b <- gsub(".*\\[", "", gsub("\\]", "", header_b))
family_b <- sub("([A-Z]).*", "\\1.", species_b)
specie_b <- sub(".*? ", "", species_b)

sp_beta <- paste(family_b, specie_b, "beta",sep = " ")

fasta_list <- c(fasta_list_alpha,fasta_list_beta)

tdp1_unlisted <- unlist(fasta_list, use.names = T)
# set species names as id
names(tdp1_unlisted) <- c(sp_alpha,sp_beta)

# msa w/ CLustal Omega

TDP_both <- msa::msa(tdp1_unlisted, "ClustalOmega", type = "protein", verbose = T)

# convert MSA to a seqinr msa
TDP_both_aln_seqinr <- msaConvert(TDP_both, type = "seqinr::alignment")

# write fasta of msa -> this will be used by MEGA

# to do guidabce
TDP_both_aln_seqinr$nam <- gsub("-", "", TDP_both_aln_seqinr$nam)

write.fasta(as.list(TDP_both_aln_seqinr$seq), names = TDP_both_aln_seqinr$nam, file.out = "msa/both/TDP1both_msa.fasta", as.string = F)

# rphast removed from CRAN so rphast::write.msa does not work
# write.msa(x = TDPalpha_aln_seqinr, format = "FASTA", file = "TDP1alpha_msa_v2.fasta")
saveRDS(TDP_both_aln_seqinr, file = "msa/both/TDP1both_msa.rds")
