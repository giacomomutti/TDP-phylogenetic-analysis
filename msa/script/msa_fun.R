########## TDP1 ALPHA MSA ##########

# this functions are created to create MSA tdp1 alpha w/ clustal Omega
# it produces:
# MSA.fasta
# picture of the MSA
# MSA.rds

# dependencies

library(seqinr)
library(msa)


# Generate the fasta vector
# read all fasta
FastaListGen <- function(path_fasta){
  
  # path_fasta="data/TDP1_alpha_fasta/"
  # list of fasta files
  raw_files <- list.files(path = path_fasta, pattern = "fasta")
  
  fasta_list <- NULL
  
  for (i in 1:length(raw_files)) {
    fasta_list[gsub(".fasta", "", raw_files[i])] <- seqinr::read.fasta(
      file = paste0(path_fasta, 
                    raw_files[i], 
                    sep = ""
      ),
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
  
  return(tdp1_unlisted)
  
}

#Generate the output for MEGA & rds file
# msa w/ CLustal Omega

PP_Meg_InpGen <- function(fasta_vector, prefix_file){
  
  input_list <- vector("list", 2)
  
  names(input_list) <- c("MEGA", "PP")
  
  TDPalpha_aln <- msa::msa(fasta_vector, "ClustalOmega", type = "protein", verbose = T)
  
  # convert MSA to a seqinr msa
  TDPalpha_aln_seqinr <- msaConvert(TDPalpha_aln, type = "seqinr::alignment")
  
  input_list[["MEGA"]] <- TDPalpha_aln_seqinr
  input_list[["PP"]] <- TDPalpha_aln
  
  
  return(input_list)
  
}


PP_AAbin_InpGen <- function(fasta_vector, prefix_file){
  
  input_list <- vector("list", 2)
  
  names(input_list) <- c("aabin", "PP")
  
  TDPalpha_aln <- msa::msa(fasta_vector, "ClustalOmega", type = "protein", verbose = T)
  
  # convert MSA to a seqinr msa
  TDPalpha_aln_seqinr <- msaConvert(TDPalpha_aln, type = "ape::AAbin")
  
  input_list[["aabin"]] <- TDPalpha_aln_seqinr
  input_list[["PP"]] <- TDPalpha_aln
  
  
  return(input_list)
  
}


