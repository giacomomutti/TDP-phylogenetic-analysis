######## TDP1 an##########

#Load the functions
source("msa/script/function/msa_fun.R")

# Subset for species to display on prettyprint (maybe think to an automated way)
sub <- vector("list", 2)
names(sub) <- c("alpha", "beta")
sub[["alpha"]] <- c(2, 85, 63, 77, 37, 71)
sub[["beta"]] <- c(1, 4, 14, 11)

#Select the dom to do the analysis
dom="beta"

#Folder and file name
folder_out=paste0("msa/",dom,"/")
prefix_MEGA=paste0("TDP1",dom,"_msa_rev")
prefix_PP=paste0("TDP1",dom,"_PP")

# Launch the functions
# Generate a vector of Fasta seq with name
fasta_vec <- FastaListGen(paste0("data/TDP1_",dom,"_fasta/"))

# Input to generate Pretty Print .pdf and MEGA
l <- PP_Meg_InpGen(fasta_vec)

# Generate the input for MEGA and the pdf for PrettyPrint

# write fasta of msa -> this will be used by MEGA
write.fasta(as.list(l$MEGA$seq), names = TDPalpha_aln_seqinr$nam, 
            file.out = paste0(folder_out,prefix_MEGA,".fasta"), 
            as.string = F)

# rphast removed from CRAN so rphast::write.msa does not work
# write.msa(x = TDPalpha_aln_seqinr, format = "FASTA", file = "TDP1alpha_msa_v2.fasta")
saveRDS(l$MEGA, file = paste0(folder_out,prefix_MEGA,".rds"))

# create image of alignment of representative species, this needs latex!
# Subset need ot be changed according to the species used
alpha_prettyprint <- msa::msaPrettyPrint(l$PP, 
                                         file = paste0(folder_out,prefix_PP,".pdf"),
                                         subset = sub[[dom]], # sig. species
                                         output = "pdf", 
                                         showNames = "left", 
                                         showLogo = "top", 
                                         logoColors = "rasmol",
                                         shadingMode = "similar", 
                                         showConsensus = "bottom", 
                                         showLegend = T, askForOverwrite = FALSE,
                                         furtherCode = c("\\defconsensus{.}{lower}{upper}", "\\showruler{1}{top}"))
