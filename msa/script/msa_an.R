######## TDP1 an##########

#Load the functions
source("/media/ubuntu/26f3b6d7-0e54-4371-9efd-01b2f0f03216/project/MuttiTDP_phyloan/TDP-phylogenetic-analysis/msa/script/msa_fun.R")

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
prettyprint <- msa::msaPrettyPrint(l$PP, 
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

# Align the two lists
fasta_vec_al <- FastaListGen(paste0("data/TDP1_","alpha","_fasta/"))
fasta_vec_bet <- FastaListGen(paste0("data/TDP1_","beta","_fasta/"))

# Unique list
names(fasta_vec_al) <- paste0(gsub("\\.","",gsub(" ", "_",names(fasta_vec_al))), "_alpha")
names(fasta_vec_bet) <- paste0(gsub("\\.", "",gsub(" ", "_",names(fasta_vec_bet))), "_beta")

# Combine the two list 
fasta_tot <- c(fasta_vec_al, fasta_vec_bet)

# Adjust manually the sequence of protein according mrBayes
fasta_tot <- gsub("J", "",gsub("B","",fasta_tot))

# Create the alignement

x <- PP_AAbin_InpGen(fasta_tot, prefix_file = "tot")

# First of all you need to create the nexus for mrBayes

library(ape)
library(pastis)
write.nexus.data(x$aabin, format = "protein", file = "data/prova.nex")


