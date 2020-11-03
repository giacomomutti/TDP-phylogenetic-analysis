######## TDP1 an##########

library(ape)
library(pastis)
#library(patchwork)

#Load the functions
source("./msa/script/msa_fun.R")

# Subset for species to display on prettyprint (maybe think to an automated way)
sub <- vector("list", 2)
names(sub) <- c("alpha", "beta")
sub[["alpha"]] <- c(2, 85, 63, 77, 37, 71)
sub[["beta"]] <- c(1, 4, 14, 11)

#Select the dom to do the analysis
dom <- c("alpha","beta")

for (el in dom){
  #Folder and file name
  folder_out=paste0("msa/",el,"/")
  prefix_MEGA=paste0("TDP1",el,"_msa_rev")
  prefix_PP=paste0("TDP1",el,"_PP")
  # Launch the functions
  # Generate a vector of Fasta seq with name
  fasta_vec <- FastaListGen(paste0("data/TDP1_",el,"_fasta/"))

  # Input to generate Pretty Print .pdf and MEGA
  l <- PP_Meg_InpGen(fasta_vec)

  # Generate the input for MEGA and the pdf for PrettyPrint

  # write fasta of msa -> this will be used by MEGA
  write.fasta(as.list(l$MEGA$seq), names = l$MEGA$nam,
              file.out = paste0(folder_out,prefix_MEGA,".fasta"),
              as.string = F)

  # write.msa(x = TDPalpha_aln_seqinr, format = "FASTA", file = "TDP1alpha_msa_v2.fasta")
  saveRDS(l$MEGA, file = paste0(folder_out,prefix_MEGA,".rds"))
  saveRDS(l$PP, file = paste0(folder_out,prefix_PP,".rds"))
  
  # create image of alignment of representative species, this needs latex!
  # Subset need ot be changed according to the species used
  prettyprint <- msa::msaPrettyPrint(l$PP,
                                     file = paste0(folder_out,prefix_PP,".tex"),
                                     subset = sub[[el]], # sig. species
                                     output = "tex",
                                     showNames = "left",
                                     showLogo = "top",
                                     logoColors = "rasmol",
                                     shadingMode = "similar",
                                     showConsensus = "bottom",
                                     showLegend = T,
                                     askForOverwrite = FALSE,
                                     furtherCode = c("\\defconsensus{.}{lower}{upper}", "\\showruler{1}{top}"))
}

alpha_msa <- readRDS('msa/alpha/TDP1alpha_PP.rds')
beta_msa <- readRDS('msa/beta/TDP1beta_PP.rds')

# p1 <- viz_msa(alpha_msa, 'Alpha')
# p2 <- viz_msa(beta_msa, 'Beta')

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
write.nexus.data(x$aabin, format = "protein", file = "Mrbayes/both.nex")

l <- PP_Meg_InpGen(fasta_tot)
# write fasta of msa -> this will be used by MEGA
write.fasta(as.list(l$MEGA$seq), names = l$MEGA$nam,
            file.out = paste0("Mrbayes/both.fasta"),
            as.string = F)

# p3 <- viz_msa(x$PP, 'Both')

### MSAs of both proteins present in both species ###

fasta_vec_alpha <- FastaListGen(paste0("data/TDP1_",'alpha',"_fasta/"))
na <- names(fasta_vec_alpha)

fasta_vec_beta <- FastaListGen(paste0("data/TDP1_",'beta',"_fasta/"))
nb <- names(fasta_vec_beta)

common_species <- nb[nb %in% na]

a_common <- PP_Meg_InpGen(fasta_vec_alpha[names(fasta_vec_alpha) %in% common_species])
b_common <- PP_Meg_InpGen(fasta_vec_beta[names(fasta_vec_beta) %in% common_species])

# p4 <- viz_msa(a_common$PP, 'Alpha in common')
# p5 <- viz_msa(b_common$PP, 'Beta in common')

# p_tot <- p1 / p2 / p3 / p4 / p5 & theme(legend.position = "bottom")
# 
# p_tot <- p_tot + plot_layout(guides = 'collect')

list_pp <- list('Alpha' = alpha_msa,
                'Beta' = beta_msa,
                'Both' = x$PP,
                'Alpha_common' = a_common$PP,
                'Beta_common' = b_common$PP)

p_tot <- viz_msa2(list_pp = list_pp)

ggsave(plot = p_tot, filename = 'msa/plots/msa_viz.pdf', device = cairo_pdf,
       width = 8.5, height = 11)
ggsave(plot = p_tot, filename = 'msa/plots/msa_viz.png', device = png,
       width = 8.5, height = 11, dpi = 300)


d_plot <- density_msas(list_pp)

ggsave(plot = d_plot, filename = 'msa/plots/densities_msa.pdf', device = cairo_pdf)
ggsave(plot = d_plot, filename = 'msa/plots/densities_msa.png', device = 'png',
       dpi = 600)

save.image(file='msa/msa_env.RData')

# check if alpha withoutt animals and algae changes

animals <- c('S. cerevisiae','H. sapiens', 'D. rerio', 'C. brachyrhynchos', 'C. picta bellii', 'X. tropicalis', 'H. glaber')
algae <- c('C. subellipsoidea C-169', "T. socialis","C. sorokiniana","M. conductrix")
ani_al <- c(algae, animals)

fasta_al_ani <- fasta_vec_al[-as.vector(sapply(ani_al, grep, names(fasta_vec_al)))]
fasta_al_ani <- gsub("J", "",gsub("B","",fasta_al_ani))

fasta_anim <- fasta_vec_al[-as.vector(sapply(animals, grep, names(fasta_vec_al)))]
fasta_anim <- gsub("J", "",gsub("B","",fasta_anim))

alpha_al_ani <- PP_Meg_InpGen(fasta_al_ani)
alpha_anim <- PP_Meg_InpGen(fasta_anim)

noanim <- list("Alpha no animals" = alpha_anim$PP, 'Alpha no algae+animals' = alpha_noani_al$PP)

all <- c(list_pp,noanim)
plt3 <- viz_msa2(all)
plt4 <- density_msas(all)

ggsave(plot = plt3, filename = 'msa/plots/msa_viz_animals.pdf', device = cairo_pdf,
       width = 8.5, height = 14)
ggsave(plot = plt4, filename = 'msa/plots/densities_msa_animals.pdf', device = cairo_pdf)


# check if this approach makes sense using cleaned msa with gblocks server
# with the most stringetn parameters

# this should act as control as it should be much more conserved
both_red <- readAAMultipleAlignment("Mrbayes/data/both_red.fasta")

new_l <- c(list_pp, 'Both_red' = both_red)

ctrl <- viz_msa2(new_l)
ctrl2 <- density_msas(new_l)

ggsave(plot = ctrl, filename = 'msa/plots/msa_viz_ctrl.pdf', device = cairo_pdf,
       width = 8.5, height = 11)
ggsave(plot = ctrl2, filename = 'msa/plots/densities_msa_ctrl.pdf', device = cairo_pdf)