########## TDP1 ALPHA MSA ##########

# this functions are created to create MSA tdp1 alpha w/ clustal Omega
# it produces:
# MSA.fasta
# picture of the MSA
# MSA.rds

# dependencies

library(seqinr)
library(msa)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

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


msa_to_df <- function(msa) {
  
  aa <-  c('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','U','O','B','J','Z','X','-')
  
  mat <- as.matrix(msa)
  
  n_species <- nrow(mat)  
  
  # count gaps
  gaps <- apply(mat, 2, function(x) length(which(x=='-')))/n_species
  # extract most common aminoacid and then count conservation
  # gaps
  
  mat_df <- t(apply(mat, 2, function(x) table(factor(x, levels = aa))))
  
  # no gaps
  mat_df_nogap <- mat_df[,-ncol(mat_df)]
  
  # check which maximum and then assess conservation
  # max_4col <- colnames(mat_df_nogap)[apply(mat_df_nogap,1,which.max)]
  
  # which is maximum
  max_4col <- apply(mat_df_nogap,1,max)/n_species
  
  df <- data.frame('gap' = gaps, 'most_conserved' = max_4col)
  
  return(df)
}

viz_msa <- function(msa, title, method = 'gam'){
  
  msa_df <- msa_to_df(msa)
  n_species <- nrow(as.matrix(msa))
  
  new_df <- msa_df %>% 
    mutate(id = row_number()) %>% 
    pivot_longer(!id, names_to = 'Type', values_to = 'Freq') 
  
  p1 <- ggplot(new_df, aes(x=id, y=Freq, fill=Type, color = Type)) + 
    geom_area(color = NA, alpha = 0.5) + 
    stat_smooth(se = FALSE, size = 1, method = method) + 
    scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
    scale_x_continuous(expand = c(0,0)) +
    labs(title=paste(title, 'n:', n_species), x = 'Site of MSA', subtitle = 'Frequency of gaps and most conserved aminoacid along the alignment') +
    scale_fill_manual(values = c('#b71540','#0a3d62')) + 
    scale_color_manual(values = c('#b71540','#0a3d62')) +
    theme_classic()  + 
    theme(legend.position = 'bottom', 
          axis.title.y = element_blank())
  
  # p2 <- ggplot(new_df, aes(x=type, y=freq, fill=type)) + 
  #   geom_violin() +
  #   scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
  #   scale_x_discrete(expand = c(0,0)) + 
  #   labs(title = title, subtitle = 'Violin of frequency') +
  #   scale_fill_manual(values = c('#b71540','#0a3d62')) + 
  #   theme_classic() + 
  #   theme(legend.position = 'none',
  #         axis.text.y = element_text(angle = 90)) 
  # 
  # p3 <- (p2 | p1) + plot_layout(widths = c(1,3)) 
  return(p1)
}

density_msas <- function(list_pp){
  
  df <- NULL
  
  for (i in seq_along(list_pp)){
    a <- msa_to_df(list_pp[[i]])
    a <- a %>% mutate(Msa = paste(names(list_pp)[i]), n_sp = nrow(list_pp[[i]])) 
    df <- rbind(df,a)
  }
  
  new_df <- df %>% 
    mutate(id = row_number()) %>% 
    pivot_longer(!c(id,Msa,n_sp), names_to = 'Type', values_to = 'Freq')
  
  n_spp <- tapply(new_df$n_sp, new_df$Msa, max)
  
  ggplot(new_df, aes(x=Freq, color=Msa)) + 
    geom_density(alpha = 0.1) + 
    facet_grid(rows = vars(Type), scales = 'free') + 
    scale_color_viridis_d(option = 'D', labels = paste0(names(n_spp),' n:', n_spp)) + 
    theme_minimal() + 
    labs(title = 'Densities of gaps and most conserved site in all MSAs',
         subtitle = 'gap: left to right -> more gaps in the MSA, most_conserved: left to right -> more conservation') +
    theme(legend.position = 'bottom')
}



viz_msa2 <- function(list_pp, method = 'gam'){
  
  df <- NULL
  
  for (i in seq_along(list_pp)){
    a <- msa_to_df(list_pp[[i]])
    a <- a %>% mutate(Msa = paste(names(list_pp)[i]), n_sp = nrow(list_pp[[i]])) 
    df <- rbind(df,a)
  }
  
  # msa_df <- msa_to_df(msa)
  # n_species <- nrow(as.matrix(msa))
  
  new_df <- df %>% 
    group_by(Msa) %>% 
    mutate(num = 1:n()) %>%   
    pivot_longer(!c(Msa,n_sp,num), names_to = 'Type', values_to = 'Freq') 
  
  n_spp <- tapply(new_df$n_sp, new_df$Msa, max)
  label_msa <- paste0(names(n_spp),', n:', n_spp)
  names(label_msa) <- names(n_spp)

  p1 <- ggplot(new_df, aes(x=num, y=Freq, fill = Type, color = Type)) + 
    geom_area(color = NA, alpha = 0.5) + 
    facet_wrap(ncol = 1, ~ Msa, scales = 'free', labeller = labeller(Msa = label_msa)) +
    #facet_grid(ncol = vars(Msa), scales = 'free') +
    stat_smooth(se = FALSE, size = 1, method = method) + 
    scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
    scale_x_continuous(expand = c(0,0)) +
    labs(x = 'Site of MSA', subtitle = 'Frequency of gaps and most conserved aminoacid along the alignment') +
    scale_fill_manual(values = c('#b71540','#0a3d62')) + 
    scale_color_manual(values = c('#b71540','#0a3d62')) +
    theme_minimal()  + 
    theme(legend.position = 'bottom')
  
  # p2 <- ggplot(new_df, aes(x=type, y=freq, fill=type)) + 
  #   geom_violin() +
  #   scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
  #   scale_x_discrete(expand = c(0,0)) + 
  #   labs(title = title, subtitle = 'Violin of frequency') +
  #   scale_fill_manual(values = c('#b71540','#0a3d62')) + 
  #   theme_classic() + 
  #   theme(legend.position = 'none',
  #         axis.text.y = element_text(angle = 90)) 
  # 
  # p3 <- (p2 | p1) + plot_layout(widths = c(1,3)) 
  return(p1)
}
