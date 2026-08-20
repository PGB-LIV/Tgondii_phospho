################################################
# this script is for the motif analysis        #
# central residue= ST (Y not major in tgondii) #  
#                                              #  
#             ---- Gold ----                   #
#                   &                          #   
#             ---- GSB ----                    #     
################################################
require(rmotifx)
require(ggpubr)
library(stringr)
library(ggplot2)
library(ggseqlogo)
library(UpSetR)
library(dplyr)

# rmotifx background - 15mers from peptides in FDR_output
background <- read.csv(file ="01_motif_analysis/inputs/motif_ST_background.txt", header = FALSE)
bg.seqs <- as.character(background$V1)
write.table(bg.seqs, file = paste0("01_motif_analysis/outputs/background_seq.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

# Gold (motif and enrichment analysis) and GSB (motif only)
file_list<-c("01_motif_analysis/inputs/gold_motif_seqs.txt", "01_motif_analysis/inputs/gsb_motif_seqs.txt")
for (f in file_list){
  foreground <- read.csv(file = f, header = FALSE)
  # file name - G vs GSB
  file_name <- strsplit(f, "/")[[1]][3]
  flr_file <- substring(file_name, 1, nchar(file_name) - 4)
  print(flr_file)
  
  ##############
  # 1. rmotifx #  
  #            #
  ##############
  
  
  # central residue in 15mer
  foreground$central <- substr(foreground$V1, 8, 8)
  background$central <- substr(background$V1, 8, 8)
  # filter for ST 15mers
  foreground_ST <- subset(foreground, central == "S" | central == "T")
  fg.seqs <- as.character(foreground_ST$V1)
  
  #find enriched motifs - ST
  mot <- motifx(fg.seqs, bg.seqs, central.res = 'ST', min.seqs = 20, pval.cutoff = 1e-6)
  write.csv(mot, file = paste0("01_motif_analysis/outputs/", flr_file, "/All_motifs_", flr_file, "_ST_01.csv"), row.names = FALSE) 
  
  # filter foreground for S 15mers
  foreground_S <- subset(foreground, central == "S")
  fg.seqs_S <- as.character(foreground_S$V1)
  # filter background for S 15mers
  background_S <- subset(background, central =="S")
  bg.seqs_S <- as.character(background_S$V1)
  
  #find enriched motifs
  mot_S <- motifx(fg.seqs_S, bg.seqs_S, central.res = 'S', min.seqs = 20, pval.cutoff = 1e-6)
  write.csv(mot_S, file = paste0("01_motif_analysis/outputs/", flr_file, "/All_motifs_", flr_file, "_S_with_02.csv"), row.names = FALSE) 

  # filter foreground for T 15mers
  foreground_T <- subset(foreground, central == "T")
  fg.seqs_T <- as.character(foreground_T$V1)
  # filter background for T 15mers
  background_T <- subset(background, central == "T")
  bg.seqs_T <- as.character(background_T$V1)
  
  #find enriched motifs
  mot_T <- motifx(fg.seqs_T, bg.seqs_T, central.res = 'T', min.seqs = 20, pval.cutoff = 1e-6)
  write.csv(mot_T, file = paste0("01_motif_analysis/outputs/", flr_file, "/All_motifs_", flr_file, "_T_03.csv"), row.names = FALSE) 
  
  # skip if no motifs
  if (is.null(mot) || nrow(mot) == 0) next 
  
  #####################################
  # 2.Compare motifs with ST vs S or T#
  # as the central residue            #
  #####################################
  
  #  remove regex to make comparison
  mot$motif_regex_T <- gsub("\\[S", "", mot$motif)
  mot$motif_regex_T <- gsub("\\]", "", mot$motif_regex_T)
  mot$motif_regex_S <- gsub("T\\]", "", mot$motif)
  mot$motif_regex_S <- gsub("\\[", "", mot$motif_regex_S)
  # combine all as each motif can either have S or T as the central residue
  motif_S_regex <- mot$motif_regex_S
  motif_T_regex <- mot$motif_regex_T
  mot_ST_central <- c(motif_S_regex, motif_T_regex)
  # motifs from 1 residue as the central residue 
  mot_T_central <- mot_T$motif
  mot_S_central <- mot_S$motif
  
  print(paste("In ", flr_file, "these motifs are the same with ST or T as central res:")) 
  print(intersect(mot_ST_central, mot_T_central))
  print(paste("In ", flr_file, "these motifs are in ST but not T:")) 
  print(setdiff(mot_ST_central, mot_T_central))
  print(paste("In ", flr_file, "these motifs are in T but not ST:")) 
  print(setdiff(mot_T_central, mot_ST_central))
  
  print(paste("In ", flr_file, "these motifs are the same with ST or S as central res:")) 
  print(intersect(mot_ST_central, mot_S_central))
  print(paste("In ", flr_file, "these motifs are in ST but not S:")) 
  print(setdiff(mot_ST_central, mot_S_central))
  print(paste("In ", flr_file, "these motifs are in S but not ST:")) 
  print(setdiff(mot_S_central, mot_ST_central))
  
  # combine all motifs and remove dups
  overall <- c(mot_S_central, mot_T_central, mot_ST_central)
  overall <- unique(overall)
  # to df
  overall <- as.data.frame(overall)
  
  # change col names for merge to identify what motifs differ in the rmotifx instances
  colnames(overall) <- "motif"
  mot_S_central <- as.data.frame(mot_S_central)
  colnames(mot_S_central) <- "motif"
  mot_T_central <- as.data.frame(mot_T_central)
  colnames(mot_T_central) <- "motif"
  mot_ST_central <- as.data.frame(mot_ST_central)
  colnames(mot_ST_central) <- "motif"
  
  # flag to indicate it is in S central res
  mot_S_central$is_in_S_central <- TRUE
  overall <- left_join(overall, mot_S_central)
  # flag to indicate it is in T central res
  mot_T_central$is_in_T_central <- TRUE
  overall <- left_join(overall, mot_T_central)
  #flag to indicate it is in ST central res
  mot_ST_central$is_in_ST_central <- TRUE
  overall <- left_join(overall, mot_ST_central)
  # NA-> 0
  overall[is.na(overall)] <- 0
  
  # upset plot to show the differences
  png(paste0("01_motif_analysis/outputs/", flr_file, "/upset_plot_", flr_file, ".png"), width = 16, height=14, units="cm", res=300)
  print(upset(overall, text.scale = 1.5))
  dev.off()

  ############################################################
  # 3. Group motifs into classes                             #
  # https://pmc.ncbi.nlm.nih.gov/articles/PMC3254672/        #
  # continue with ST as central res                          #
  ############################################################
  # Group 1: Proline directed
  mot$pro <- str_sub(mot$motif, 12, 12) == "P"
  
  # Group2: Acidic: ≥5 Glu/Asp at +1 to +6 (acidic)
  mot[c("minus","pos")] <- str_split_fixed(mot$motif, "\\[ST\\]", 2)
  mot$pos_loc <- substr(mot$pos, 1, 6)
  mot$de <- str_count(mot$pos_loc, "D") + str_count(mot$pos_loc, "E")
  mot$acidic_1 <- ifelse(mot$de >=5, TRUE, FALSE)
  
  # Group 3: Basophilic  Arg/Lys at −3 (basic)
  mot$basic <- (str_sub(mot$motif, 5, 5) == "R") |
    (str_sub(mot$motif, 5, 5) == "K")
  
  # Group 2: acidophilic Glu/Asp at +1/+2 or +3 (acidic),
  mot$acidic_2<- 
    (str_sub(mot$motif, 12, 12) %in% c("E", "D"))| 
    (str_sub(mot$motif, 13, 13) %in% c("E", "D"))| 
    (str_sub(mot$motif, 14, 14) %in% c("E", "D"))
  
  # Group 3: Basophilic ≥2 Arg/Lys at −6 to −1 (basic)
  mot$minus_loc <- substr(mot$minus, 2, 7)
  mot$rk <- str_count(mot$minus_loc, "R") + str_count(mot$minus_loc, "K")
  mot$basic_2 <- ifelse(mot$rk > 1, TRUE, FALSE)
  
  mot$motif_class <- ifelse(mot$pro,"Proline-directed",
                    ifelse(mot$acidic_1,"Acidophilic",
                           ifelse(mot$basic,"Basophilic",
                                  ifelse(mot$acidic_2,"Acidophilic",
                                         ifelse(mot$basic_2,"Basophilic",
                                                "Other")))))
  counter <- 1
  # seqlogo
  seq_plots <- list()
  seq_plots_bits <- list()
  groups <- c("Acidophilic", "Basophilic", "Proline-directed", "Other") 
  for (motif_group in groups){
    ######################################### 
    # 4. PTM count in each motif class      #
    #########################################
    # filter rmotifx result
    # Get the motifs that fall into a group
    mot_filtered <- subset(mot, motif_class == motif_group)
    motifs_in_group <- mot_filtered$motif
    if (length(motifs_in_group) == 0) next
    print(motif_group)
    print(mot_filtered$motif)
    
    
    # Filter foreground by motifs in the motif group (ST foreground)
    foreground_motif <- filter(foreground_ST, grepl(paste(motifs_in_group, collapse='|'),V1))
    # in each motif class, get the proteins in the foreground
    # unique Protein PTM that have motif
    print("unique PTM sites in class:")
    print(foreground_motif %>% 
            distinct(V2, V3) %>% 
            count())
    #unique proteins
    print("unique proteins in class:")
    print(foreground_motif %>% 
            distinct(V2) %>% 
            count())
    ################################################
    #5. motif seqlogo of enriched motif sequences  #
    ################################################
    pos <- c("-7", "-6", "-5", "-4", "-3", "-2", "-1", "p", "1", "2", "3", "4", "5", "6", "7")
    break_list <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)

    # seq
    foreground_motif_seq <- foreground_motif$V1
    if (length(foreground_motif_seq) == 0) next
    # csv
    write.table(foreground_motif_seq, file = paste0("01_motif_analysis/outputs/",flr_file,"/",motif_group, "_15mer_", flr_file, ".txt"), quote = FALSE, row.names = FALSE,col.names = FALSE)
    
    # seq logo
    # probability (y axis)
    motif_plot <- ggplot() + geom_logo(foreground_motif_seq, method = "probability") + scale_x_continuous(labels = pos, breaks = break_list) + 
    ggtitle(motif_group)  + theme_logo()
    # bits (yaxis)
    motif_plot_bits <- ggplot() + geom_logo(foreground_motif_seq) + scale_x_continuous(labels = pos, breaks = break_list) + 
      ggtitle(motif_group) + theme_logo() +
      guides(fill = guide_legend(title = "Chemistry"))
    # results per motif group
    seq_plots[[counter]] <- motif_plot
    seq_plots_bits[[counter]] <- motif_plot_bits
    counter <- counter + 1
    }
  # create plot for all motif classes
  overall_seq_plot <- ggpubr::ggarrange(plotlist=seq_plots, common.legend = TRUE)
  ggsave(file = paste0("01_motif_analysis/outputs/", flr_file, "/rmotifx_all_", flr_file, ".png"), overall_seq_plot, height = 10, width = 12, dpi = 330)
  overall_seq_plot_bits <- ggpubr::ggarrange(plotlist=seq_plots_bits, common.legend = TRUE)
  ggsave(file = paste0("01_motif_analysis/outputs/", flr_file, "/rmotifx_all_", flr_file, "_bits.png"), overall_seq_plot_bits, height = 8, width = 12, dpi = 330)
  
  ##############################################
  # 6. For each motif, what PTM sites have it? #
  #                                            #
  ##############################################
  all_motifs <- mot$motif
  counter <- 0
  motif_ptm_mappings <- list()
  for (motif in all_motifs){
    counter <- counter + 1
    print(motif)
    # filter to get only proteins in foreground that contain the motif
    protein_with_motif <- filter(foreground_ST, grepl(motif, V1))
    # new column with protein id, site and residue
    protein_with_motif$protein_with_ptm_site <- paste0(protein_with_motif$V2, "_", protein_with_motif$central, protein_with_motif$V3)
    # protein id with ptm pos
    protein_with_motif <- protein_with_motif$protein_with_ptm_site
    # keep unique
    protein_with_motif <- unique(protein_with_motif)
    print(length(protein_with_motif))
    all_protein_with_motif <- paste(protein_with_motif, collapse = ":")
    
    # collect proteins for each motif in df
    motif_ptm_mappings[[counter]] <- data.frame(motif = motif, PTM_sites_with_motif = all_protein_with_motif)
  }
  # join results for all motifs in FLR cat
  motif_proteins_df <- do.call(rbind, motif_ptm_mappings)
  # add the primary protein name back to rmotifx results
  flr_motif_proteins <- left_join(mot, motif_proteins_df, by = "motif")
  flr_motif_proteins <- flr_motif_proteins %>% select(-c(motif_regex_T:basic_2))
  # output
  write.csv(flr_motif_proteins, paste("01_motif_analysis/outputs/", flr_file, "/All_motifs_",flr_file,"_with_proteins.csv", sep = ""), row.names = F)
  
}
