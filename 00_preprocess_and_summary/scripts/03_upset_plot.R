#########################################
# Overlap of peptidoforms in Tgondii    #
# strains                               #
######################################### 
library(UpSetR)
library(dplyr)
library(stringr)

# https://www.analytics-tuts.com/upset-plot-in-r-2/
peptidoform_df <-read.csv("data/all_datasets_merged_Site_Peptidoform_centric_UniProt_tgondii.tsv",sep="\t")
# FLR filter
peptidoform_df <- filter(peptidoform_df, Site.Passes.Threshold..0.05. == 1)
# remove from Proteins any protein decoys (rev_)
peptidoform_df <- peptidoform_df %>%
  mutate(proteins_col = sapply(str_split(Proteins, ":"), function(x) {
    paste(x[!str_starts(x, "rev_")], collapse = ":")
  }))
# keep t gondii only
peptidoform_df <- peptidoform_df %>%
  mutate(proteins_tg = sapply(str_split(proteins_col, ":"), function(x) {
    paste(x[str_starts(x, regex("^(tg|kaf)", ignore_case = TRUE))], collapse = ":")
  }))

# remove rows that only had contam
peptidoform_df <- filter(peptidoform_df, proteins_tg != "")
# remove decoy peptidoform sites
peptidoform_df <- filter(peptidoform_df, Decoy.Modification.Site == 0)
# distinct peptidoforms sites
peptidoform_df <- distinct(peptidoform_df, Peptidoform, Peptide.Modification.Position, .keep_all = TRUE)
peptidoform_df$ME49 <- as.integer(as.logical(str_detect(peptidoform_df$proteins_tg, "TGME49")))
peptidoform_df$RH88 <- as.integer(as.logical(str_detect(peptidoform_df$proteins_tg, "KAF")))
peptidoform_df$GT1 <- as.integer(as.logical(str_detect(peptidoform_df$proteins_tg, "TGGT1")))
peptidoform_df$VEG <- as.integer(as.logical(str_detect(peptidoform_df$proteins_tg, "TGVEG")))
peptido_upset<- peptidoform_df %>% select(ME49,RH88,GT1, VEG)
# Plot and save
png("00_preprocess_and_summary/outputs/tgondii_upset.png", res=330, height=15, width=22, units="cm")
upset(peptido_upset, mainbar.y.label= "Intersection size", sets.x.label= "Set size (peptidoform sites)", sets.bar.color="#44AA99",
      main.bar.color="#44AA99",matrix.color="#44AA99", text.scale=1.5, point.size = 3, line.size = 1)
dev.off()

