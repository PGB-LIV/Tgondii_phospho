########################
# Tgondii manuscript   #  
# Figure 1             #
######################## 
# https://github.com/PGB-LIV/Rice_Phospho_Manuscript/blob/main/Figs/Rice_fig1.R

library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)
library(dplyr)
library(stringr)
library(patchwork)

safe_colorblind_palette <- c("#6699CC","#117733", "#AA4499","#CC6677","#DDCC77")

#############################################
# Collapse experiments to PXD level (sum?)  #
#############################################

##########################
# A. Peptidoform counts  #  
##########################
FLR<-read.csv("data/all_datasets_merged_Site_Peptidoform_centric_Uniprot_tgondii.tsv",sep="\t")
table(FLR$Source.Dataset.Identifier)

# remove contam
prefixes <- c("KAF", "TGVEG", "TGME49", "TGGT1")
pattern <- paste0("^(", paste(prefixes, collapse = "|"), ")")

FLR <- FLR %>%
  filter(
    str_split(Proteins, ":") %>%
      sapply(function(x) any(str_detect(str_trim(x), regex(pattern, ignore_case = TRUE))))
  )

# DF of all peptidoforms per PXD (including decoy sites)
overall<-as.data.frame(table(FLR$Source.Dataset.Identifier))
colnames(overall)<-c("Source.Dataset.Identifier", "sum")
overall$Count<-"Overall"

# Remove decoy sites
FLR<-subset(FLR,Decoy.Modification.Site==0)

# DF of counts of peptidoforms at 5% FLR
peptido_5<-FLR%>% group_by(Source.Dataset.Identifier)%>%
  summarize(sum=sum(Site.Passes.Threshold..0.05.))
peptido_5$Count<-"0.05 FLR"

# DF of counts of peptidoforms at 1% FLR
peptido_1<-FLR%>% group_by(Source.Dataset.Identifier)%>%
  summarize(sum=sum(Site.Passes.Threshold..0.01.))
peptido_1$Count<-"0.01 FLR"

# Combine peptidoform counts
peptidoform_overall<-rbind(overall,peptido_5,peptido_1)

peptido<-ggplot2::ggplot(peptidoform_overall, aes(fill=Count, y=as.numeric(sum), x=Source.Dataset.Identifier)) + geom_bar(position='dodge', stat='identity')+
  theme(axis.text.x = element_text(angle = 45 , hjust=1))+ ggtitle(expression(paste("A. Overall ", italic("Toxoplasma gondii"), " phosphorylation build")))+
  scale_fill_manual(values=safe_colorblind_palette[3:7], name="Peptidoform-site count")+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+xlab("Data set")+ylab("Count")+theme(text = element_text(size=18))+
  scale_y_continuous(labels=scales::comma)

#https://tidytales.ca/snippets/2022-12-22_patchwork-shared-axis-labels/#shared-y-axis-labels
#https://www.data-imaginist.com/posts/2024-01-05-patchwork-1-2-0/

######################################################
#C Gold silver bronze counts per site, with alanines #
# reference ME49                                     #  
######################################################
counts_A<-read.csv("data/Data_S1_G2S1B_0.05_protein_pos_all_prot_mapping_28062026_ptmx.csv")
me49_counts<-counts_A %>% 
  filter(str_detect(Protein, "TGME49"))
safe_colorblind_palette <- c("#888888","#6699CC","#117733", "#CC6677")
cat_counts<-subset(me49_counts,select=c("PTM_FLR_category","PTM_residue"))
cat_counts<-melt(table(cat_counts))

level_order<-c("Bronze","Silver","Gold")

me49 <-ggplot(cat_counts,aes(fill=PTM_residue, y=value, x=factor(PTM_FLR_category,level=level_order)))+ 
  geom_bar(position='dodge', stat='identity')+ ggtitle("C. GSB phosphosite counts in ME49")+
  geom_text(aes(label = value),size = 4, vjust = -0.4, position = position_dodge(0.9))+
  ylab("Count of sites")+
  xlab("Category")+theme(text = element_text(size=18))+
  scale_fill_manual(values=safe_colorblind_palette,name="Phosphosite residue")+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+scale_y_continuous(labels=scales::comma)

######################################################
#Gold silver bronze counts per site, with alanines   #
# facet by strain                                    #  
######################################################
counts_A<-read.csv("data/Data_S1_G2S1B_0.05_protein_pos_all_prot_mapping_28062026_ptmx.csv")
counts_A$strain <- ifelse(startsWith(counts_A$Protein, "KAF"), "KAF",
                          ifelse(startsWith(counts_A$Protein, "TGME49"), "TGME49",
                                 ifelse(startsWith(counts_A$Protein, "TGGT1"), "TGGT1",
                                        ifelse(startsWith(counts_A$Protein, "TGVEG"), "TGVEG", NA))))

safe_colorblind_palette <- c("#888888","#6699CC","#117733", "#CC6677")
cat_counts<-subset(counts_A,select=c("PTM_FLR_category","PTM_residue", "strain"))
cat_counts<-cat_counts %>% 
  group_by(PTM_FLR_category,PTM_residue, strain) %>% 
  summarise(n = n())
level_order<-c("Bronze","Silver","Gold")

gsb_by_strain<-ggplot(cat_counts,aes(fill=PTM_residue, y=n, x=factor(PTM_FLR_category,level=level_order)))+ 
  geom_bar(position='dodge', stat='identity')+
  geom_text(aes(label = n),size = 4, vjust = -0.4, position = position_dodge(0.9))+
  ylab("Count of sites")+
  xlab("Category")+theme(text = element_text(size=18))+
  scale_fill_manual(values=safe_colorblind_palette,name="Phosphosite residue")+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+scale_y_continuous(labels=scales::comma)+facet_wrap(~strain)

######################################################
#B Gold silver bronze counts per strain              #
######################################################
counts_A<-read.csv("data/Data_S1_G2S1B_0.05_protein_pos_all_prot_mapping_28062026_ptmx.csv")
counts_noA<-filter(counts_A,PTM_residue!="A")
counts_noA$strain <- ifelse(startsWith(counts_noA$Protein, "KAF"), "RH88",
                            ifelse(startsWith(counts_noA$Protein, "TGME49"), "ME49",
                                   ifelse(startsWith(counts_noA$Protein, "TGGT1"), "GT1",
                                          ifelse(startsWith(counts_noA$Protein, "TGVEG"), "VEG", NA))))


safe_colorblind_palette <- c("#888888","#6699CC","#117733", "#CC6677")
cat_counts<-subset(counts_noA,select=c("PTM_FLR_category", "strain"))
#cat_counts<-melt(table(cat_counts))
cat_counts<-cat_counts %>% 
  group_by(PTM_FLR_category, strain) %>% 
  summarise(n = n())
level_order<-c("Bronze","Silver","Gold")

gsb_summary<-ggplot(cat_counts,aes(y=n, x=factor(PTM_FLR_category,level=level_order),fill=factor(PTM_FLR_category,level=level_order)))+ 
  geom_bar(position='dodge', stat='identity')+ ggtitle("B. Summary per strain\n   ")+
  geom_text(aes(label = n),size = 4, vjust = -0.4, position = position_dodge(0.9))+
  ylab("Count of sites")+
  xlab("Category") +theme(text = element_text(size=18))+
  scale_fill_manual(values=c("#b05e27", "#a8a9ad", "#d4af37"), name ="Confidence category")+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+scale_y_continuous(labels=scales::comma)+
  facet_wrap(~factor(strain, levels= c("ME49", "GT1", "RH88", "VEG")), ncol=4)



#################################
#D log10(PSM count) by PTM res  #
# GSB                           #
#################################
counts_A<-read.csv("data/Data_S1_G2S1B_0.05_protein_pos_all_prot_mapping_28062026_ptmx.csv")
me49_counts<-counts_A %>% 
  filter(str_detect(Protein, "TGME49"))
# remove ala
counts_noA<-filter(me49_counts,PTM_residue!="A")
counts_noA$log_PSM<-log10(counts_noA$Sum_of_PSM_counts.5.FLR.)
level_order<-c("Bronze","Silver","Gold")
safe_colorblind_palette <- c("#6699CC","#117733", "#CC6677")

psm_counts<-ggplot(data=counts_noA,aes(y=log_PSM, x=factor(PTM_residue), fill=PTM_residue))+ 
  geom_boxplot()+facet_wrap(~factor(PTM_FLR_category,level=level_order)) + ggtitle("D. PSM counts in ME49")+
  scale_fill_manual(values=safe_colorblind_palette,name="Phosphosite residue")+ylab("log10(PSM count)")+ xlab("Residue")+theme(text = element_text(size=18))+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+scale_y_continuous(labels=scales::comma)


##############
# final plot #
##############

patch_plot<-peptido/gsb_summary/(me49+psm_counts)
ggsave("00_preprocess_and_summary/outputs/tgondii_summary.png",plot=patch_plot,dpi=330, width=20,height=16)

#overall_plot<- ggarrange(peptido, me49,gsb_summary ,psm_counts, ncol=2, nrow=2,labels=c("A.", "B.", "C.", "D.") )
#ggsave("outputs/tgondii_summary.png",plot=overall_plot,dpi=330, width=18,height=16)

#####################
# STY distribution  #
#####################

# Distribution of STY in GSB
table(counts_noA$PTM_FLR_category)
table(counts_noA$PTM_FLR_category)/24360

# gold 
gold <- subset(counts_noA, PTM_FLR_category  == "Gold")
table(gold$PTM_residue)
table(gold$PTM_residue)/8021 #total in gold  

# silver 
silver <- subset(counts_noA, PTM_FLR_category  == "Silver")
table(silver$PTM_residue)
table(silver$PTM_residue)/7756  #total in silver

# bronze 
bronze <- subset(counts_noA, PTM_FLR_category  == "Bronze")
table(bronze$PTM_residue)
table(bronze$PTM_residue)/8583    #total in bronze
