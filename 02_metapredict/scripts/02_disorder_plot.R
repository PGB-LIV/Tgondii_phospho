##################################################
# Disorder plots                                 #
# input: from 01_metapredict_processing_ME49.py  # 
##################################################
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(dplyr)
safe_colorblind_palette <- c("#6699CC","#117733", "#CC6677")


df <- read.csv(paste("03_metapredict/inputs/disorder_processed_tgondii.csv",sep=""))

# GSB sites
phospho <- subset(df,Phospho=="True")
# Column with the amino acid (not with "pS")
phospho$aa<-phospho$PTM_residue
# Add "p" to show phospho
phospho$PTM_residue<-paste("p",phospho$PTM_residue,sep="")
phospho<-subset(phospho,select=c(Protein,PTM_residue,disorder,PTM_FLR_category,aa))

# All other sites not phosphorylated
other_ST<-subset(df, Phospho!= "True")
# ST
other_ST<-subset(other_ST, PTM_residue=="S"|PTM_residue=="T")
other_ST<-subset(other_ST,select=c(Protein,PTM_residue,disorder,PTM_FLR_category))
other_ST$aa <- other_ST$PTM_residue

######################## 
# disorder score stats #
########################

# median
pS<-subset(phospho,PTM_residue=="pS")
quantile(pS$disorder)
pT<-subset(phospho,PTM_residue=="pT")
quantile(pT$disorder)

S_results<-subset(other_ST,PTM_residue=="S")
quantile(S_results$disorder)
T_results<-subset(other_ST,PTM_residue=="T")
quantile(T_results$disorder)


#Kolmogorov-Smirnov 
ks.test(pS$disorder, S_results$disorder, exact = TRUE)
ks.test(pT$disorder, T_results$disorder, exact = TRUE)


##############
# Plots      #
##############

# combine
overall<-rbind(phospho,other_ST)
n_labels <- overall %>%
  group_by(PTM_residue) %>%
  summarise(n = n(), y_pos = max(disorder, na.rm = TRUE) + 0.25, .groups = "drop")


p1<-ggplot(data=overall, mapping=aes(x=PTM_residue, y=disorder,fill=aa))+geom_boxplot()+theme_bw()+labs(title="A. pST versus ST")+
  scale_fill_manual(values=safe_colorblind_palette,name="")+ylab("Disorder score") + xlab("Residue")+ geom_hline(yintercept=0.5, linetype="dashed", color = "#888888",linewidth=1.5)+
  geom_signif(comparisons = list(c("S", "pS"),c("T", "pT")),annotation = c("D= 0.188***","D= 0.275***"), y_position=c(1.0,1.1), textsize = 5) + theme(text = element_text(size=18))+
  geom_text(data = n_labels,
            aes(x = PTM_residue, y = y_pos, label = paste0("n=", n)),
            position = position_dodge(width = 0.75),
            inherit.aes = FALSE, size = 4, vjust = 0) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) 

# Gold-Silver-Bronze
phospho$PTM_FLR_category<-factor(phospho$PTM_FLR_category,levels=c("Gold","Silver","Bronze"))
n_labels <- phospho %>%
  group_by(PTM_FLR_category, PTM_residue, aa) %>%
  summarise(n = n(), y_pos = max(disorder, na.rm = TRUE) + 0.05, .groups = "drop")


p2<-ggplot(data=phospho, mapping=aes(x=PTM_residue, y=disorder,fill=aa))+geom_boxplot()+theme_bw()+labs(title="B. pST in Gold-Silver-Bronze")+facet_wrap(~PTM_FLR_category)+
  scale_fill_manual(values=safe_colorblind_palette,name="")+ylab("Disorder score")+xlab("Residue")+geom_hline(yintercept=0.5, linetype="dashed", color = "#888888",linewidth=1.5)+
  geom_text(data = n_labels,
            aes(x = PTM_residue, y = y_pos, label = paste0("n=", n), group = aa),
            position = position_dodge(width = 0.75),
            inherit.aes = FALSE, size = 4, vjust = 0) +
  theme(text = element_text(size = 18))

ggarrange(p1,p2,ncol=1, common.legend = TRUE)
ggsave(paste("03_metapredict/outputs/tgondii_disorder_pST.png",sep=""),dpi=330, height=12, width=14)


