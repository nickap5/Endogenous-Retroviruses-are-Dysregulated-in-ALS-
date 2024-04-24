library(readxl)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(ggpubr)

######################################
# Read meta data

MetaData = readxl::read_xlsx("meta_data.xlsx")


##############################################
options(digits = 3)

# Format MetaData3 to have distinct individuals
MetaData3 <- MetaData %>% 
  dplyr::filter(grepl("Cortex|Spinal_Cord",Source)) %>% 
  dplyr::distinct(Subject, .keep_all = TRUE) 


ggplot( MetaData3, aes( x=RIN, y=`Post Mortem Interval in Hours` ))+
  geom_point(pch=10)+
  theme_bw() +
  ggtitle("Spearman Correlation of RIN and PMI") +
  geom_smooth(method='lm', formula= y~x) +
  stat_cor(method = "spearman", label.x = 2.5, label.y = 125) 

ggsave("IndividSamp_RINvPMI_Spearman.jpeg")


ggplot( MetaData3, aes( x=RIN, y=pH ))+
  geom_point(pch=10)+
  theme_bw() +
  ggtitle("Spearman Correlation of RIN and pH") +
  geom_smooth(method='lm', formula= y~x) +
  stat_cor(method = "spearman", label.x = 2.5, label.y = 8) 

ggsave("IndividSamp_RINvpH_Spearman.jpeg")


ggplot( MetaData3, aes( y=`Age at Symptom Onset`, x=`Age at Death` ))+
  geom_point(pch=10)+
  theme_bw() +
  ggtitle("Spearman Correlation of Age at Death and Symptom Onset") +
  geom_smooth(method='lm', formula= y~x) +
  stat_cor(method = "spearman", label.x = 30, label.y = 100) 


ggsave("IndividSamp_DeathvSxOnset_Spearman.jpeg")


ggplot( MetaData3, aes( x=`Age at Death`, y=`Disease Duration in Months` ))+
  geom_point(pch=10)+
  theme_bw() +
  ggtitle("Spearman Correlation of Age at Death and Disease Duration") +
  geom_smooth(method='lm', formula= y~x) +
  stat_cor(method = "spearman", label.x = 30, label.y = 180) 


ggsave("IndividSamp_DeathvDzDur_Spearman.jpeg")

library(Hmisc)
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

res2 <- MetaData3 %>% 
  dplyr::select(`Age at Symptom Onset`,`Age at Death`, RIN,
           pH, `Post Mortem Interval in Hours`, `Disease Duration in Months`) %>% 
  as.matrix()


corr_res <- rcorr(res2,type="spearman")
res3 <- flattenCorrMatrix(corr_res$r, corr_res$P) 
res3 <- res3 %>% 
  mutate(padj= p.adjust(p,method = "fdr")) %>% 
  arrange(padj)
write.csv(res3,"Spearman_Corr_Matrix.csv")


# Take only those samples with complete information
res2 <- MetaData3 %>% 
  dplyr::select(`Age at Symptom Onset`,`Age at Death`, RIN,
         pH, `Post Mortem Interval in Hours`, `Disease Duration in Months`) %>% 
  dplyr::filter(complete.cases(.)) %>% 
  as.matrix()


corr_res <- rcorr(res2,type="spearman")
res3 <- flattenCorrMatrix(corr_res$r, corr_res$P) 
res3 <- res3 %>% 
  mutate(padj= p.adjust(p,method = "fdr")) %>% 
  arrange(padj)
write.csv(res3,"Spearman_Corr_Matrix_CompleteCasesOnly.csv")


Complete_Sites <- MetaData3 %>% 
  dplyr::select(`Age at Symptom Onset`,`Age at Death`, RIN,
                pH, `Post Mortem Interval in Hours`, 
                `Disease Duration in Months`, `Site Specimen Collected`) %>% 
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::count(`Site Specimen Collected`)
write.csv(Complete_Sites,"Complete_Site_Table.csv")



All_Sites <- MetaData3 %>% 
  dplyr::select(`Age at Symptom Onset`,`Age at Death`, RIN,
                pH, `Post Mortem Interval in Hours`, 
                `Disease Duration in Months`, `Site Specimen Collected`) %>% 
  dplyr::count(`Site Specimen Collected`)

write.csv(All_Sites,"All_Site_Table.csv")



####################################################
# Reveal Site Bias
######################################################


tmp <- MetaData3 %>% 
  dplyr::select(c("ExternalSubjectId","Site Specimen Collected")) %>% 
  dplyr::mutate(Site_Abbreviation = abbreviate(`Site Specimen Collected`,1,
                                                      method="both.sides",strict = TRUE))

tmp_2 <- MetaData3 %>% 
  dplyr::select(ExternalSubjectId,`Age at Symptom Onset`,`Age at Death`, RIN,
                pH, `Post Mortem Interval in Hours`, 
                `Disease Duration in Months`, `Site Specimen Collected`) %>% 
  dplyr::filter(complete.cases(.)) %>% 
  dplyr::select(c("ExternalSubjectId","Site Specimen Collected")) %>% 
  dplyr::mutate(Site_Abbreviation = abbreviate(`Site Specimen Collected`,1,
                                                      method="both.sides",strict = TRUE))


library(RColorBrewer)
myColors <- brewer.pal(9,"Set1")
names(myColors) <- unique(tmp$Site_Abbreviation)
colScale <- scale_fill_manual(name = "Site_Abbreviation",values = myColors)

ggplot(tmp, aes(x=Site_Abbreviation,fill = Site_Abbreviation)) + 
  geom_bar() + theme_bw(base_size=14) + colScale +
  ggtitle("Number of Patients per Site for all Patients")
ggsave("Site_Dist_AllPts.jpeg")


ggplot(tmp_2, aes(x=Site_Abbreviation,fill = Site_Abbreviation)) + 
  geom_bar() + theme_bw(base_size=14) + colScale +
  ggtitle("Number of Patients per Site for Patients with Complete Data")
ggsave("Site_Dist_SubsetPts.jpeg")


