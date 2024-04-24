library(tidyverse)
library(fgsea)
library(biomaRt)
library(epitools)

### Read Data

metaData <- readxl::read_xlsx("ALS_Consortium_RNA_Metadata_01042021.xlsx")
metaData <- metaData %>% 
  mutate(ExternalSampleId = gsub("-", "_", ExternalSampleId))
MC_DEA_res <- read.csv("15Sep22_MC_HighHERVPtsvCtrls_CompModel_DEA.csv")
MC_PtOnly_DEA_res <- read.csv("15Sep22_HighHERVPtsvOtherPts_CompModel_DEA.csv")
MC_DEA_meta <- read.csv("15Sep22_MC_HERVKPts_MetaData.csv")

AllCtx_DEA_res <- read.csv("18Sep22_AllCtx_HighHERVPtsvCtrls_CompModel_DEA.csv")
AllCtx_DEA_meta <- read.csv("18Sep22_AllCtx_HERVKPts_MetaData.csv")

####### Do GSEA

all_genes <- AllCtx_DEA_res$X

HERVK_All_Genes <- grep("ERVK", all_genes, value = TRUE)

HERVK_Env_Genes <- grep("env", all_genes, value = TRUE)


ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507


synaptic_vesicle_go <- getBM(attributes=c('ensembl_gene_id'),
                             filters = 'go', values = 'GO:0008021', mart = ensembl)[,1]

synaptic_vesicle_go_rs <- sample(synaptic_vesicle_go, size = 10)

synaptic_vesicle_go_rs_50 <- sample(synaptic_vesicle_go, size = 50)

embryonic_limb_go <- getBM(attributes=c('ensembl_gene_id'),
                           filters = 'go', values = 'GO:0030326', mart = ensembl)[,1]

ltp_go <- getBM(attributes=c('ensembl_gene_id'),
                filters = 'go', values = 'GO:0060291', mart = ensembl)[,1]

ltd_go <- getBM(attributes=c('ensembl_gene_id'),
                filters = 'go', values = 'GO:0060292', mart = ensembl)[,1]

inhibitory_synapse_go <- getBM(attributes=c('ensembl_gene_id'),
                               filters = 'go', values = 'GO:0060077', mart = ensembl)[,1]

excitatory_synapse_go <- getBM(attributes=c('ensembl_gene_id'),
                               filters = 'go', values = 'GO:0060076', mart = ensembl)[,1]

presynaptic_az_go <- getBM(attributes=c('ensembl_gene_id'),
                           filters = 'go', values = 'GO:0048786', mart = ensembl)[,1]

postsynaptic_membrane_go <- getBM(attributes=c('ensembl_gene_id'),
                                  filters = 'go', values = 'GO:0045211', mart = ensembl)[,1]

neuron_apoptosis_go <- getBM(attributes=c('ensembl_gene_id'),
                             filters = 'go', values = 'GO:0051402', mart = ensembl)[,1]

neuron_death_go <- getBM(attributes=c('ensembl_gene_id'),
                         filters = 'go', values = 'GO:1901214', mart = ensembl)[,1]



pathway_list <- list(
  synaptic_vesicle_go = synaptic_vesicle_go,
  HERVK_All = HERVK_All_Genes,
  HERVK_Env = HERVK_Env_Genes,
  neuro_apoptosis = neuron_apoptosis_go,
  neuro_death = neuron_death_go,
  presynaptic = presynaptic_az_go,
  postsynaptic = postsynaptic_membrane_go
)

# MC GSEA
data <- MC_DEA_res %>% 
  dplyr::mutate(Ranks = as.numeric(log2FoldChange * (-log(pvalue, base = 10)))) %>% 
  dplyr::arrange(Ranks) 
        
Ranks <- data$Ranks

names(Ranks) = data$X


set.seed(1)

fgseaRes = fgsea(pathways = pathway_list, 
                 stats = Ranks)
fgseaRes$sig = ifelse(fgseaRes$padj<0.05,"significant","not significant")



fgseaRes$NES = as.numeric(fgseaRes$NES)


ggplot(fgseaRes, aes(x=pathway,y=NES)) +
  geom_col(width=0.8, position='dodge', aes(fill = sig)) +
  labs(title="FGSEA Results Motor Cortex HERV Pts v Ctrls",
       x="Gene Set", y = "Normalized Enrichment Score") +
  
  scale_color_gradient(low="black", high="orange") +
  theme_classic() 
  
ggsave("MC_GSEA_PtsvCtrl.jpeg", width = 12, height = 8,
       device="jpeg")

# Pts Only
data <- MC_PtOnly_DEA_res %>% 
  dplyr::mutate(Ranks = as.numeric(log2FoldChange * (-log(pvalue, base = 10)))) %>% 
  dplyr::arrange(Ranks) 

Ranks <- data$Ranks

names(Ranks) = data$X

fgseaRes = fgsea(pathways = pathway_list, 
                 stats = Ranks)
fgseaRes$sig = ifelse(fgseaRes$padj<0.05,"significant","not significant")



fgseaRes$NES = as.numeric(fgseaRes$NES)


ggplot(fgseaRes, aes(x=pathway,y=NES)) +
  geom_col(width=0.8, position='dodge', aes(fill = sig)) +
  labs(title="FGSEA Results Motor Cortex HERV Pts v Other Pts",
       x="Gene Set", y = "Normalized Enrichment Score") +
  
  scale_color_gradient(low="black", high="orange") +
  theme_classic() 

ggsave("MC_GSEA_PtsOnly.jpeg", width = 12, height = 8,
       device="jpeg")

# All Ctx GSEA

data <- AllCtx_DEA_res %>% 
  dplyr::mutate(Ranks = as.numeric(log2FoldChange * (-log(pvalue, base = 10)))) %>% 
  dplyr::arrange(Ranks) 

Ranks <- data$Ranks

names(Ranks) = data$X


set.seed(1)

fgseaRes = fgsea(pathways = pathway_list, 
                 stats = Ranks)
fgseaRes$sig = ifelse(fgseaRes$padj<0.05,"significant","not significant")



fgseaRes$NES = as.numeric(fgseaRes$NES)


ggplot(fgseaRes, aes(x=pathway,y=NES)) +
  geom_col(width=0.8, position='dodge', aes(fill = sig)) +
  labs(title="FGSEA Results All Cortex",
       x="Gene Set", y = "Normalized Enrichment Score") +
  
  scale_color_gradient(low="black", high="orange") +
  theme_classic() 

ggsave("All_Ctx_GSEA.jpeg", width = 12, height = 8,
       device="jpeg")

########## Check Meta Data 
metaDataPts <- read.csv("HERVK_Pt_Meta_data.csv")


Sex_HERV_Tab <- metaDataPts %>% 
  dplyr::count(Sex,IsHERVK) %>% 
  tidyr::spread(key = IsHERVK, value = n) %>% 
  as.data.frame()

dd <- Sex_HERV_Tab %>% 
  dplyr::select(Other_Pt,HERV_Pt) %>% 
  as.matrix()
rownames(dd) <- Sex_HERV_Tab[,1]

epitools::oddsratio(x = dd)

Onset_HERV_Tab <- metaDataPts %>% 
  dplyr::count(`Site of Motor Onset`,IsHERVK) %>% 
  tidyr::spread(key = IsHERVK, value = n) %>% 
  na.omit() %>% 
  as.data.frame()

dd <- Onset_HERV_Tab %>% 
  dplyr::select(Other_Pt,HERV_Pt) %>% 
  as.matrix()
rownames(dd) <- Onset_HERV_Tab[,1]

epitools::oddsratio(x = dd)

Age_HERV_Tab <- metaDataPts %>% 
  dplyr::count(Age_Factor,IsHERVK) %>% 
  tidyr::spread(key = IsHERVK, value = n) %>% 
  na.omit() %>% 
  as.data.frame()

dd <- Age_HERV_Tab %>% 
  dplyr::select(Other_Pt,HERV_Pt) %>% 
  as.matrix()
rownames(dd) <- Age_HERV_Tab[,1]

epitools::oddsratio(x = dd)

RIN_HERV_Tab <- metaDataPts %>% 
  dplyr::count(RIN_Factor,IsHERVK) %>% 
  tidyr::spread(key = IsHERVK, value = n) %>% 
  na.omit() %>% 
  as.data.frame()

dd <- RIN_HERV_Tab %>% 
  dplyr::select(Other_Pt,HERV_Pt) %>% 
  as.matrix()
rownames(dd) <- RIN_HERV_Tab[,1]

epitools::oddsratio(x = dd)

Hx_HERV_Tab <- metaDataPts %>% 
  dplyr::count(`Family History of ALS/FTD?`,IsHERVK) %>% 
  tidyr::spread(key = IsHERVK, value = n) %>% 
  na.omit() %>% 
  as.data.frame()

dd <- Hx_HERV_Tab %>% 
  dplyr::select(Other_Pt,HERV_Pt) %>% 
  as.matrix()
rownames(dd) <- Hx_HERV_Tab[,1]

epitools::oddsratio(x = dd)

Dementia_HERV_Tab <- metaDataPts %>% 
  dplyr::count(`MND with Dementia?`,IsHERVK) %>% 
  tidyr::spread(key = IsHERVK, value = n) %>% 
  na.omit() %>% 
  as.data.frame()

dd <- Dementia_HERV_Tab %>% 
  dplyr::select(Other_Pt,HERV_Pt) %>% 
  as.matrix()
rownames(dd) <- Dementia_HERV_Tab[,1]

epitools::oddsratio(x = dd)
