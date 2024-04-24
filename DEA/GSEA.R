library(tidyverse)
library(fgsea)
library(biomaRt)

### Read Data

metaData <- readxl::read_xlsx("ALS_Consortium_RNA_Metadata_01042021.xlsx")
metaData <- metaData %>% 
  mutate(ExternalSampleId = gsub("-", "_", ExternalSampleId))


SC_All_DEA_res <- read.csv("allsc.AllSamps.SexFilt.CompModel.csv")
SC_Females_DEA_res <- read.csv("allsc.Females.CompModel.csv")
SC_Males_DEA_res <- read.csv("Males.CompModel.csv")
VSTCounts <- read.table("29Sep22_ERVK_Names.tsv",sep = "\t")

Ctx_All_DEA_res <- read.csv("allctx.AllSamps.SexFilt.CompModel.csv")
Ctx_Females_DEA_res <- read.csv("allctx.Females.CompModel.csv")
Ctx_Males_DEA_res <- read.csv("allctx.Males.CompModel.csv")



####### Do GSEA

all_genes <- rownames(VSTCounts)

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

# All Pts GSEA
data <- Ctx_All_DEA_res %>% 
  dplyr::mutate(Ranks = as.numeric(log2FoldChange * (-log(pvalue, base = 10)))) %>% 
  dplyr::arrange(Ranks) 

Ranks <- data$Ranks

names(Ranks) = data$X

# No Sig
# fgseaRes = fgsea(pathway_list, Ranks, nperm=10000, maxSize=500)
# Lots Sig for all HERVK
set.seed(1)
# fgseaRes = fgseaMultilevel(pathway_list, Ranks)
fgseaRes = fgsea(pathways = pathway_list, 
                 stats = Ranks)
fgseaRes$sig = ifelse(fgseaRes$padj<0.05,"significant","not significant")



fgseaRes$NES = as.numeric(fgseaRes$NES)


ggplot(fgseaRes, aes(x=pathway,y=NES)) +
  geom_col(width=0.8, position='dodge', aes(fill = sig)) +
  labs(title="FGSEA Results Cortex All Patients vs Controls",
       x="Gene Set", y = "Normalized Enrichment Score") +
  
  scale_color_gradient(low="black", high="orange") +
  theme_classic()  +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face = "bold"),
        plot.title = element_text(face = "bold",
                                  size = 16,
                                  hjust = 0.5))

ggsave("Ctx_GSEA_AllPtsvCtrl.jpeg", width = 14, height = 8,
       device="jpeg")

# Females Only
data <- Ctx_Females_DEA_res %>% 
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
  labs(title="FGSEA Results All Cortex Females Only",
       x="Gene Set", y = "Normalized Enrichment Score") +
  
  scale_color_gradient(low="black", high="orange") +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face = "bold"),
        plot.title = element_text(face = "bold",
                                  size = 16,
                                  hjust = 0.5))

ggsave("AllCtx_GSEA_FemalesOnly.jpeg", width = 14, height = 8,
       device="jpeg")

# Males Only

data <- Ctx_Males_DEA_res %>% 
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
  labs(title="FGSEA Results Cortex Males Only",
       x="Gene Set", y = "Normalized Enrichment Score") +
  
  scale_color_gradient(low="black", high="orange") +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face = "bold"),
        plot.title = element_text(face = "bold",
                                  size = 16,
                                  hjust = 0.5))

ggsave("All_Ctx_MalesOnly_GSEA.jpeg", width = 14, height = 8,
       device="jpeg")


###################################################
# SC

# All Pts GSEA
data <- SC_All_DEA_res %>% 
  dplyr::mutate(Ranks = as.numeric(log2FoldChange * (-log(pvalue, base = 10)))) %>% 
  dplyr::arrange(Ranks) 

Ranks <- data$Ranks

names(Ranks) = data$X

# No Sig
# fgseaRes = fgsea(pathway_list, Ranks, nperm=10000, maxSize=500)
# Lots Sig for all HERVK
set.seed(1)
# fgseaRes = fgseaMultilevel(pathway_list, Ranks)
fgseaRes = fgsea(pathways = pathway_list, 
                 stats = Ranks)
fgseaRes$sig = ifelse(fgseaRes$padj<0.05,"significant","not significant")



fgseaRes$NES = as.numeric(fgseaRes$NES)


ggplot(fgseaRes, aes(x=pathway,y=NES)) +
  geom_col(width=0.8, position='dodge', aes(fill = sig)) +
  labs(title="FGSEA Results Spinal Cord All Patients vs Controls",
       x="Gene Set", y = "Normalized Enrichment Score") +
  
  scale_color_gradient(low="black", high="orange") +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face = "bold"),
        plot.title = element_text(face = "bold",
                                  size = 16,
                                  hjust = 0.5))

ggsave("SC_GSEA_AllPtsvCtrl.jpeg", width = 14, height = 8,
       device="jpeg")

# Females Only
data <- SC_Females_DEA_res %>% 
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
  labs(title="FGSEA Results Spinal Cord Females Only",
       x="Gene Set", y = "Normalized Enrichment Score") +
  
  scale_color_gradient(low="black", high="orange") +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face = "bold"),
        plot.title = element_text(face = "bold",
                                  size = 16,
                                  hjust = 0.5))

ggsave("AllSC_GSEA_FemalesOnly.jpeg", width = 14, height = 8,
       device="jpeg")

# Males Only

data <- SC_Males_DEA_res %>% 
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
  labs(title="FGSEA Results Spinal Cord Males Only",
       x="Gene Set", y = "Normalized Enrichment Score") +
  
  scale_color_gradient(low="black", high="orange") +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14,face = "bold"),
        plot.title = element_text(face = "bold",
                                  size = 16,
                                  hjust = 0.5))

ggsave("All_SC_MalesOnly_GSEA.jpeg", width = 16, height = 8,
       device="jpeg")

