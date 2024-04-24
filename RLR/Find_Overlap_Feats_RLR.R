library(tidyverse)
library(biomaRt)

# Read in and format data
SC_dat <- read.csv("SC_Norm_CaretVarImp_Top_SiglambdaMin.csv")

CTX_Dat <- read.csv("Norm_CaretVarImp_Top_SiglambdaMin.csv")

CTX_DEA_Dat <- read.csv("allctx.AllSamps.SexFilt.CompModel.csv") %>% 
  dplyr::select(c("X","log2FoldChange","padj"))
SC_DEA_Dat <- read.csv("allsc.AllSamps.SexFilt.CompModel.csv") %>% 
  dplyr::select(c("X","log2FoldChange","padj"))

Overlap_dat <- CTX_Dat %>% 
  dplyr::filter(X %in% SC_dat$X) %>% 
  dplyr::arrange(desc(Overall))
Overlap_dat <- merge(Overlap_dat,SC_dat, by="X", 
                     suffixes = c("_CTX","_SC"))

# reformat rownames
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- Overlap_dat$X
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)

# Find Overlap
x <- merge(Overlap_dat,G_list,by.x = "X", by.y = "ensembl_gene_id") %>% 
  merge(CTX_DEA_Dat, by = "X") %>% 
  dplyr::rename(log2FoldChange_CTX = log2FoldChange) %>% 
  dplyr::rename(padj_CTX = padj) %>% 
  merge(SC_DEA_Dat, by = "X") %>% 
  dplyr::rename(log2FoldChange_SC = log2FoldChange) %>% 
  dplyr::rename(padj_SC = padj)
 

write.csv(x,"Overlap_Top_Sig_LambdaMin_CTX_SC.csv")


