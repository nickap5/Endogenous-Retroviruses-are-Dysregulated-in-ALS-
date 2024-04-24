library(readxl)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)


############################################
# Plot Function
library(scales)
numColors <- 2 # How many colors you need
getColors <- scales::brewer_pal('qual', palette = "Set1") 
names(myPalette) <- c("env", "other")

##############################################
# Get All Counts

AllData <- read.table("count_data.txt", 
                      sep = "\t", header = TRUE)


## Read in data
MetaData = read.csv("meta_data.csv")



## Round to whole numbers
rawCounts <- AllData
rawCounts[,-1] <- round(rawCounts[,-1], digits = 0)

## Filter Low Expressed Genes
keep <- rowSums(rawCounts[,-1]) > 5
table(keep, useNA = 'always')

filtCounts <- rawCounts[keep, ]

countDataMatrix <- as.matrix(filtCounts)



##########################################################
# ColInfo



Age_Breaks <- c(0,45,55,65,75, 200)
Age_Labels <- c("Younger_than_45", "45_to_55", "55_to_65", "65_to_75", "Older_than_75")


RIN_Breaks <- c(0,20,30,40,50, 200)
RIN_Labels <- c("Less_than_20", "20_to_30", "30_to_40", "40_to_50", "More_than_50")

MetaData2 <- MetaData %>% 
  dplyr::mutate(Age_Factor = as.character(cut(as.numeric(.$Age), Age_Breaks,Age_Labels))) %>% 
  dplyr::mutate(Age_Factor = as.factor(tidyr::replace_na(Age_Factor, "Unknown"))) %>% 
  dplyr::mutate(RIN_Factor = as.character(cut(as.numeric(.$RIN), RIN_Breaks,RIN_Labels))) %>% 
  dplyr::mutate(RIN_Factor = as.factor(tidyr::replace_na(RIN_Factor, "Unknown")))


###########################################
# Sub CTX

MetaData2 <- MetaData2 %>% 
  dplyr::filter(grepl("cortex",Tissue_subtype)) %>% 
  dplyr::rename(Status = Class_2)

countDataMatrix2 <- countDataMatrix[,colnames(countDataMatrix) %in% MetaData2$Sample]

MetaData2<-MetaData2[MetaData2$Sample %in% colnames(countDataMatrix2),]

################################################################
# Filter out biological sex-associated features 


dds = DESeqDataSetFromMatrix(countData=countDataMatrix2,
                             colData=MetaData2,design =~Gender)
dds$Gender = relevel(dds$Gender,ref="Male")
dds = DESeq(dds,betaPrior=T)

## Obtain genes significantly different in female (i.e. gender associated) with FDR < 0.05
res = results(dds)
sig = res[! is.na(res$padj) & res$padj<0.05,]

vsd = vst(dds, blind = FALSE)
vstCounts <- assay(vsd)

DESeq2::plotPCA(vsd, intgroup = "Status")
ggsave(filename = "PCA_NotBlind_Status_CTX.jpeg")


DESeq2::plotPCA(vsd, intgroup = "Gender")
ggsave(filename = "PCA_NotBlind_Sex_CTX.jpeg")



countDataMatrix2 = countDataMatrix2[! (rownames(countDataMatrix2) %in% rownames(sig)),]


##################################################################
# Create a DESeq Data object from count matrix looking at ALS v Ctrl

MetaData2$C9ORF72 = as.factor(MetaData2$C9ORF72)
dds <- DESeqDataSetFromMatrix(countData=countDataMatrix2, 
                              colData=MetaData2, 
                              design=~Site + Tissue + C9ORF72 + Age_Factor + RIN_Factor + Status)

dds$Status = relevel(factor(dds$Status),ref="control")
keep <- rowSums(counts(dds)) > 5 
dds.filt <- dds[keep,]


ddsObj <- DESeq(dds.filt)
ddsObj <- estimateSizeFactors(ddsObj)

vsd = vst(ddsObj, blind = FALSE)
vstCounts <- assay(vsd)

DESeq2::plotPCA(vsd, intgroup = "Status")
ggsave(filename = "PCA_NotBlind_SexFilt_Status_CTX_.jpeg")


DESeq2::plotPCA(vsd, intgroup = "Gender")
ggsave(filename = "PCA_NotBlind_SexFilt_Sex_CTX_.jpeg")



