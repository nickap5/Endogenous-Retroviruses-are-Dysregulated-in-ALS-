library(readxl)
library(dplyr)
library(DESeq2)
library(ggplot2)


############################################
# Plot Function
library(scales)
numColors <- 2 # How many colors you need
getColors <- scales::brewer_pal('qual', palette = "Set1") 
myPalette <- getColors(numColors)
names(myPalette) <- c("env", "other")

##############################################
# Get All Counts

## Read in data
AllData <-  read.table("all_data.txt", 
                      sep = "\t", header = TRUE)


## Read in data
MetaData = read.csv("meta_data.csv")



## Round to whole numbers
rawCounts <- AllData
rawCounts[,-1] <- round(rawCounts[,-1], digits = 0)

#### Filter
## Filter Low Expressed Genes
keep <- rowSums(rawCounts[,-1]) > 5
table(keep, useNA = 'always')

filtCounts <- rawCounts[keep, ]

countDataMatrix <- as.matrix(filtCounts)



##########################################################
# ColInfo

#Get info on regions

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
# Sub Ctx

MetaData2 <- MetaData2 %>% 
  dplyr::filter(Tissue_subtype=="cortex") %>% 
  dplyr::rename(Status = Class_2)

countDataMatrix2 <- countDataMatrix[,colnames(countDataMatrix) %in% MetaData2$Sample]

MetaData2<-MetaData2[MetaData2$Sample %in% colnames(countDataMatrix2),]
MetaData2$C9ORF72 = as.factor(MetaData2$C9ORF72)

################################################################
# Separate Females from Males
FemaleSamps <- MetaData2$Sample[MetaData2$Gender=="Female"]
countDataFemale <- countDataMatrix2[,colnames(countDataMatrix2) %in% FemaleSamps]
MetaDataFemale <- MetaData2[MetaData2$Sample %in% FemaleSamps,]

MaleSamps <- MetaData2$Sample[MetaData2$Gender=="Male"]
countDataMale <- countDataMatrix2[,colnames(countDataMatrix2) %in% MaleSamps]
MetaDataMale <- MetaData2[MetaData2$Sample %in% MaleSamps,]


##################################################################
# Create a DESeq Data object from count matrix looking at ALS v Ctrl
# Females

dds <- DESeqDataSetFromMatrix(countData=countDataFemale, 
                                     colData=MetaDataFemale, 
                                     design=~Site + Tissue + C9ORF72 + Age_Factor + RIN_Factor + Status)

dds$Status = relevel(factor(dds$Status),ref="control")
keep <- rowSums(counts(dds)) > 5 
dds.filt <- dds[keep,]

# Find DEGs
ddsObj <- DESeq(dds.filt)
ddsObj <- estimateSizeFactors(ddsObj)
results.dds <- results(ddsObj, alpha = 0.05)
results.dds.sig <- subset(results.dds,!is.na(results.dds$padj))
results.dds.sig.ALS <- subset(results.dds, results.dds$padj  < 0.05)
results.dds.sig.ALS.ENS <- results.dds.sig.ALS[grep("ENS", rownames(results.dds.sig.ALS)),]


write.csv(results.dds.sig.ALS, 
          "allctx.Females.CompModel.csv")

write.csv(results.dds.sig.ALS.ENS, 
          "allctx.Females.LFTab.ENSOnly.CompModel.csv")
results.dds.sig.ALS.ENS.IPA <- results.dds.sig.ALS.ENS %>% 
  as.data.frame() %>% 
  dplyr::arrange(padj) %>% 
  dplyr::slice_head(n = 1550)
write.csv(results.dds.sig.ALS.ENS.IPA, 
          "allctx.Females.ForIPA.LFTab.csv")

# Shrink vals for downstream analysis
library(ashr)
ddsShrink <- lfcShrink(ddsObj,
                       res = results.dds,
                       type = "ashr")

vsd = vst(ddsObj, blind = FALSE)
vstCounts <- assay(vsd)
write.csv(vstCounts, "DispNotBlind_VSTCounts_allctx_females.csv")

jpeg("PCA_Females_CompModel_ntop500.jpg")
DESeq2::plotPCA(vsd,
                intgroup = "Status",
                ntop = 500,
                returnData = FALSE)
dev.off()


resLFC <- lfcShrink(ddsObj, coef="Status_als_vs_control", type="ashr")
jpeg("MAregularresults_Females_CompleteModel.jpg")
plotMA(results.dds, ylim=c(-1,1))
dev.off()

jpeg("MAlfcShrunk_Females_CompleteModel.jpg")
plotMA(resLFC, ylim=c(-1,1))
dev.off()


library("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
ens <- rownames(resLFC)

annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_transcript_id", "ensembl_gene_id",
               "gene_biotype", "external_gene_name"),
  filter="ensembl_gene_id",
  values=ens,
  uniqueRows=TRUE)

annotLookup2 <- annotLookup %>% 
  as.data.frame() %>% 
  distinct(ensembl_gene_id, .keep_all = TRUE) %>% 
  arrange(match(ensembl_gene_id, ens))
ens2 <- substr(ens,1,15)  
ens2[ens2 %in% annotLookup2$ensembl_gene_id] = annotLookup2$external_gene_name


library(EnhancedVolcano)
jpeg("VolcanoPlot_ShrunkFC_Females_CompleteModel.jpg")
  
EnhancedVolcano(resLFC,
                lab = ens2,
                x = 'log2FoldChange',
                title = 'Volcano Plot Female ALS v Control',
                y = 'pvalue',
                xlim = c(-5,5),
                colAlpha = 1)
dev.off()

ERVKFeats <- as.data.frame(results.dds.sig.ALS[grep("ERVK", rownames(results.dds.sig.ALS)),])
feat_of_interest <- grep("env", rownames(ERVKFeats),value = TRUE)

ERVKFeats <- dplyr::mutate(ERVKFeats,
                           axis_col = ifelse(rownames(ERVKFeats) %in% feat_of_interest, "env", "other"))

axis_cols = myPalette[ERVKFeats$axis_col]
axis_cols = axis_cols[order(rownames(ERVKFeats))]

write.csv(ERVKFeats, "AllCtxDEA_Females_FeatPlot_HERVKOnly.csv")
ggplot(ERVKFeats, aes(x = substr(rownames(ERVKFeats),1,10), y = log2FoldChange, col = padj)) +
  geom_point(size = 4) +
  labs(title="HERVK DE Features in Female Cortical Samples",
       x="Gene", y = "Expression (log2FC)") +
  scale_color_gradient(low="black", high="orange") +
  ylim(-1.5,1.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, 
                                   colour = axis_cols),
        axis.text = element_text(size=11),
        plot.title = element_text(size=16))
ggsave("AllCtxDEA_Females_NAP_HERVKOnly.jpeg", 
       device="jpeg")






########################################################
# Males Only
#####################################################

dds <- DESeqDataSetFromMatrix(countData=countDataMale, 
                              colData=MetaDataMale, 
                              design=~Site + Tissue + C9ORF72 + Age_Factor + RIN_Factor + Status)

dds$Status = relevel(factor(dds$Status),ref="control")
keep <- rowSums(counts(dds)) > 5 
dds.filt <- dds[keep,]

# Find DEGs
ddsObj <- DESeq(dds.filt)
ddsObj <- estimateSizeFactors(ddsObj)
results.dds <- results(ddsObj, alpha = 0.05)
results.dds.sig <- subset(results.dds,!is.na(results.dds$padj))
results.dds.sig.ALS <- subset(results.dds, results.dds$padj  < 0.05)
results.dds.sig.ALS.ENS <- results.dds.sig.ALS[grep("ENS", rownames(results.dds.sig.ALS)),]

write.csv(results.dds.sig.ALS, 
          "allctx.Males.CompModel.csv")

write.csv(results.dds.sig.ALS.ENS, 
          "allctx.Males.LFTab.ENSOnly.CompModel.csv")
results.dds.sig.ALS.ENS.IPA <- results.dds.sig.ALS.ENS %>% 
  as.data.frame() %>% 
  dplyr::arrange(padj) %>% 
  dplyr::slice_head(n = 1550)
write.csv(results.dds.sig.ALS.ENS.IPA, 
          "allctx.Males.ForIPA.LFTab.csv")

# Shrink vals for downstream analysis
library(ashr)
ddsShrink <- lfcShrink(ddsObj,
                       res = results.dds,
                       type = "ashr")

vsd = vst(ddsObj, blind = FALSE)
vstCounts <- assay(vsd)
write.csv(vstCounts, "DispNotBlind_VSTCounts_allctx_Males.csv")

jpeg("PCA_Males_CompModel_ntop500.jpg")
DESeq2::plotPCA(vsd,
                intgroup = "Status",
                ntop = 500,
                returnData = FALSE)
dev.off()


resLFC <- lfcShrink(ddsObj, coef="Status_als_vs_control", type="ashr")
jpeg("MAregularresults_Males_CompleteModel.jpg")
plotMA(results.dds, ylim=c(-1,1))
dev.off()

jpeg("MAlfcShrunk_Males_CompleteModel.jpg")
plotMA(resLFC, ylim=c(-1,1))
dev.off()


library("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
ens <- rownames(resLFC)

annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_transcript_id", "ensembl_gene_id",
               "gene_biotype", "external_gene_name"),
  filter="ensembl_gene_id",
  values=ens,
  uniqueRows=TRUE)

annotLookup2 <- annotLookup %>% 
  as.data.frame() %>% 
  distinct(ensembl_gene_id, .keep_all = TRUE) %>% 
  arrange(match(ensembl_gene_id, ens))
ens2 <- substr(ens,1,15)  
ens2[ens2 %in% annotLookup2$ensembl_gene_id] = annotLookup2$external_gene_name


library(EnhancedVolcano)
jpeg("VolcanoPlot_ShrunkFC_Males_CompleteModel.jpg")

EnhancedVolcano(resLFC,
                lab = ens2,
                x = 'log2FoldChange',
                title = 'Volcano Plot Male ALS v Control',
                y = 'pvalue',
                xlim = c(-5,5),
                colAlpha = 1)
dev.off()

ERVKFeats <- as.data.frame(results.dds.sig.ALS[grep("ERVK", rownames(results.dds.sig.ALS)),])
feat_of_interest <- grep("env", rownames(ERVKFeats),value = TRUE)

ERVKFeats <- dplyr::mutate(ERVKFeats,
                           axis_col = ifelse(rownames(ERVKFeats) %in% feat_of_interest, "env", "other"))

axis_cols = myPalette[ERVKFeats$axis_col]
axis_cols = axis_cols[order(rownames(ERVKFeats))]

write.csv(ERVKFeats, "AllCtxDEA_Males_FeatPlot_HERVKOnly.csv")
ggplot(ERVKFeats, aes(x = substr(rownames(ERVKFeats),1,10), y = log2FoldChange, col = padj)) +
  geom_point(size = 4) +
  labs(title="HERVK DE Features in All Cortical Male Samples",
       x="Gene", y = "Expression (log2FC)") +
  scale_color_gradient(low="black", high="orange") +
  ylim(-1.5,1.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, 
                                   colour = axis_cols),
        axis.text = element_text(size=11),
        plot.title = element_text(size=16))
ggsave("AllCtxDEA_Males_NAP_HERVKOnly.jpeg", 
       device="jpeg")