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
myPalette <- getColors(numColors)
names(myPalette) <- c("env", "other")

##############################################

MetaData2 <- read.csv("MetaData2.csv", row.names = 1)
countDataMatrix2 <- read.csv("CountDataMatrix2.csv", row.names = 1)
HERVKPtMeta <- read.csv("HERVKPts_MetaData.csv", row.names = 1)
HERVKCtrlMeta <- read.csv("HERVKCtrls_MetaData.csv", row.names = 1)
CtrlFreqCounts <- read.csv("FreqCounts_Ctrls.csv", row.names = 1)
PtFreqCounts <- read.csv("FreqCounts_Pts.csv", row.names = 1)

HERVKMediatedALSptsFinal <- HERVKPtMeta$Sample
HERVKMediatedCtrlsFinal <- HERVKCtrlMeta$Sample
# Meta Data Split
metaData = MetaData2
metaDataPts <- metaData %>% 
  dplyr::mutate(IsHERVK = if_else(metaData$Sample %in% HERVKMediatedALSptsFinal,
                                  "High_HERVK","Other")) %>% 
  dplyr::filter(Status == "als")  

metaDataCtrls <- metaData %>% 
  dplyr::mutate(IsHERVK = if_else(metaData$Sample %in% HERVKMediatedCtrlsFinal,
                                  "High_HERVK","Other")) %>% 
  dplyr::filter(Status == "control") 



############################################
# What's the best percent of HERVs for COntrols
########################################
# Proportion = 0.3
delta_m_list = list()
pt_median <- median(PtFreqCounts$freq)
for (s in seq(0.05,0.5, by = 0.05)) {
  temp1 <- CtrlFreqCounts %>% 
    dplyr::slice_max(freq,prop=s) %>% 
    dplyr::pull(freq)
  delta_m <- median(temp1) - pt_median
  p_val <- wilcox.test(temp1, PtFreqCounts$freq)$p.value
  s_id <- as.character(s)
  delta_m_list[[s_id]] <- c(delta_median = delta_m,p_value = p_val)
}

delta_m_df <- as.data.frame(dplyr::bind_rows(delta_m_list, .id = s_id))
colnames(delta_m_df)[1] <- "proportion"
delta_m_df$q_val <- p.adjust(delta_m_df$p_value,method = "fdr")


delta_m_df %>% 
  ggplot(aes(proportion,group=1)) +
  geom_point(aes(y = delta_median),size=6) +
  geom_line(aes(y = delta_median, colour = "Delta Medians"),size=2) +
  geom_point(aes(y = -log(q_val,base = 5)),size=6) +
  geom_line(aes(y = -log(q_val,base = 5), colour = "-log5(q-value)"),size=2) +
  labs(title="Delta Median Frequency Values Vs Proportion Controls", y = "") +
  theme_classic()
ggsave("OptimizationOfProp_FreqVal_Ctrls.jpeg",device="jpeg")




######################################################
# Compare HERV Ctrls to Pts
HervCtrls <- subset(MetaData2$Sample, MetaData2$Status == "control" & MetaData2$Sample %in% HERVKMediatedCtrlsFinal)
MetaData3 <- subset(MetaData2, MetaData2$Sample %in% HervCtrls | MetaData2$Status=="als")
countDataMatrix3 <- countDataMatrix2[,colnames(countDataMatrix2) %in% MetaData3$Sample]

dds <- DESeqDataSetFromMatrix(countData=countDataMatrix3, 
                              colData=MetaData3, 
                              design=~Site + Tissue + C9ORF72 + Age_Factor + RIN_Factor + Status)

dds$Status = relevel(factor(dds$Status),ref="als")
keep <- rowSums(counts(dds)) > 5 
dds.filt <- dds[keep,]

# Find DEGs
ddsObj <- DESeq(dds.filt)
ddsObj <- estimateSizeFactors(ddsObj)
results.dds <- results(ddsObj, alpha = 0.05)
results.dds.sig <- subset(results.dds,!is.na(results.dds$padj))
results.dds.sig.ALS <- subset(results.dds, results.dds$padj  < 0.05)
results.dds.sig.ALS.ENS <- results.dds.sig.ALS[grep("ENS", rownames(results.dds.sig.ALS)),]


results.dds.sig.ALS.ENS.IPA <- results.dds.sig.ALS.ENS %>% 
  as.data.frame() %>% 
  dplyr::arrange(padj) %>% 
  dplyr::slice_head(n = 1550)


write.csv(results.dds.sig.ALS, "HERVCtrlsvPts_CompModel_DEA.csv")
write.csv(results.dds.sig.ALS.ENS.IPA, "ForIPA_HERVCtrlsvPts_CompModel_DEA.csv")



# Shrink vals for downstream analysis
library(ashr)
ddsShrink <- lfcShrink(ddsObj,
                       res = results.dds,
                       type = "ashr")

vsd = vst(ddsObj, blind = FALSE)
vstCounts <- assay(vsd)



resLFC <- lfcShrink(ddsObj, coef="Status_control_vs_als", type="ashr")


# Volc Stuff
library(biomaRt)
mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "useast")
ens <- rownames(resLFC)[!grepl("ERV",rownames(resLFC))]

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
ens2 <- substr(rownames(resLFC),1,15)  
ens2[ens2 %in% annotLookup2$ensembl_gene_id] = annotLookup2$external_gene_name



library(EnhancedVolcano)
jpeg("VolcanoPlot_ShrunkFC_SexFilt_CompleteModel_HERVCtrlsvPts.jpg")

EnhancedVolcano(resLFC,
                lab = ens2,
                x = 'log2FoldChange',
                title = 'Volcano Plot HERVK Controls v All Patients',
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

write.csv(ERVKFeats, "AllCtx_DEACompMod_AllSamps_HERVCtrlsvPts_LFTab.csv")
ggplot(ERVKFeats, aes(x = substr(rownames(ERVKFeats),1,10), y = log2FoldChange, col = padj)) +
  geom_point(size = 4) +
  labs(title="HERVK DE Features in Cortical HERV Controls vs Pts",
       x="Gene", y = "Expression (log2FC)") +
  scale_color_gradient(low="black", high="orange") +
  ylim(-3,3) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, 
                                   colour = axis_cols),
        axis.text = element_text(size=11),
        plot.title = element_text(size=16))
ggsave("AllCtx_DEACompMod_AllSamps_HERVCtrlsvPts_Volc.jpeg", 
       device="jpeg")

######################################################
# Remove HERV Ctrls
NoHervCtrls <- subset(MetaData2$Sample, MetaData2$Status == "control" & !(MetaData2$Sample %in% HERVKMediatedCtrlsFinal))
MetaData3 <- subset(MetaData2, MetaData2$Sample %in% NoHervCtrls | MetaData2$Status=="als")
countDataMatrix3 <- countDataMatrix2[,colnames(countDataMatrix2) %in% MetaData3$Sample]

dds <- DESeqDataSetFromMatrix(countData=countDataMatrix3, 
                              colData=MetaData3, 
                              design=~Site + Tissue + C9ORF72 + Age_Factor + RIN_Factor + Status)

dds$Status = relevel(factor(dds$Status),ref="als")
keep <- rowSums(counts(dds)) > 5 
dds.filt <- dds[keep,]

# Find DEGs
ddsObj <- DESeq(dds.filt)
ddsObj <- estimateSizeFactors(ddsObj)
results.dds <- results(ddsObj, alpha = 0.05)
results.dds.sig <- subset(results.dds,!is.na(results.dds$padj))
results.dds.sig.ALS <- subset(results.dds, results.dds$padj  < 0.05)
results.dds.sig.ALS.ENS <- results.dds.sig.ALS[grep("ENS", rownames(results.dds.sig.ALS)),]


results.dds.sig.ALS.ENS.IPA <- results.dds.sig.ALS.ENS %>% 
  as.data.frame() %>% 
  dplyr::arrange(padj) %>% 
  dplyr::slice_head(n = 1550)


write.csv(results.dds.sig.ALS, "NoHERVCtrlsvPts_CompModel_DEA.csv")
write.csv(results.dds.sig.ALS.ENS.IPA, "ForIPA_NoHERVCtrlsvPts_CompModel_DEA.csv")



# Shrink vals for downstream analysis
library(ashr)
ddsShrink <- lfcShrink(ddsObj,
                       res = results.dds,
                       type = "ashr")

vsd = vst(ddsObj, blind = FALSE)
vstCounts <- assay(vsd)



resLFC <- lfcShrink(ddsObj, coef="Status_control_vs_als", type="ashr")


# Volc Stuff
library(biomaRt)
mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "useast")
ens <- rownames(resLFC)[!grepl("ERV",rownames(resLFC))]

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
ens2 <- substr(rownames(resLFC),1,15)  
ens2[ens2 %in% annotLookup2$ensembl_gene_id] = annotLookup2$external_gene_name



library(EnhancedVolcano)
jpeg("VolcanoPlot_ShrunkFC_SexFilt_CompleteModel_NoHERVCtrlsvPts.jpg")

EnhancedVolcano(resLFC,
                lab = ens2,
                x = 'log2FoldChange',
                title = 'Volcano Plot Non-HERVK Controls v All Pts',
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

write.csv(ERVKFeats, "AllCtx_DEACompMod_AllSamps_NoHERVCtrlsvPts_LFTab.csv")
ggplot(ERVKFeats, aes(x = substr(rownames(ERVKFeats),1,10), y = log2FoldChange, col = padj)) +
  geom_point(size = 4) +
  labs(title="HERVK DE Features in Cortical Non-HERV Controls vs Pts",
       x="Gene", y = "Expression (log2FC)") +
  scale_color_gradient(low="black", high="orange") +
  ylim(-3,3) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, 
                                   colour = axis_cols),
        axis.text = element_text(size=11),
        plot.title = element_text(size=16))
ggsave("AllCtx_DEACompMod_AllSamps_Non_HERVCtrlsvPts_Volc.jpeg", 
       device="jpeg")



######################################################
# Compare High HERV Ctrls to Other Ctrls
countDataMatrix3 <- countDataMatrix2[,colnames(countDataMatrix2) %in% metaDataCtrls$Sample]

dds <- DESeqDataSetFromMatrix(countData=countDataMatrix3, 
                              colData=metaDataCtrls, 
                              design=~Site + Tissue + C9ORF72 + Age_Factor + RIN_Factor + IsHERVK)

dds$IsHERVK = relevel(factor(dds$IsHERVK),ref="Other")
keep <- rowSums(counts(dds)) > 5 
dds.filt <- dds[keep,]

# Find DEGs
ddsObj <- DESeq(dds.filt)
ddsObj <- estimateSizeFactors(ddsObj)
results.dds <- results(ddsObj, alpha = 0.05)
results.dds.sig <- subset(results.dds,!is.na(results.dds$padj))
results.dds.sig.ALS <- subset(results.dds, results.dds$padj  < 0.05)
results.dds.sig.ALS.ENS <- results.dds.sig.ALS[grep("ENS", rownames(results.dds.sig.ALS)),]


results.dds.sig.ALS.ENS.IPA <- results.dds.sig.ALS.ENS %>% 
  as.data.frame() %>% 
  dplyr::arrange(padj) %>% 
  dplyr::slice_head(n = 1550)


write.csv(results.dds.sig.ALS, "HighHERVCtrlsvOtherCtrls_CompModel_DEA.csv")
write.csv(results.dds.sig.ALS.ENS.IPA, "ForIPA_HighHERVCtrlsvOtherCtrls_CompModel_DEA.csv")



# Shrink vals for downstream analysis
library(ashr)
ddsShrink <- lfcShrink(ddsObj,
                       res = results.dds,
                       type = "ashr")

vsd = vst(ddsObj, blind = FALSE)
vstCounts <- assay(vsd)



resLFC <- lfcShrink(ddsObj, coef="IsHERVK_High_HERVK_vs_Other", type="ashr")


# Volc Stuff
library(biomaRt)
mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "useast")
ens <- rownames(resLFC)[!grepl("ERV",rownames(resLFC))]

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
ens2 <- substr(rownames(resLFC),1,15)  
ens2[ens2 %in% annotLookup2$ensembl_gene_id] = annotLookup2$external_gene_name



library(EnhancedVolcano)
jpeg("VolcanoPlot_ShrunkFC_SexFilt_CompleteModel_HERVCtrlsvOtherCtrls.jpg")

EnhancedVolcano(resLFC,
                lab = ens2,
                x = 'log2FoldChange',
                title = 'Volcano Plot HERVK Controls v Other Controls',
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

write.csv(ERVKFeats, "AllCtx_DEACompMod_AllSamps_HERVCtrlsvOtherCtrls_LFTab.csv")
ggplot(ERVKFeats, aes(x = substr(rownames(ERVKFeats),1,10), y = log2FoldChange, col = padj)) +
  geom_point(size = 4) +
  labs(title="HERVK DE Features in Cortical HERV Controls vs Other Controls",
       x="Gene", y = "Expression (log2FC)") +
  scale_color_gradient(low="black", high="orange") +
  ylim(-3,3) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, 
                                   colour = axis_cols),
        axis.text = element_text(size=11),
        plot.title = element_text(size=16))
ggsave("AllCtx_DEACompMod_AllSamps_Non_HERVCtrlsvOtherCtrls_Volc.jpeg", 
       device="jpeg")



######################################################
# Compare HERV Ctrls to HERV Pts
HervCtrls <- subset(MetaData2$Sample, MetaData2$Status == "control" & MetaData2$Sample %in% HERVKMediatedCtrlsFinal)
HervPts <- subset(MetaData2$Sample, MetaData2$Status == "als" & MetaData2$Sample %in% HERVKMediatedALSptsFinal)
MetaData3 <- subset(MetaData2, MetaData2$Sample %in% HervCtrls | MetaData2$Sample %in% HervPts)
countDataMatrix3 <- countDataMatrix2[,colnames(countDataMatrix2) %in% MetaData3$Sample]

dds <- DESeqDataSetFromMatrix(countData=countDataMatrix3, 
                              colData=MetaData3, 
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


results.dds.sig.ALS.ENS.IPA <- results.dds.sig.ALS.ENS %>% 
  as.data.frame() %>% 
  dplyr::arrange(padj) %>% 
  dplyr::slice_head(n = 1550)


write.csv(results.dds.sig.ALS, "HERVCtrlsvHERVPts_CompModel_DEA.csv")
write.csv(results.dds.sig.ALS.ENS.IPA, "ForIPA_HERVCtrlsvHERVPts_CompModel_DEA.csv")



# Shrink vals for downstream analysis
library(ashr)
ddsShrink <- lfcShrink(ddsObj,
                       res = results.dds,
                       type = "ashr")

vsd = vst(ddsObj, blind = FALSE)
vstCounts <- assay(vsd)



resLFC <- lfcShrink(ddsObj, coef="Status_als_vs_control", type="ashr")


# Volc Stuff
library(biomaRt)
mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "useast")
ens <- rownames(resLFC)[!grepl("ERV",rownames(resLFC))]

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
ens2 <- substr(rownames(resLFC),1,15)  
ens2[ens2 %in% annotLookup2$ensembl_gene_id] = annotLookup2$external_gene_name



library(EnhancedVolcano)
jpeg("VolcanoPlot_ShrunkFC_SexFilt_CompleteModel_HERVCtrlsvHERVPts.jpg")

EnhancedVolcano(resLFC,
                lab = ens2,
                x = 'log2FoldChange',
                title = 'Volcano Plot HERVK Controls v HERVK Patients',
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

write.csv(ERVKFeats, "AllCtx_DEACompMod_AllSamps_HERVCtrlsvHERVPts_LFTab.csv")
ggplot(ERVKFeats, aes(x = substr(rownames(ERVKFeats),1,10), y = log2FoldChange, col = padj)) +
  geom_point(size = 4) +
  labs(title="HERVK DE Features in Cortical HERV Controls vs HERVK Pts",
       x="Gene", y = "Expression (log2FC)") +
  scale_color_gradient(low="black", high="orange") +
  ylim(-3,3) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, 
                                   colour = axis_cols),
        axis.text = element_text(size=11),
        plot.title = element_text(size=16))
ggsave("AllCtx_DEACompMod_AllSamps_HERVCtrlsvHERVPts_Volc.jpeg", 
       device="jpeg")
