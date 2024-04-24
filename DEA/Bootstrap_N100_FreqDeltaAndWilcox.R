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
# Filter out sex-associated features since PCA separates by sex if you don't
cat = data.frame(MetaData2$Gender,MetaData2$Site,MetaData2$Tissue)
colnames(cat) = c("gender","group","tissue")
rownames(cat) = colnames(countDataMatrix2)

dds = DESeqDataSetFromMatrix(countData=countDataMatrix2,
                             colData=cat,design =~gender)
dds$gender = relevel(dds$gender,ref="Male")
dds = DESeq(dds,betaPrior=T)

## Obtain genes significantly different in female (i.e. gender associated) with FDR < 0.05
res = results(dds)
sig = res[! is.na(res$padj) & res$padj<0.05,]


countDataMatrix2 = countDataMatrix2[! (rownames(countDataMatrix2) %in% rownames(sig)),]


##################################################################
# Create a DESeq Data object from count matrix  looking at ALS v Ctrl

MetaData2$C9ORF72 = as.factor(MetaData2$C9ORF72)
dds <- DESeqDataSetFromMatrix(countData=countDataMatrix2, 
                              colData=MetaData2, 
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





##################################################################

# Meta Data Split
metaData = MetaData2
ALSPtsOnlyinfo <- metaData[metaData$Status=="als", ]
ALSPtsOnlyinfo <- tibble::rownames_to_column(ALSPtsOnlyinfo, "X")


# Ctrls

ControlsOnlyinfo <- metaData[metaData$Status=="control", ]
ControlsOnlyinfo <- tibble::rownames_to_column(ControlsOnlyinfo, "X")


# Get HERVKs
ERVKFeats <- as.data.frame(results.dds[grep("ERVK", rownames(results.dds)),])
# HERVK_List <- grep("env", rownames(ERVKFeats),value = TRUE)
HERVK_List <- rownames(ERVKFeats)



# VST Counts
vsd <- vst(ddsObj, blind=FALSE)

vst_counts <- assay(vsd)
write.csv(vst_counts, "VST_NotBlind_CompModel.csv")

# ALS Pts
metaData = MetaData2
ALSPtsOnlyinfo <- metaData[metaData$Status=="als", ]
ALSPtsOnlyinfo <- tibble::rownames_to_column(ALSPtsOnlyinfo, "X")
ALSPtsOnlyCounts <- counts(ddsObj, normalize=T)[, metaData$Status=="als"]

PtMean <- as.data.frame(rowMeans(ALSPtsOnlyCounts)) 
colnames(PtMean) = "Mean"

# Ctrls

ControlsOnlyinfo <- metaData[metaData$Status=="control", ]
ControlsOnlyinfo <- tibble::rownames_to_column(ControlsOnlyinfo, "X")
ControlsOnlyCounts <- counts(ddsObj, normalize=T)[, metaData$Status=="control"]

CtrlMean <- as.data.frame(rowMeans(ControlsOnlyCounts)) 
colnames(CtrlMean) = "Mean"
##################################################################
# Figure out  who the high HERV samples are

HERVK_DEA_Tab <- as.data.frame(results.dds) %>% 
  tibble::rownames_to_column("GeneID") %>% 
  filter(rownames(results.dds) %in% HERVK_List) %>%
  dplyr::rename(logFC = log2FoldChange) 


# Make a table of patients and genes above average
temp2 <- data.frame(matrix(ncol = 2, nrow = 0))
for (p in ALSPtsOnlyinfo$Sample) {
  
  # For every ALS Pt
  for (g in HERVK_List) {
    base_mean = abs(subset(HERVK_DEA_Tab$logFC,HERVK_DEA_Tab$GeneID == g))
    Fold_Change = 2^abs(subset(HERVK_DEA_Tab$logFC,HERVK_DEA_Tab$GeneID == g))
    FC_SE = 4*(2^subset(HERVK_DEA_Tab$lfcSE,HERVK_DEA_Tab$GeneID == g))
    # For all the ERV genes in our curated list
    if (abs(ALSPtsOnlyCounts[g,p]) > (base_mean*Fold_Change*FC_SE)){
      # if the expression value for the individual patient 
      # is greater than the threshold we are interested
      
      temp2 <- rbind(temp2, data.frame(patient = p, gene = g))
    }         
  }
}

HERVKPtCounts <- plyr::count(temp2$patient)
HERVKptGeneCounts <- plyr::count(temp2$gene)


# Make a table of controls and genes above average
temp2 <- data.frame(matrix(ncol = 2, nrow = 0))
for (p in ControlsOnlyinfo$Sample) {
  
  # For every control Pt
  for (g in HERVK_List) {
    base_mean = abs(subset(HERVK_DEA_Tab$logFC,HERVK_DEA_Tab$GeneID == g))
    Fold_Change = 2^abs(subset(HERVK_DEA_Tab$logFC,HERVK_DEA_Tab$GeneID == g))
    FC_SE = 4*(2^subset(HERVK_DEA_Tab$lfcSE,HERVK_DEA_Tab$GeneID == g))
    
    # For all the ERV genes in our curated list
    if (abs(ControlsOnlyCounts[g,p]) > (base_mean*Fold_Change*FC_SE)){
      # if the expression value for the individual patient 
      # is greater than the threshold we are interested
      
      temp2 <- rbind(temp2, data.frame(patient = p, gene = g))
    }
  }
}
ControlHERVKPtCount <- plyr::count(temp2$patient)
ControlHERVKGenes <- plyr::count(temp2$gene)

#########################################
# Bootstrap

ggdat <- rbind(data.frame(value = ControlHERVKPtCount$freq, variable = "Ctrl", 
                          Id = ControlHERVKPtCount$x), 
               data.frame(value = HERVKPtCounts$freq, variable = "Pt",
                          Id = HERVKPtCounts$x))

ggdat2 <- ggdat %>% 
  dplyr::arrange(Id) %>% 
  dplyr::mutate(Region = metaData$Tissue[order(metaData$Sample)])

boot_list = list()
for (n in 1:100) {
  temp7 <- ggdat2 %>% 
    dplyr::group_by(variable) %>% 
    dplyr::slice_sample(n=200)
  Ctrl_Freq = temp7$value[temp7$variable=="Ctrl"]
  Pt_Freq = temp7$value[temp7$variable=="Pt"]
  delta_med = (median(Ctrl_Freq) - median(Pt_Freq))
  wilcox_pval = wilcox.test(Ctrl_Freq,Pt_Freq)$p.value
  
  boot_list[[n]] = c(delta_med,wilcox_pval)
}
names(boot_list) = 1:100

boot_df <- as.data.frame(t(dplyr::bind_rows(boot_list)))
colnames(boot_df) <- c("delta_median","Wilcoxon_pval")
boot_df$Wilcoxon_qval <- p.adjust(boot_df$Wilcoxon_pval,method = "fdr")

write.csv(boot_df,"Bootstrapped_FreqVal_CtrlvPt.csv")
# None are Sig
table(boot_df$Wilcoxon_qval<0.05)

boot_df$significance = ifelse(boot_df$Wilcoxon_qval<0.05,"significant","not significant")
ggplot(boot_df, aes(x=as.numeric(rownames(boot_df)),y=delta_median)) +
  geom_point(size = 3,aes(colour = significance)) +
  geom_hline(yintercept = 2.5,colour="red") +
  ylim(-3,7) +
  labs(title = "Difference in Median Frequency Value Between ALS and Control", 
       x = "Iteration Number", y = "Median Controls - ALS") +
  theme_bw()

ggsave("BootstrapPlot.jpeg")


  
  
  
  