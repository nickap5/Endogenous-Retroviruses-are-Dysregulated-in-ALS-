library(tidyverse)
library(biomaRt)
library(ggpubr)

# Get DEA Info
Ctx_DEA_Res <- read.csv("allctx.AllSamps.SexFilt.CompModel.csv")

Ctx_DEA_Res <- Ctx_DEA_Res %>% 
  dplyr::mutate(Ranks = abs(log2FoldChange)*(-log10(pvalue))) %>% 
  dplyr::mutate(Significance = (-log10(pvalue))) %>% 
  dplyr::arrange(desc(Ranks)) %>% 
  dplyr::rename(name = X)

# Read in MADs and add to DF
file_names <- list.files()
rmads <- read.csv(grep("Row_MADs",file_names,value = TRUE))
colnames(Ctx_DEA_Res)
colnames(rmads)[1] = "name"
Ctx_DEA_Res2 <- dplyr::inner_join(Ctx_DEA_Res,rmads,by="name")


##################################
# Correlations with Caret Var Imp 1 Min


# Read in Log Reg Lambda Min Data and add HGNCs
top_sig_lr_VImpMin <- read.csv(grep("CaretVarImp_Top_SiglambdaMin",file_names,value = TRUE))
colnames(top_sig_lr_VImpMin) = c("name", "Feature_Importance")

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- top_sig_lr_VImpMin$name

G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
colnames(G_list)[1] <- "name"
top_sig_lr_VImpMin_new <- top_sig_lr_VImpMin %>% 
  dplyr::left_join(G_list,by="name") %>% 
  dplyr::distinct(name, .keep_all = TRUE)

# Add rank to LR Coeffs

top_sig_lr_VImpMin_final <- full_join(top_sig_lr_VImpMin_new,Ctx_DEA_Res2,by="name")
top_sig_lr_VImpMin_final <- dplyr::filter(top_sig_lr_VImpMin_final, Feature_Importance != 0)


top_sig_lr_coeffs <- read.csv(grep("Lambda1SE.*Top_Sig",file_names,value = TRUE))
top_sig_lr_VImpMin_final2 <- right_join(top_sig_lr_VImpMin_new,top_sig_lr_coeffs,by="name")

sp <- ggscatter(top_sig_lr_VImpMin_final2, x = "coefficient", y = "Feature_Importance",
                add = "reg.line",  # Add regressin line
                xlim = c(0,0.2),
                ylim = c(0,0.2),
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + ggtitle("Spearman Correlation Between Logistic Regression Coefficient \nand Logistic Regression Feature Importance Lamda Min Cortex") +
  stat_cor(method = "spearman", label.x = 0.1, label.y = 0.2)
ggsave("Norm_Spearman_LRVImpMinvCoeff_CTX.jpeg")



sp <- ggscatter(top_sig_lr_VImpMin_final, x = "Ranks", y = "Feature_Importance",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + ggtitle("Spearman Correlation Between DEA Statistic \nand Logistic Regression Feature Importance Cortex") +
  stat_cor(method = "spearman", label.x = -25, label.y = 0.5) + xlab("|DEA Rank|") + ylab("LR Feature Importance")
ggsave("Spearman_LRVImpMinvDEARanks.jpeg")

sp <- ggscatter(top_sig_lr_VImpMin_final, x = "rmads", y = "Feature_Importance",
                add = "reg.line",  # Add regressin line
                xlim = c(0,2000),
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + ggtitle("Spearman Correlation Between MAD \nand Logistic Regression Feature Importance Cortex") +
  stat_cor(method = "spearman", label.x = 0, label.y = 0.5) + xlab("MAD") + ylab("LR Feature Importance")
ggsave("Spearman_LRVImpMinvMAD.jpeg")
