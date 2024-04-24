library(tidyverse)
library(DESeq2)
library(glmnet)
library(caret)
library(ROCR)
library(caTools)

set.seed(1) 

#########################################
# Function originally from:
# https://www.jihongzhang.org/post/2019-02-19-lasso-regression-with-glmnet/
##########################################
coeff2dt <- function(fitobject, s) {
  coeffs <- coef(fitobject, s) 
  coeffs.dt <- data.frame(name = coeffs@Dimnames[[1]][coeffs@i + 1], coefficient = coeffs@x) 
  
  # reorder the variables in term of coefficients
  return(coeffs.dt[order(coeffs.dt$coefficient, decreasing = T),])
}

########################################
# Read ENSG, TE and All
########################################

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


###########################################
# Sub Ctx

MetaData <- MetaData %>% 
  dplyr::filter(Tissue_subtype=="cortex") %>% 
  dplyr::rename(Status = Class_2)

countDataMatrixAll <- countDataMatrix[,colnames(countDataMatrix) %in% MetaData$Sample]

MetaData<-MetaData[MetaData$Sample %in% colnames(countDataMatrixAll),]
MetaData$C9ORF72 = as.factor(MetaData$C9ORF72)


######################################
# Find Most significant DEGs for Logit
DEA_LFTab <- read.csv("allctx.AllSamps.SexFilt.CompModel.csv")

ENSG_Sig_Feats <- DEA_LFTab %>% 
  dplyr::mutate(Ranks = abs(log2FoldChange)*(-log10(pvalue))) %>% 
  dplyr::arrange(desc(Ranks)) %>% 
  dplyr::filter(grepl("ENSG", X)) %>% 
  dplyr::slice_head(n = 500) %>% 
  dplyr::select(X)

ERV_Sig_Feats <- DEA_LFTab %>% 
  dplyr::mutate(Ranks = abs(log2FoldChange)*(-log10(pvalue))) %>% 
  dplyr::arrange(desc(Ranks)) %>% 
  dplyr::filter(grepl("ERV", X)) %>% 
  dplyr::slice_head(n = 500) %>% 
  dplyr::select(X)

All_Sig_Feats <- DEA_LFTab %>% 
  dplyr::mutate(Ranks = abs(log2FoldChange)*(-log10(pvalue))) %>% 
  dplyr::arrange(desc(Ranks)) %>% 
  dplyr::slice_head(n = 500) %>% 
  dplyr::select(X)

Bal_Sig_Feats <- rbind(head(ENSG_Sig_Feats,250),
                       head(ERV_Sig_Feats,250))

#Get MAD 500 feats
rmads <- apply(countDataMatrixAll, 1, mad) 
rmads <- sort(rmads,decreasing = TRUE)
MAD_Rows <- names(head(rmads,500))



# Balanced
temp1 <- rmads[grepl("ENSG",names(rmads))]
temp1 <- names(head(temp1,250))
temp2 <- rmads[grepl("ERV",names(rmads))]
temp2 <- names(head(temp2,250))
MAD_Bal_Rows <- c(temp1,temp2)


Rand_Feats <- sample(rownames(countDataMatrixAll),500)


target_row_n <- sum(grepl("ENSG", rownames(countDataMatrixAll)))
Rand_ERVs <- sample(grep("ERV",rownames(countDataMatrixAll), value = TRUE),target_row_n)

Bal_Rand_Feats <- subset(rownames(countDataMatrixAll),
                         grepl("ENSG", rownames(countDataMatrixAll)) | rownames(countDataMatrixAll) %in% Rand_ERVs)


##################
# Make List
list_of_feats = list()
list_of_feats[["Top_Sig"]] = All_Sig_Feats$X

list_of_feats[["Top_MAD"]] = MAD_Rows

list_of_feats[["Bal_Sig"]] = Bal_Sig_Feats$X

list_of_feats[["Bal_MAD"]] = MAD_Bal_Rows

list_of_feats[["Rand"]] = Rand_Feats

for (n in 1:10) {
  # n =1
  set.seed(n)
  name = paste0("Balanced_Random_",n)
  temp1 <- sample(grep("ERV",rownames(countDataMatrixAll), value = TRUE),250)
  temp2 <- sample(grep("ENSG",rownames(countDataMatrixAll), value = TRUE),250)
  temp3 <- c(temp1,temp2)
  
  list_of_feats[[name]] = temp3
}



########################################
# Logit Regression
########################################


df <- read.csv("MoR_NormCounts_AllCTX.csv",row.names = 1)

# Write MADs
df2 <- as.data.frame(rmads)
write.csv(df2, "Row_MADs_MoRNorm.csv")

for (n in names(list_of_feats)) {
  # Use right feats for matching
  
  row_names = as.character(list_of_feats[[n]])
  # Subset
  df2 <- dplyr::filter(df,rownames(df) %in% row_names)
  
  
  df2<-as.data.frame(t(df2))
  
  df2$Is_Pt = ifelse(MetaData$Class_2 == "control",0,1)
  
  df2$Is_Pt = factor(df2$Is_Pt)
  df2$Is_Pt = relevel(df2$Is_Pt,ref = "0")
  
  
  
  set.seed(1)
  
  df3 <- df2
  df2$Is_Pt<-NULL
  
  sample_split <- sample.split(Y = df3$Is_Pt, SplitRatio = 0.7)
  train_set <- subset(x = df2, sample_split == TRUE)
  train_set <- lapply(train_set, function(x) if(is.factor(x)) factor(x) else x)
  test_set <- subset(x = df2, sample_split == FALSE)
  test_set <- lapply(test_set, function(x) if(is.factor(x)) factor(x) else x)
  
  
  a <- list(rownames(df2)[sample_split==TRUE],colnames(df2))
  b <- list(rownames(df2)[sample_split==FALSE],colnames(df2))
  # Remember x is your predictors and y is your outcome
  x.train <- matrix(unlist(train_set), 
                    ncol = ncol(df2), 
                    nrow = sum(sample_split,na.rm = TRUE),
                    dimnames = a)
  x.test <- matrix(unlist(test_set), 
                   ncol = ncol(df2), 
                   nrow = length(df3$Is_Pt) - sum(sample_split,na.rm = TRUE),
                   dimnames = b)
  y.train <- subset(df3$Is_Pt,sample_split==TRUE)
  y.test <- subset(df3$Is_Pt,sample_split!=TRUE)
  
  cv_fit <- cv.glmnet(x.train, y.train, alpha = 0, nfolds = 5, 
                      type.measure="mse", family = "binomial")
  
  VImp <- caret::varImp(cv_fit$glmnet.fit, lambda = cv_fit$lambda.1se, scale = TRUE) 
  write.csv(VImp, paste0("Norm_CaretVarImp_",n,"lambda1se.csv"))
  VImp <- caret::varImp(cv_fit$glmnet.fit, lambda = cv_fit$lambda.min, scale = TRUE) 
  write.csv(VImp, paste0("Norm_CaretVarImp_",n,"lambdaMin.csv"))
  
  jpeg(paste0("Norm_CTX_CV_LogReg_GLMNet",n,".jpeg"))
  plot(cv_fit,label = TRUE)
  dev.off()
  
  write.csv(coeff2dt(fitobject = cv_fit, s = "lambda.min") %>% head(20),
            paste0("Norm_CTX_Coeffs_Lambda1SE_",as.character(cv_fit$lambda.min),"_",
                   n,".csv"))
  
  
  
  probs <- predict(cv_fit, newx = x.test, type = "response")
  
  pred <- factor(ifelse(probs > 0.7, "als", "control"))
  pred <- relevel(pred, ref = "control")
  
  cm <- confusionMatrix(factor(pred), 
                        factor(ifelse(y.test == "1", "als", "control")), 
                        positive = "als")
  
  write.csv(cm$byClass, paste0("Norm_CTX_Logit_Performance_03_",n,".csv"))
  cfm = cm$table
  
  save(file=paste0("Norm_CTX_GLM_Obj_03_",n),cv_fit)
  
  jpeg(paste0("Norm_CTX_FourFold_03_",n,".jpeg"))
  fourfoldplot(cfm, color = c("cyan", "pink"),
               conf.level = 0, margin = 1, main = "Confusion Matrix")
  dev.off()
  
  
  
  pred <- prediction(probs, y.test)
  
  # calculate probabilities for TPR/FPR for predictions
  jpeg(paste0("Norm_CTX_03_",n,"_ROCplot.jpeg"))
  perf <- performance(pred,"tpr","fpr")
  auc_ROCR <- performance(pred, measure = "auc")
  auc_ROCR <- round(auc_ROCR@y.values[[1]],3)
  plot(perf,colorize=FALSE, col="red", main = paste0("ROC Plot ",n)) # plot ROC curve
  lines(c(0,1),c(0,1),col = "black", lty = 1 )
  mtext(paste0("AUC = ",auc_ROCR), side=1, line=3.5, at=0.9)
  dev.off()
  
  
  
}

