library(tidyverse)
library(epitools)

#######################################
# Concordance of High HERV Patients
#######################################


########
# Ctx Data
###########

Ctx_Freq_Pts_Vals <- read.csv("FreqCounts_Pts.csv")
Ctx_Freq_Ctrls_Vals <- read.csv("FreqCounts_Ctrls.csv")
Ctx_Freq_Pts_Samps <- Ctx_Freq_Pts_Vals %>% 
  dplyr::slice_max(freq,prop=0.20) %>% 
  dplyr::pull(x)


################
# Meta Data
###############

MetaData <- readxl::read_xlsx("ALS_Consortium_Metadata.xlsx")
MetaData$ExternalSampleId = gsub("-","_",MetaData$ExternalSampleId)

Ctx_HERV_Meta <- MetaData %>% 
  dplyr::filter(ExternalSampleId %in% Ctx_Freq_Pts_Samps) %>% 
  dplyr::distinct(ExternalSubjectId,.keep_all = TRUE) 

Ctx_OtherPt_Meta <- MetaData %>% 
  dplyr::filter(grepl("ALS Spectrum MND", `Subject Group`)) %>% 
  dplyr::filter(grepl("Cortex", `Sample Source`)) %>% 
  dplyr::filter(!(ExternalSampleId %in% Ctx_Freq_Pts_Samps)) %>% 
  dplyr::distinct(ExternalSubjectId,.keep_all = TRUE) 


###################################################
##################################################
# Ctx ORs
###################################################
###################################################
# Replace with relevant values from meta data for all variables of interest


# Sex Female
# table(Ctx_HERV_Meta$Sex=="Male")
# table(Ctx_HERV_Meta$Sex=="Female")
dd <- matrix (c(10,20,30,40),
              nrow=2, ncol=2, byrow=TRUE)
outc <- c("PoI", "Other Pt")
meas <- c("Male", "Female")
dimnames(dd) <- list("Sex" = meas, "Type of Pt" = outc)
epitools::oddsratio.fisher(dd, rev = "c")





############################################
# Age at Death
##########################################
# quantile(as.numeric(MetaDataPoI$`Age at Death`))
table(as.numeric(Ctx_HERV_Meta$`Age at Death`) <= 55)
table(as.numeric(Ctx_HERV_Meta$`Age at Death`) > 65 & as.numeric(Ctx_HERV_Meta$`Age at Death`) < 70)
table(as.numeric(Ctx_HERV_Meta$`Age at Death`) >= 75)

table(as.numeric(Ctx_OtherPt_Meta$`Age at Death`) <= 55)
table(as.numeric(Ctx_OtherPt_Meta$`Age at Death`) > 65 & as.numeric(Ctx_OtherPt_Meta$`Age at Death`) < 70)
table(as.numeric(Ctx_OtherPt_Meta$`Age at Death`) >= 75)


# Young

dd <- matrix (c(10,20,30,40),
              nrow=2, ncol=2, byrow=TRUE)
outc <- c("PoI", "Other Pt")
meas <- c("No", "Yes")
dimnames(dd) <- list("Young" = meas, "Type of Pt" = outc)
epitools::oddsratio.fisher(dd, rev = "c")


# Middle

dd <- matrix (c(10,20,30,40),
              nrow=2, ncol=2, byrow=TRUE)
outc <- c("PoI", "Other Pt")
meas <- c("No", "Yes")
dimnames(dd) <- list("Middle" = meas, "Type of Pt" = outc)
epitools::oddsratio.fisher(dd, rev = "c")


# Old

dd <- matrix (c(10,20,30,40),
              nrow=2, ncol=2, byrow=TRUE)
outc <- c("PoI", "Other Pt")
meas <- c("No", "Yes")
dimnames(dd) <- list("Old" = meas, "Type of Pt" = outc)
epitools::oddsratio.fisher(dd, rev = "c")




############################################
# Age at Sx Onset
##########################################
# quantile(as.numeric(MetaDataPoI$`Age at Symptom Onset`), na.rm = TRUE)
table(as.numeric(Ctx_HERV_Meta$`Age at Symptom Onset`) <= 55)
table(as.numeric(Ctx_HERV_Meta$`Age at Symptom Onset`) > 65 & as.numeric(Ctx_HERV_Meta$`Age at Symptom Onset`) < 66)
table(as.numeric(Ctx_HERV_Meta$`Age at Symptom Onset`) >= 75)


table(as.numeric(Ctx_OtherPt_Meta$`Age at Symptom Onset`) <= 55)
table(as.numeric(Ctx_OtherPt_Meta$`Age at Symptom Onset`) > 65 & as.numeric(Ctx_OtherPt_Meta$`Age at Symptom Onset`) < 66)
table(as.numeric(Ctx_OtherPt_Meta$`Age at Symptom Onset`) >= 75)



# Young

dd <- matrix (c(10,20,30,40),
              nrow=2, ncol=2, byrow=TRUE)
outc <- c("PoI", "Other Pt")
meas <- c("No", "Yes")
dimnames(dd) <- list("Young" = meas, "Type of Pt" = outc)
epitools::oddsratio.fisher(dd, rev = "c")


# Middle

dd <- matrix (c(10,20,30,40),
              nrow=2, ncol=2, byrow=TRUE)
outc <- c("PoI", "Other Pt")
meas <- c("No", "Yes")
dimnames(dd) <- list("Middle" = meas, "Type of Pt" = outc)
epitools::oddsratio.fisher(dd, rev = "c")


# Old

dd <- matrix (c(10,20,30,40),
              nrow=2, ncol=2, byrow=TRUE)
outc <- c("PoI", "Other Pt")
meas <- c("No", "Yes")
dimnames(dd) <- list("Old" = meas, "Type of Pt" = outc)
epitools::oddsratio.fisher(dd, rev = "c")

# Replace with relevant values to correct for multiple comparisons
ps <- c(0.01,0.02,0.03)
p.adjust(ps,method = "fdr")
