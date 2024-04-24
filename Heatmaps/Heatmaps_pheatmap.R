library(readxl)
library(dplyr)
library(tidyverse)
library(janitor)
library(pheatmap)
library(grid)

# Function based on: https://stackoverflow.com/questions/43051525/how-to-draw-pheatmap-plot-to-screen-and-also-save-to-file
save_pheatmap_pdf <- function(x, filename, width=12, height=8) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  setHook("grid.newpage", NULL, "replace")
  grid.text("Samples", y=-0.07, gp=gpar(fontsize=16))
  grid.text("MAD 5000 Features", x=-0.07, rot=90, gp=gpar(fontsize=16))
  dev.off()
}

## Read in data
MetaData = read.csv("meta_data.csv")

mad5000 = read.table("allRegions_MAD5000.txt")

mat = as.matrix(mad5000)



cat = data.frame(MetaData$Gender,MetaData$Class_2,MetaData$Tissue)
colnames(cat) = c("gender","status", "tissue")
rownames(cat) = MetaData$Sample

cat$tissue_type = MetaData$Tissue_subtype
mat2 = mat[,order(cat$tissue_type)]

# Only interested in CBM, CTX, and SC
CNS_regions <- c("cerebellum", "cortex", "spinal_cord")

cat <- subset(cat, cat$tissue_type %in% CNS_regions)
mat2 <- mat2[,colnames(mat2) %in% rownames(cat)]

# Center but do not scale matrix as it is VST
centered_mat <- t(apply(mat2, 1, function(x) scale(x, scale=FALSE)))
rownames(centered_mat) = rownames(mat2)
colnames(centered_mat) = colnames(mat2)


# Tissue Type

col_ann = data.frame(Tissue_Type = cat$tissue_type)
rownames(col_ann) = rownames(cat)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
phmap <- pheatmap(centered_mat,main = "CNS Region Tissue Type Heatmap",cluster_rows = FALSE, scale = "none",
                  annotation_col = col_ann, show_colnames = FALSE, 
                  show_rownames = FALSE)
setHook("grid.newpage", NULL, "replace")
grid.text("Samples", y=-0.07, gp=gpar(fontsize=16))
grid.text("MAD 5000 Features", x=-0.07, rot=90, gp=gpar(fontsize=16))
save_pheatmap_pdf(phmap, "Heatmap_Global_CNSRegion_AllMAD.pdf")


# Biological Sex

col_ann = data.frame(Biological_Sex = cat$gender)
rownames(col_ann) = rownames(cat)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
phmap <- pheatmap(centered_mat,main = "CNS Region Biological Sex Heatmap",cluster_rows = FALSE, scale = "none",
                  annotation_col = col_ann, show_colnames = FALSE, 
                  show_rownames = FALSE)
setHook("grid.newpage", NULL, "replace")
grid.text("Samples", y=-0.07, gp=gpar(fontsize=16))
grid.text("MAD 5000 Features", x=-0.07, rot=90, gp=gpar(fontsize=16))
save_pheatmap_pdf(phmap, "Heatmap_Global_Sex_AllMAD.pdf")

# Pt Status

col_ann = data.frame(Status = cat$status)
rownames(col_ann) = rownames(cat)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
phmap <- pheatmap(centered_mat,main = "CNS Region Patient Status Heatmap",cluster_rows = FALSE, scale = "none",
                  annotation_col = col_ann, show_colnames = FALSE, 
                  show_rownames = FALSE)
setHook("grid.newpage", NULL, "replace")
grid.text("Samples", y=-0.07, gp=gpar(fontsize=16))
grid.text("MAD 5000 Features", x=-0.07, rot=90, gp=gpar(fontsize=16))
save_pheatmap_pdf(phmap, "Heatmap_Global_Status_AllMAD.pdf")



##############################################################
# Just HERVK Now

mat2 = mat2[grepl("ERVK", rownames(mat2)),]
cat <- subset(cat, rownames(cat) %in% colnames(mat2))

centered_mat <- t(apply(mat2, 1, function(x) scale(x, scale=FALSE)))
rownames(centered_mat) = rownames(mat2)
colnames(centered_mat) = colnames(mat2)


# Tissue Type

col_ann = data.frame(Tissue_Type = cat$tissue_type)
rownames(col_ann) = rownames(cat)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
phmap <- pheatmap(centered_mat,main = "CNS Region Tissue Type Heatmap HERV-K Features Only",cluster_rows = FALSE, scale = "none",
                  annotation_col = col_ann, show_colnames = FALSE, 
                  show_rownames = FALSE)
setHook("grid.newpage", NULL, "replace")
grid.text("Samples", y=-0.07, gp=gpar(fontsize=16))
grid.text("MAD 5000 HERV-K Features", x=-0.07, rot=90, gp=gpar(fontsize=16))
save_pheatmap_pdf(phmap, "Heatmap_Global_CNSRegion_HERVKMAD.pdf")


# Biological Sex

col_ann = data.frame(Biological_Sex = cat$gender)
rownames(col_ann) = rownames(cat)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
phmap <- pheatmap(centered_mat,main = "CNS Region Biological Sex Heatmap HERV-K Features Only",cluster_rows = FALSE, scale = "none",
                  annotation_col = col_ann, show_colnames = FALSE, 
                  show_rownames = FALSE)
setHook("grid.newpage", NULL, "replace")
grid.text("Samples", y=-0.07, gp=gpar(fontsize=16))
grid.text("MAD 5000 HERV-K Features", x=-0.07, rot=90, gp=gpar(fontsize=16))
save_pheatmap_pdf(phmap, "Heatmap_Global_Sex_HERVKMAD.pdf")

# Pt Status

col_ann = data.frame(Status = cat$status)
rownames(col_ann) = rownames(cat)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
phmap <- pheatmap(centered_mat,main = "CNS Region Patient Status Heatmap HERV-K Features Only",cluster_rows = FALSE, scale = "none",
                  annotation_col = col_ann, show_colnames = FALSE, 
                  show_rownames = FALSE)
setHook("grid.newpage", NULL, "replace")
grid.text("Samples", y=-0.07, gp=gpar(fontsize=16))
grid.text("MAD 5000 HERV-K Features", x=-0.07, rot=90, gp=gpar(fontsize=16))
save_pheatmap_pdf(phmap, "Heatmap_Global_Status_HERVKMAD.pdf")



###################################################################
# All Features Cortex Only
################################################################


mad5000 = read.table("CTX_MAD5000.txt")

mat = as.matrix(mad5000)

cat = data.frame(MetaData$Gender,MetaData$Class_2,MetaData$Tissue)
colnames(cat) = c("gender","status", "tissue")
rownames(cat) = MetaData$Sample

cat$tissue_type = MetaData$Tissue_subtype

# Only interested in CTX
CNS_regions <- c("cortex")

cat <- subset(cat, cat$tissue_type %in% CNS_regions)
mat2 <- mat[,colnames(mat) %in% rownames(cat)]


centered_mat <- t(apply(mat2, 1, function(x) scale(x, scale=FALSE)))
rownames(centered_mat) = rownames(mat2)
colnames(centered_mat) = colnames(mat2)


# Tissue Type

col_ann = data.frame(Tissue_Type = cat$tissue)
rownames(col_ann) = rownames(cat)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
phmap <- pheatmap(centered_mat,main = "Cortex Region Tissue Type Heatmap",cluster_rows = FALSE, scale = "none",
                  annotation_col = col_ann, show_colnames = FALSE, 
                  show_rownames = FALSE)
setHook("grid.newpage", NULL, "replace")
grid.text("Samples", y=-0.07, gp=gpar(fontsize=16))
grid.text("MAD 5000 Features", x=-0.07, rot=90, gp=gpar(fontsize=16))
save_pheatmap_pdf(phmap, "Heatmap_CTX_CNSRegion_AllMAD.pdf")


# Sex

col_ann = data.frame(Biological_Sex = cat$gender)
rownames(col_ann) = rownames(cat)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
phmap <- pheatmap(centered_mat,main = "Cortex Region Biological Sex Heatmap",cluster_rows = FALSE, scale = "none",
                  annotation_col = col_ann, show_colnames = FALSE, 
                  show_rownames = FALSE)
setHook("grid.newpage", NULL, "replace")
grid.text("Samples", y=-0.07, gp=gpar(fontsize=16))
grid.text("MAD 5000 Features", x=-0.07, rot=90, gp=gpar(fontsize=16))
save_pheatmap_pdf(phmap, "Heatmap_CTX_Sex_AllMAD.pdf")

# Pt Status

col_ann = data.frame(Status = cat$status)
rownames(col_ann) = rownames(cat)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
phmap <- pheatmap(centered_mat,main = "Cortex Region Patient Status Heatmap",cluster_rows = FALSE, scale = "none",
                  annotation_col = col_ann, show_colnames = FALSE, 
                  show_rownames = FALSE)
setHook("grid.newpage", NULL, "replace")
grid.text("Samples", y=-0.07, gp=gpar(fontsize=16))
grid.text("MAD 5000 Features", x=-0.07, rot=90, gp=gpar(fontsize=16))
save_pheatmap_pdf(phmap, "Heatmap_CTX_Status_AllMAD.pdf")

