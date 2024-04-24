library(epitools)

# Check HERV-K feature association 

# Replace with relevant values
dd <- matrix (c(10,20,30,40),
              nrow=2, ncol=2, byrow=TRUE)
outc <- c("herv_als", "other_als")
herv <- c("down", "up")
dimnames(dd) <- list("HERV feats" = herv, "Outcome" = outc)

epitools::oddsratio.fisher(dd, rev = "c")

dd <- matrix (c(10,20,30,40),
              nrow=2, ncol=2, byrow=TRUE)
outc <- c("herv_ctrl", "other_als")
herv <- c("down", "up")
dimnames(dd) <- list("HERV feats" = herv, "Outcome" = outc)

epitools::oddsratio.fisher(dd, rev = "c")


# env


dd <- matrix (c(10,20,30,40),
              nrow=2, ncol=2, byrow=TRUE)
outc <- c("herv_als", "other_als")
herv <- c("down", "up")
dimnames(dd) <- list("env feats" = herv, "Outcome" = outc)
epitools::oddsratio.fisher(dd, rev = "c")


dd <- matrix (c(10,20,30,40),
              nrow=2, ncol=2, byrow=TRUE)
outc <- c("herv_als", "other_als")
herv <- c("down", "up")
dimnames(dd) <- list("env feats" = herv, "Outcome" = outc)
epitools::oddsratio.fisher(dd, rev = "c")

# Check chromosomal location associations
dd <- matrix (c(10,20,30,40),
              nrow=2, ncol=2, byrow=TRUE)
outc <- c("herv_als", "other_als")
chrom <- c("notCh3", "Ch3")
dimnames(dd) <- list("Chromosome 3" = chrom, "Outcome" = outc)

epitools::oddsratio.fisher(dd, rev = "c")


dd <- matrix (c(10,20,30,40),
              nrow=2, ncol=2, byrow=TRUE)
outc <- c("herv_control", "other_als")
chrom <- c("notCh3", "Ch3")
dimnames(dd) <- list("Chromosome 3" = chrom, "Outcome" = outc)

epitools::oddsratio.fisher(dd, rev = "c")



dd <- matrix (c(10,20,30,40),
              nrow=2, ncol=2, byrow=TRUE)
outc <- c("herv_als", "other_als")
chrom <- c("notCh8", "Ch8")
dimnames(dd) <- list("Chromosome 8" = chrom, "Outcome" = outc)

epitools::oddsratio.fisher(dd, rev = "c")

# Adjust p-values for multiple comparisons
# Replace with relevant values
p.adjust(c(0.01,0.02,0.03,0.04),method = "fdr")