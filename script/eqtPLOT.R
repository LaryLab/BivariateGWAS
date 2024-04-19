Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)
devtools::install_github("RitchieLab/eQTpLot")
library(eQTpLot)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("snpStats")
library(LDheatmap)
install.packages("LDheatmap")
source("https://bioconductor.org/biocLite.R")
biocLite(c("snpStats","rtracklayer","GenomicRanges","GenomInfoDb","IRanges"))
devtools::install_github("SFUStatgen/LDheatmap")



###exammple 

data(GWAS.df.example)
data(eQTL.df.example)
head(eQTL.df.example)
data(LD.df.example)
head(LD.df.example)


eQTpLot(GWAS.df = GWAS.df.example, eQTL.df = eQTL.df.example, gene = c("BBS1", "ACTN3"), 
        gbuild = "hg19",  trait = "LDL", tissue =  "all", CollapseMethod = "min", 
        GeneList = T) 

eQTpLot(GWAS.df = GWAS.df.example, eQTL.df = eQTL.df.example, gene = "BBS1", 
        gbuild = "hg19",  trait = "LDL", tissue =  "all", CollapseMethod = "min")
eQTpLot(GWAS.df = GWAS.df.example, eQTL.df = eQTL.df.example, gene = "BBS1", 
        gbuild = "hg19",  trait = "LDL", tissue =  "all", CollapseMethod = "min")


######## READ MY DATA 

df2<-read.table("/work/larylab/NAYEMA/BIVARIATE_GWAS/Bivariate_gwas/data/analyticfile_full_ELX_AD_FNBMD.txt",header=TRUE)
head(df2)
AD<- df2[, c("chromosome.outcome", "position.outcome","SNP", "pval.exposure","beta.exposure")]
FNBMD<-df2[, c("chromosome.outcome", "position.outcome","SNP","pval.outcome","beta.outcome")]

colnames(AD)<-c("CHR","BP","SNP","P","BETA")
colnames(FNBMD)<-c("CHR","BP","SNP","P","BETA")

AD_data <- AD %>%
  mutate(PHE = "Alzheimers")
FNBMD_data <- FNBMD %>%
  mutate(PHE = "FNBMD")

#### SELECT YOUR POSITION 

position <-124186714
# Range
range <- 500000 

# Calculate upstream and downstream positions
upstream <- position - range
downstream <- position + range 


AD_PLEKHA1ed_data[merged_data$chromosome.outcome ==17& 
                               merged_data$position.outcome >= upstream & 
                               merged_data$position.outcome <= downstream, ]