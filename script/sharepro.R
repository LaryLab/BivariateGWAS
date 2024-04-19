######### read your data
# my data

setwd("/work/larylab/NAYEMA/BIVARIATE_GWAS/Bivariate_gwas/")
df<-read.table("/work/larylab/NAYEMA/BIVARIATE_GWAS/Bivariate_gwas/result/ELX_AD_FNBMD1_result",header=TRUE)

df2<-read.table("/work/larylab/NAYEMA/BIVARIATE_GWAS/Bivariate_gwas/data/analyticfile_full_ELX_AD_FNBMD.txt",header=TRUE)
head(df)
head(df2)

data_ad<-read.table("/work/larylab/NAYEMA/BIVARIATE_GWAS/AD.txt",fill = TRUE, header = TRUE)
head(data_ad)
data_ad2<-data_ad[,c("variant_id","chr","pos","effect_allele","other_allele","n_cases","n_controls")]
colnames(data_ad2)<-c("SNP","CHR","POS","A1","A2","N_CASE","N_CONTROL")
library(dplyr)

# Assuming merged_data is the merged dataset and data_ad2 is the dataset containing N_CASE and N_CONTROL

# Add N column as the sum of N_CASE and N_CONTROL
data_ad2 <- mutate(data_ad2, N = data_ad2$N_CASE + data_ad2$N_CONTROL)





merged_data <- merge(df, df2, by = "SNP")
head(merged_data)

AD_data <- merged_data[, c("SNP", "chromosome.outcome", "position.outcome","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure","se.exposure","pval.exposure")]
FNBMD_data<-merged_data[, c("SNP", "chromosome.outcome", "position.outcome","effect_allele.exposure","other_allele.exposure","eaf.outcome","beta.outcome","se.outcome","pval.outcome")]

colnames(AD_data)<-c("SNP","CHR","POS","A1","A2","EAF","BETA","SE","P")
colnames(FNBMD_data)<-c("SNP","CHR","POS","A1","A2","EAF","BETA","SE","P")

head(AD_data)
head(FNBMD_data)

FNBMD_data <- FNBMD_data %>%
  mutate(N = 32965)

merged_data_ad <- merge(AD_data, data_ad2, by = c("SNP", "CHR", "POS", "A1", "A2"))
head(merged_data_ad)
AD_data<-merged_data_ad[,c("SNP","CHR","POS","A1","A2","EAF","BETA","SE","P","N")]



library(LDlinkR)

# Given position
position <-1053524
# Range
range <- 500000 

# Calculate upstream and downstream positions
upstream <- position - range
downstream <- position + range


AD_filtered_data <- AD_data[AD_data$CHR ==19& 
                               AD_data$POS >= upstream & 
                              AD_data$POS <= downstream, ]

FNBMD_filtered_data <- FNBMD_data[FNBMD_data$CHR ==19& 
                                             FNBMD_data $POS >= upstream & 
                                             FNBMD_data $POS <= downstream, ]

snp_list <- AD_filtered_data$SNP

##### GET LD per pairwise SNPs 

corr2 <-LDmatrix(snps = snp_list, 
                 pop = "EUR", 
                 token = "2c8724750485",
                 genome_build = "grch37")


#corr2[is.na(corr2)] <- 0.001



corr2_snps <- corr2$RS_number

# Check for SNPs in snp_list that are not present in corr2
missing_snps <- setdiff(snp_list, corr2_snps)
AD_filtered_data <- AD_filtered_data[!(AD_filtered_data$SNP %in% missing_snps), ]
FNBMD_filtered_data <- FNBMD_filtered_data[!(FNBMD_filtered_data$SNP %in% missing_snps), ]

# Remove row names
rownames(corr2) <- NULL

# Remove column names
colnames(corr2) <- NULL
corr3<- corr2[, -1]




write.table(AD_filtered_data, file = "/work/larylab/NAYEMA/Sharepro/SharePro_coloc/AD_chr19_3.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(FNBMD_filtered_data, file = "/work/larylab/NAYEMA/Sharepro/SharePro_coloc/FNBMD_chr19_3.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(corr3, file = "/work/larylab/NAYEMA/Sharepro/SharePro_coloc/Chr19_3.ld", sep = "\t", quote = FALSE, row.names = FALSE)

