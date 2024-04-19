
##################################

library(devtools) 

######### Install coloc and dependencies

if(!require("remotes"))
  install.packages("remotes") # if necessary
library(remotes)
install_github("chr1swallace/coloc@main",build_vignettes=TRUE)
vignette("a06_SuSiE",package="coloc")
############### loading the birary
library(coloc)
library(snpStats)

if(requireNamespace("susieR",quietly=TRUE)) {
  susie.res=coloc.susie(S3,S4)
  print(susie.res$summary)
}
install.packages("susieR")

######### install gwasglue 
suppressPackageStartupMessages(suppressWarnings({
  library(gwasglue)
  library(dplyr)
  library(gassocplot)
  library(coloc)
}))

install.packages("LDlinkR")

remotes::install_github("CBIIT/LDlinkR")

library(LDlinkR)

###########test with gassocplot 

### assoc_plot  
markers <- gassocplot::test_assoc_plot  
head(markers)  
corr <- gassocplot::test_corr # this is correlation not correlation squared and has to be ordered in the same way as the markers data frame  
plot <- assoc_plot(markers, corr)   
assoc_plot_save(plot, "assoc_plot_test.png")  

### stack_assoc_plot  
markers <- gassocplot::test_stack_assoc_plot_markers  
head(markers)  
z <- gassocplot::test_stack_assoc_plot_associations  
head(z)  
corr <- gassocplot::test_corr # this is correlation not correlation squared and has to be ordered in the same way as the markers data frame  
plot <- stack_assoc_plot(markers, z, corr, traits=c("Trait 1", "Trait 2"))  
stack_assoc_plot_save(plot, "stack_assoc_plot_test.png", 2)



################# example glasscoplot with my data 

corr2 <-LDmatrix(snps = snp_list, 
                 pop = "EUR", 
                 token = "bf8048960eda",
                 genome_build = "grch37")
str(corr)
str(corr2)
corr2_snps <- corr2$RS_number

# Check for SNPs in snp_list that are not present in corr2
missing_snps <- setdiff(snp_list, corr2_snps)

corr2_values <- as.matrix(corr2[, -1])

# Assign row and column names
rownames(corr2_values) <- corr2$RS_number
colnames(corr2_values) <- corr2$RS_number

# Check the structure of the new matrix
str(corr2_values)
corr2_values_sqrt <- sqrt(corr2_values)

MARKER_filtered <- MARKER[!(MARKER$marker %in% c("rs10210125", "rs6722104","rs79845707")), ]
Z
# Check the head of the filtered dataset
head(MARKER_filtered)


marker_names <- MARKER_filtered$marker

# Reordering columns of corr2_values_sqrt based on marker_names
corr2_values_ordered <- corr2_values_sqrt[, marker_names]

# Displaying the ordered correlation matrix
print(corr2_values_ordered)

Z_filtered <- Z[!(rownames(Z) %in% c("rs12032780", "rs17186848")), ]

### stack_assoc_plot  
markers <-   MARKER_filtered
z <- Z_filtered 
corr <- corr2_values_ordered # this is correlation not correlation squared and has to be ordered in the same way as the markers data frame  
plot <- stack_assoc_plot(markers, z, corr, traits=c("Trait 1", "Trait 2"))  
stack_assoc_plot_save(plot, "stack_assoc_plot_test.png", 2)

######################################################





######### read your data
# my data

setwd("/work/larylab/NAYEMA/BIVARIATE_GWAS/Bivariate_gwas/")
df<-read.table("/work/larylab/NAYEMA/BIVARIATE_GWAS/Bivariate_gwas/result/ELX_AD_FNBMD1_result",header=TRUE)

df2<-read.table("/work/larylab/NAYEMA/BIVARIATE_GWAS/Bivariate_gwas/data/analyticfile_full_ELX_AD_FNBMD.txt",header=TRUE)
head(df)
head(df2)

merged_data <- merge(df, df2, by = "SNP")
head(merged_data)


######################################### colocalization tests based on hit position using coloc

# Given position
position <-2053920
# Range
range <- 500000 

# Calculate upstream and downstream positions
upstream <- position - range
downstream <- position + range 


filtered_data <- merged_data[merged_data$chromosome.outcome ==17& 
                               merged_data$position.outcome >= upstream & 
                               merged_data$position.outcome <= downstream, ]
head(filtered_data)

snp_list <- filtered_data$SNP

# View the first few SNPs
head(snp_list)

##### GET LD per pairwise SNPs 

corr2 <-LDmatrix(snps = snp_list, 
                 pop = "EUR", 
                 token = "bf8048960eda",
                 genome_build = "grch37")

str(corr2)
corr2_snps <- corr2$RS_number

# Check for SNPs in snp_list that are not present in corr2
missing_snps <- setdiff(snp_list, corr2_snps)

########## remove those missing SNPs from filtered data 

snps_to_remove <- c("rs10210125", "rs6722104", "rs79845707")

# Remove rows with specified SNP IDs
filtered_data <- subset(filtered_data, !(SNP %in% snps_to_remove))

# View the updated filtered_data
head(filtered_data)


#########CONVERT AS Matrix 
rownames(corr2) <- corr2$RS_number
corr2$RS_number <- NULL  # Remove RS_number column

# Convert data frame to matrix
corr_matrix <- as.matrix(corr2)

str(corr_matrix)


corr_matrix <- corr_matrix[match(filtered_data$SNP, rownames(corr_matrix)), 
                           match(filtered_data$SNP, rownames(corr_matrix))]




####################
coloc

#######################

position <-43898887
# Range
range <- 500000 

# Calculate upstream and downstream positions
upstream <- position - range
downstream <- position + range 


filtered_data <- merged_data[merged_data$chromosome.outcome ==17& 
                               merged_data$position.outcome >= upstream & 
                               merged_data$position.outcome <= downstream, ]
head(filtered_data)


AD_data <- filtered_data[, c("SNP", "chromosome.outcome", "position.outcome","effect_allele.exposure","other_allele.exposure","beta.exposure","se.exposure","pval.exposure","eaf.exposure")]
FNBMD_data<-filtered_data[, c("SNP", "chromosome.outcome", "position.outcome","effect_allele.exposure","other_allele.exposure","beta.outcome","se.outcome","pval.outcome","eaf.outcome")]





AD_data2 <- AD_data %>%
  mutate(N = 1404,
         sdY = 1.0,
         type = "cc")
head(AD_data2)

FNBMD_data2 <- FNBMD_data %>%
  mutate(N = 1404,
         sdY = 1.0,
         type = "quant")
head(FNBMD_data2)

colnames(FNBMD_data2)
names(AD_data2) <- c("snp", "chr", "position", "ea", "nea", "beta", "se", "pvalues", "MAF","N","sdY","type")
names(FNBMD_data2) <- c("snp", "chr", "position", "ea", "nea", "beta", "se", "pvalues", "MAF","N","sdY","type")



AD_data3 <- list(
  beta = AD_data2$beta,
  varbeta=AD_data2$se,
  N=length(AD_data2$beta),  # Assuming 'N' is the number of observations
  sdY = 1,  # Assuming 'sdY' has a constant value of 1
  type = "quant",  # Assuming 'type' has a constant value of "quant"
  MAF = AD_data2$MAF,
  snp = AD_data2$snp,
  position = AD_data2$position,
  chr=AD_data2$chr,
  ea=AD_data2$ea,
  nea=AD_data2$nea,
  se=AD_data2$se,
  pvalues=AD_data2$pvalues
)

str(D1)
check_dataset(AD_data3,warn.minp=1e-2)

FNBMD_data2 <- FNBMD_data2[FNBMD_data2$MAF > 0 & FNBMD_data2$MAF < 1, ]

FNBMD_data3 <- list(
  beta = FNBMD_data2$beta,
  varbeta=FNBMD_data2$se,
  N = length(FNBMD_data2$beta),  # Assuming 'N' is the number of observations
  sdY = 1,  # Assuming 'sdY' has a constant value of 1
  type = "quant",  # Assuming 'type' has a constant value of "quant"
  MAF = FNBMD_data2$MAF,
  snp = FNBMD_data2$snp,
  position = FNBMD_data2$position,
  chr=FNBMD_data2$chr,
  ea=FNBMD_data2$ea,
  nea=FNBMD_data2$nea,
  se=FNBMD_data2$se,
  pvalues=FNBMD_data2$pvalues
)

str(AD_data3)
check_dataset(FNBMD_data3,warn.minp=1e-2)

#AD_data3$LD <- corr_matrix
#FNBMD_data3$LD<- corr_matrix
str(FNBMD_data3)

plot_dataset(FNBMD_data3)

plot_dataset(AD_data3)
my.res <- coloc.abf(dataset1=AD_data3,
                    dataset2=FNBMD_data3)
#res <- coloc.signals(AD_data3,FNBMD_data3,method="cond")


#res <- coloc.signals(AD_data3,FNBMD_data3,method="cond",p12=1e-6)
#print(res) 


o <- order(my.res$results$SNP.PP.H4,decreasing=TRUE)
cs <- cumsum(my.res$results$SNP.PP.H4[o])
w <- which(cs > 0.95)[1]
my.res$results[o,][1:w,]$snp


subset(my.res$results,SNP.PP.H4= 0.00422)



write.csv(my.res$results, file = "/work/larylab/NAYEMA/colocalization/results_for_bivariategwas_sig/chr17_0.42%", row.names = FALSE, quote = FALSE)


###########################

library(coloc)
data(coloc_test_data)
attach(coloc_test_data)



# Find the index of the SNP "rs559133256" in the snp column of AD_data3
snp_index <- which(AD_data3$snp == "rs559133256")

# Extract beta and varbeta values for the SNP
beta_value <- AD_data3$beta[snp_index]
varbeta_value <- AD_data3$varbeta[snp_index]
