
##### Date: 27/02/2024
#### Author: Nayema
#### Merging file containing rsid from dbsnp with bmd summary statistic files that does not have rsid


##################################################
#1. Set environment
###############################
getwd()
setwd("/work/larylab/NAYEMA/BIVARIATE_GWAS/")


######################################################################################
#2. Load Packages
######################################################################################
#First make sure you have the biomaRt and devtools packages installed:

#Load packages
#install.packages("devtools")
library(devtools)
#install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("knitr")
library(knitr)
#install_github("MRCIEU/MRInstruments")
#The current results were achieved with using version TwoSampleMR_0.4.11
library(MRInstruments)
#install.packages("reshape")
library("reshape")
library(dplyr)


#####################
#3. read the summary statistics and the matched file from dbsnp for getting RSID 
#####################

#### femoral neck -fnbmd summary statistics data for european people downloaded from musculoskeleta knowledge portal 
#### link http://www.gefos.org/?q=content/data-release-2015
## from dbsnp data where the values in the CHROM and POS columns match with the corresponding values in fn2stu.MAF0_.005.pos_.out_, fa2stu.MAF0_.005.pos_.out_, ls2stu.MAF0_.005.pos.out_

############ fnbmd

fn2stu <- read.table("fn2stu.MAF0_.005.pos.out_", header = TRUE)
matched_lines <- read.table("matched_lines_fnbmd_no_indel.txt", header=FALSE) # lines of fnbmd based on chr and pos with dbsnp build 37 data

#check the data
head(fn2stu)
head(matched_lines)

####### give them column names

colnames(matched_lines) <- c("chromosome","position","snp","other_allele", "reference_allele",".", "..1", "..2")
matched_lines <- matched_lines[, 1:5]# i need only these 5 columns

# Merge matched line with fnbmd based on matching values in specified columns to get rsid

merged_datafn <- merge(matched_lines, fn2stu, by = c("chromosome", "position"))
head(merged_datafn)

##### extract only required columns 

merged_datafn2 <- merged_datafn[, c(1:3, 6:11, 15)]
head(merged_datafn2)

#change column names
colnames(merged_datafn2) <- c("chromosome","position","variant_id", "rs_number","reference_allele.fnbmd","other_allele.fnbmd","eaf.fnbmd","beta.fnbmd","se.fnbmd","p.value.fnbmd")

######### save the fnbmd data 
# Write the data frame to a file without row names
write.table(merged_datafn2, file = "/work/larylab/NAYEMA/BIVARIATE_GWAS/Bivariate_gwas/data/fnbmd_data_with_rsid.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#############################
############ Alzheimers data , it was converted from build 38 to 37 using liftover

data_ad<-read.table("AD.txt",fill = TRUE, header = TRUE)
head(data_ad)


save.image(file='myEnvironment.RData')

 ####################################################
# 4. merge fnbmd with AD 
####################################################

merged_datafnbmd_ad <- merge( merged_datafn2, data_ad, by = "variant_id", all.x = TRUE)
head(merged_datafnbmd_ad)


# Remove rows with NA values which means this variant_id is not present in AD_data
subset_merged_datafnbmd_ad <- merged_datafnbmd_ad[complete.cases(merged_datafnbmd_ad), ]
head(subset_merged_datafnbmd_ad)
############## save merged data 
write.table(subset_merged_datafnbmd_ad, file = "/work/larylab/NAYEMA/BIVARIATE_GWAS/Bivariate_gwas/data/merged_fnbmd_ad_common_variant_id.txt", sep = "\t", quote = FALSE, row.names = FALSE)

###############################
## Extract only columns required
subset_merged_datafnbmd_ad2 <- subset_merged_datafnbmd_ad[, c(1:16, 20:23)]
head(subset_merged_datafnbmd_ad2)

### Rename the columns for harmonization 

colnames(subset_merged_datafnbmd_ad2) <- c("SNP", "chromosome.outcome", "position.outcome", "rs_number", "effect_allele.outcome", 
                                           "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome","pval.exposure",
                                           "chromosome.exposure", "position.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure","beta.exposure", 
                                           "se.exposure","n.cases","n.control")



############################
### The alleles for AD and fnbmd are not harmonized therefore need to harmonize them 




################################
#LOAD IN AD DATA WITH FNBMD FOR COMMON SNPS

AD_data <- subset_merged_datafnbmd_ad2[, c("SNP", "chromosome.exposure", "position.exposure", "effect_allele.exposure", "other_allele.exposure","beta.exposure", 
                                           "se.exposure","eaf.exposure","pval.exposure")]
head(AD_data)
str(AD_data)
AD_data$exposure<-"AD"
AD_data$id.exposure<-"AD.1"


#LOAD IN FNBMD DATA WITH AD FOR COMMON SNPS

FNBMD_data <- subset_merged_datafnbmd_ad2[, c("SNP", "chromosome.outcome", "position.outcome", "effect_allele.outcome", "other_allele.outcome","beta.outcome", 
                                           "se.outcome","eaf.outcome","pval.outcome")]
head(FNBMD_data)
str(FNBMD_data)
FNBMD_data$outcome<-"FNBMD"
FNBMD_data$id.outcome<-"FNBMD.1"




###############################################
# 5. HARMONISE

##############################################

dat1 <- harmonise_data( 
  exposure_dat = AD_data,
  outcome_dat = FNBMD_data
)
head(dat1)

######## save the harmonised data
write.table(dat1, file = "/work/larylab/NAYEMA/BIVARIATE_GWAS/Bivariate_gwas/data/harmonised_fnbmd_ad_common_variant_id.txt", sep = "\t", quote = FALSE, row.names = FALSE)
###############




########### MAKE DATA FOR ELX

ELX_AD_FNBMD <- dat1[, c("SNP", "chromosome.exposure", "position.exposure","effect_allele.exposure","other_allele.exposure","beta.exposure","se.exposure","pval.exposure","eaf.exposure","beta.outcome","se.outcome","pval.outcome","eaf.outcome")]

######## Calculate z score

ELX_AD_FNBMD$Trait1.Z <-ELX_AD_FNBMD$beta.exposure / ELX_AD_FNBMD$se.exposure
ELX_AD_FNBMD$Trait2.Z <- ELX_AD_FNBMD$beta.outcome / ELX_AD_FNBMD$se.outcome

############# save ELX analytic full file 

write.table(ELX_AD_FNBMD, file = "/work/larylab/NAYEMA/BIVARIATE_GWAS/Bivariate_gwas/data/analyticfile_full_ELX_AD_FNBMD.txt", sep = "\t", quote = FALSE, row.names = FALSE)


####### subset data

ELX_AD_FNBMD1 <- ELX_AD_FNBMD[, c("SNP", "Trait1.Z", "Trait2.Z","pval.exposure","pval.outcome")]

# Filter the data frame based on the condition
ELX_AD_FNBMD2 <- ELX_AD_FNBMD1[ELX_AD_FNBMD1$pval.exposure < 0.05 & ELX_AD_FNBMD1$pval.outcome < 0.05, ]

# Keep only the selected columns
ELX_AD_FNBMD3 <- ELX_AD_FNBMD2[, c("SNP", "Trait1.Z", "Trait2.Z")]
colnames(ELX_AD_FNBMD3) <- c(" SNPname", "Trait1.Z", "Trait2.Z")




# Write the data frame to a file without row names
write.table(ELX_AD_FNBMD3, file = "ELX_AD_FNBMD1.txt", sep = "\t", quote = FALSE, row.names = FALSE)



##### now use"ELX_AD_FNBMD1.txt" data to run bivariate GWAS in terminal using sbatch AND THE SCRIPT IS NAMED AS 

#result

ad_fnbmd_result<- read.table("ELX_AD_FNBMD1_result", header=TRUE)
head(ad_fnbmd_result)

# Assuming ad_fnbmd_result is your data frame
ad_fnbmd_result <- ad_fnbmd_result[order(ad_fnbmd_result$dLC.Pval), ]





#### fore arm bmd -fabmd summary statistics data for european people downloaded from musculoskeleta knowledge portal 
#### link http://www.gefos.org/?q=content/data-release-2015


fa2stu<- read.table("fa2stu.MAF0_.005.pos_.out_", header=TRUE)
mached_linesfa<-read.table("matched_lines_fabmd_no_indel.txt", header=TRUE)
head(fa2stu)
head(mached_linesfa)
colnames(mached_linesfa) <- c("chromosome","position","snp","other_allele", "reference_allele",".", "..1", "..2")
mached_linesfa <- mached_linesfa[, 1:5]
##merged_datafa <- merge( mached_linesfa, fa2stu, by = c("chromosome", "position","reference_allele"))
merged_datafa <- merge( mached_linesfa, fa2stu, by = c("chromosome", "position"))
head(merged_datafa)
merged_datafa <- merged_datafa[, 1:11,14,15]



save.image(file='myEnvironment.RData')










#### lumber spine bmd-lsbmd summary statistics data for european people downloaded from musculoskeleta knowledge portal 
#### link http://www.gefos.org/?q=content/data-release-2015

ls2stu<- read.table("ls2stu.MAF0_.005.pos.out_", header=TRUE)
matched_linesls<- read.table("matched_lines_lsbmd_no_indel.txt",header=TRUE)
head(ls2stu)
head(matched_linesls)
colnames(matched_linesls) <- c("chromosome","position","snp","other_allele", "reference_allele",".", "..1", "..2")
matched_linesls <- matched_linesls[, 1:5]
merged_datals <- merge(matched_linesls, ls2stu, by = c("chromosome", "position","reference_allele"))
merged_datals <- merge(matched_linesls, ls2stu, by = c("chromosome", "position"))
head(merged_datals)

dir()

load('myEnvironment.RData')



