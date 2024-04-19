# Install and load the CMplot package
install.packages("CMplot")
library(CMplot)




# my data

setwd("/work/larylab/NAYEMA/BIVARIATE_GWAS/Bivariate_gwas/")
df<-read.table("/work/larylab/NAYEMA/BIVARIATE_GWAS/Bivariate_gwas/result/ELX_AD_FNBMD1_result",header=TRUE)

df2<-read.table("/work/larylab/NAYEMA/BIVARIATE_GWAS/Bivariate_gwas/data/analyticfile_full_ELX_AD_FNBMD.txt",header=TRUE)
head(df)
head(df2)

merged_data <- merge(df, df2, by = "SNP")
head(merged_data)


#ad_gwas
# Subset the merged_data dataframe to include only the specified columns
all_gwas <- merged_data[, c("SNP","chromosome.outcome", "position.outcome", "pval.exposure", "pval.outcome", "dLC.Pval")]

# Rename the columns for clarity
colnames(all_gwas) <- c("SNP", "Chromosome", "Position", "AD","FNBMD","BIVARIATE")

snp_row <- all_gwas[all_gwas$SNP == "rs7175038", ]

all_gwas <- all_gwas[order(all_gwas$Chromosome, all_gwas$Position), ]


CMplot(all_gwas,type="p",plot.type="c",chr.labels=paste("Chr",c(1:22),sep=""),r=0.4,cir.axis=TRUE,
         outward=FALSE,cir.axis.col="black",cir.chr.h=1.3,chr.den.col="black",file="jpg",
         file.name="MIAMI9",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10,cir.axis.grid=TRUE,signal.line=1,verbose=TRUE)
# to remove the grid line in circles, add parameter cir.axis.grid=FALSE
# file.name: specify the output file name, the default is corresponding column name

CMplot(all_gwas,type="p",plot.type="c",r=0.4,chr.labels=paste("Chr",c(1:22),sep=""),
       threshold=c(5e-8,1e-5),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red",
                                                                                              "blue"),signal.line=1,signal.col=c("red","green"),outward=FALSE,file="jpg",file.name="miami",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10,cir.axis.grid=FALSE)
