# LiftOver to Change Genome Build

## UCSC Tool with Associated Chain File

As I am working on Linux and want to change build from hg19 to hg38:

1. [Liftover Tool](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver)
2. [Chain File](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz)

### Inside `larylab/software` Directory 

LiftOver and useful chain files are already installed, so you can skip these downloading steps.

### Download GWAS Summary Statistics 

- Check the header.
  - We need only chromosome and basepair position (start-end) information for LiftOver.
  - But want all the other columns intact.

### Make the Input File for LiftOver Using AWK and Carry Over Rest of the Columns 

Navigate to your own work path.

#### For Zipped File

```bash
zcat /work/larylab/NAYEMA/BIVARIATE_GWAS/MA_MA_meta/Meta.Analysis/GCST90027158_buildGRCh38.tsv.gz | awk 'BEGIN {OFS="\t"} {print $3, $4-1, $4, $0}' > ucsc.input.bed
```

### for unzipped file 
```bash
awk 'BEGIN {OFS="\t"} {print $3, $4-1, $4, $0}' GWAS_summary_stat_file > ucsc.input.bed
```
### Add chr infront of chromosome numbers and delete the header 
 
 ```bash
 awk 'BEGIN {OFS="\t"} NR>1 {$1 = "chr" $1; print}' ucsc.input.bed > new_file.bed
 ```
 
 ## Your file is now ready to be executable by liftOver:
 
 ```bash
 ./liftOver -bedPlus=3 -tab new_file.bed hg38ToHg19.over.chain ucsc.output.bed ucsc.unmapped.bed
 ```