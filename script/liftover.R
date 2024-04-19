library(rtracklayer)
library(optparse)
library(bigsnpr)
require(dplyr)
require(data.table)
require(dplyr)

library(data.table)

# File path
file_path <- "/work/larylab/NAYEMA/BIVARIATE_GWAS/MA_MA_meta/Meta.Analysis/GCST90027158_buildGRCh38.tsv.gz"

# Read the compressed TSV file
data <- fread(file_path)

# Display the first few rows of the data
head(data)
setnames(data, old = c("chromosome", "base_pair_location"), new = c("chr", "pos"))






data <- fread("/work/larylab/NAYEMA/BIVARIATE_GWAS/fn2stu.MAF0_.005.pos.out_",header=T)
head(data)
setnames(data, old = c("chromosome", "position"), new = c("chr", "pos"))
data <- as.data.frame(data)

# Print the data without row indices
head(data)



################################################################################

snp_modifyBuild <- function(info_snp, liftOver,
                            from = "hg38", to = "hg19",
                            check_reverse = TRUE) {
  
  if (!all(c("chr", "pos") %in% names(info_snp)))
    stop("Expecting variables 'chr' and 'pos' in input 'info_snp'.")
  
  # Need BED UCSC file for liftOver
  info_BED <- with(info_snp, data.frame(
    chrom = paste0("chr", sub("^0", "", chr)),
    start = pos - 1L, end = pos,
    id = seq_along(pos)))
  
  BED <- tempfile(fileext = ".BED")
  bigreadr::fwrite2(stats::na.omit(info_BED),
                    BED, col.names = FALSE, sep = " ", scipen = 50)
  
  # Need chain file
  liftOver_path <- "/work/larylab/software/liftOver"
  
  # Make the liftOver executable
  #system(paste("chmod +x", shQuote(liftOver_path)))
  
  # Construct the URL using 'paste0'
  url <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
  chain <- tempfile(fileext = ".over.chain.gz")
 utils::download.file(url, destfile = chain, quiet = TRUE)
#chain <- "/work/larylab/NAYEMA/liftover/hg19ToHg38.over.chain.gz"
  # Run liftOver (usage: liftOver oldFile map.chain newFile unMapped)
  lifted <- tempfile(fileext = ".BED")
  system2(liftOver_path, c(BED, chain, lifted, tempfile(fileext = ".txt")))
  
  # Read the mapped positions and perform QC
  new_pos <- bigreadr::fread2(lifted, nThread = 1)
  is_bad <- vctrs::vec_duplicate_detect(new_pos$V4) |
    (new_pos$V1 != info_BED$chrom[new_pos$V4])
  new_pos <- new_pos[which(!is_bad), ]
  
  # Update SNP positions
  pos0 <- info_snp$pos
  info_snp$pos <- NA_integer_
  info_snp$pos[new_pos$V4] <- new_pos$V3
  
  # Reverse liftover if specified
  if (check_reverse) {
    pos2 <- suppressMessages(
      Recall(info_snp, liftOver_path, from = to, to = from, check_reverse = FALSE)$pos)
    info_snp$pos[pos2 != pos0] <- NA_integer_
  }
  
  # Print message indicating variants not mapped
  message("%d variants have not been mapped.", sum(is.na(info_snp$pos)))
  
  info_snp
}
################################################################################


################################################################################
lifted <- snp_modifyBuild(data, "/work/larylab/software/liftOver", from = "hg38",
                          to = "hg19",
                          check_reverse = TRUE)

fwrite(lifted, file = "/work/larylab/NAYEMA/BIVARIATE_GWAS/Bivariate_gwas/data/liffted", col.names = TRUE, sep = "\t")
head(lifted)
head(lifted, 100)




##################

################################################################################

#' Modify genome build
#'
#' Modify the physical position information of a data frame
#' when converting genome build using executable *liftOver*.
#'
#' @param info_snp A data frame with columns "chr" and "pos".
#' @param liftOver Path to liftOver executable. Binaries can be downloaded at
#'   \url{https://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/liftOver} for Mac
#'   and at \url{https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver}
#'   for Linux.
#' @param from Genome build to convert from. Default is `hg18`.
#' @param to Genome build to convert to. Default is `hg19`.
#' @param check_reverse Whether to discard positions for which we cannot go back
#'   to initial values by doing 'from -> to -> from'. Default is `TRUE`.
#'
#' @references
#' Hinrichs, Angela S., et al. "The UCSC genome browser database: update 2006."
#' Nucleic acids research 34.suppl_1 (2006): D590-D598.
#'
#' @return Input data frame `info_snp` with column "pos" in the new build.
#' @export
#'
snp_modifyBuild <- function(info_snp, liftOver,
                            from = "hg19", to = "hg38",
                            check_reverse = TRUE) {
  
  if (!all(c("chr", "pos") %in% names(info_snp)))
    stop2("Expecting variables 'chr' and 'pos' in input 'info_snp'.")
  
  # Make sure liftOver is executable
  liftOver <- make_executable(normalizePath(liftOver))
  
  # Need BED UCSC file for liftOver
  info_BED <- with(info_snp, data.frame(
    # sub("^0", "", c("01", 1, 22, "X")) -> "1"  "1"  "22" "X"
    chrom = paste0("chr", sub("^0", "", chr)),
    start = pos - 1L, end = pos,
    id = seq_along(pos)))
  
  BED <- tempfile(fileext = ".BED")
  bigreadr::fwrite2(stats::na.omit(info_BED),
                    BED, col.names = FALSE, sep = " ", scipen = 50)
  
  # Need chain file
  url <- paste0("ftp://hgdownload.cse.ucsc.edu/goldenPath/", from, "/liftOver/",
                from, "To", tools::toTitleCase(to), ".over.chain.gz")
  chain <- tempfile(fileext = ".over.chain.gz")
  utils::download.file(url, destfile = chain, quiet = TRUE)
  
  # Run liftOver (usage: liftOver oldFile map.chain newFile unMapped)
  lifted <- tempfile(fileext = ".BED")
  system2(liftOver, c(BED, chain, lifted, tempfile(fileext = ".txt")))
  
  # Read the ones lifter + some QC
  new_pos <- bigreadr::fread2(lifted, nThread = 1)
  is_bad <- vctrs::vec_duplicate_detect(new_pos$V4) |
    (new_pos$V1 != info_BED$chrom[new_pos$V4])
  new_pos <- new_pos[which(!is_bad), ]
  
  pos0 <- info_snp$pos
  info_snp$pos <- NA_integer_
  info_snp$pos[new_pos$V4] <- new_pos$V3
  
  if (check_reverse) {
    pos2 <- suppressMessages(
      Recall(info_snp, liftOver, from = to, to = from, check_reverse = FALSE)$pos)
    info_snp$pos[pos2 != pos0] <- NA_integer_
  }
  
  message2("%d variants have not been mapped.", sum(is.na(info_snp$pos)))
  
  info_snp
}
