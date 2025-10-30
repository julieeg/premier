# Rscript function to run liftOver 
# Written by Julie E. Gervis


## Note:
# hg19(GRCh37)->hg38(GRCh38)
# hg38(GRCh38)->hg19(GRCh37)

## Command args
args=commandArgs(trailingOnly=T)
path_to_file=args[1] 
liftFrom="hg19" #args[2]
liftTo="hg38" #args[3]

file_to_lift <- basename(path_to_file)
chain <- paste0(liftFrom, "To", gsub("h","H",liftTo))
file_lifted=gsub(".txt", paste0("_",liftTo, ".txt"), file_to_lift)


## Load required packages

library(data.table)
library(tidyverse)

library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(R.utils)



## Download chain.links to opt folder

#hg19ToHg38=https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz ; gzip -d opt/hg19ToHg38.over.chain.gz
#hg38ToHg19=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz ; gzip -d opt/hg19ToHg38.over.chain.gz


## Load file to lift
#file <- fread(path_to_file)


## Write R function to perform liftOver

liftOver <- function(path_to_file, liftFrom, liftTo, file_lifted, chain, snp_col="SNP", chr_col="CHR", pos_col="BP", ea_col="EA", nea_col="NEA", beta_col="BETA", se_col="SE") {
  
  dat.all <- fread(path_to_file) %>% select(SNP=snp_col, CHR=chr_col, POS=pos_col, EA=ea_col, NEA=nea_col, BETA=beta_col, SE=se_col)
  
  granges_input <-GRanges(
    seqnames = Rle(paste0("chr", dat.all$CHR)),  # Add 'chr' prefix
    ranges = IRanges(start = dat.all$POS, 
                     end = dat.all$POS + 1),  # Use end = start + 1 for SNPs
    rsid = dat.all$SNP
)
  
  seqlevelsStyle(granges_input) = "UCSC"
  
  chain_file <- rtracklayer::import.chain(paste0("../../opt/", chain,".over.chain"))  # Chain file for liftFrom to liftTo
  ranges.LiftTo <- rtracklayer::liftOver(granges_input, chain_file)
  
  # Extract the converted coordinates
  dat.liftTo <- as.data.frame(unlist(ranges.LiftTo)) %>%
    mutate(CHR=as.numeric(gsub("chr", "", seqnames)), POS=as.numeric(start), SNP=rsid)
  
  # Merge lifted data 
  dat.lifted <- dat.all %>% select(SNP, EA, NEA, BETA, SE) %>%
    left_join(dat.liftTo, by = "SNP", relationship = "many-to-many") %>%
    select(CHR, SNP, POS, POS, EA, NEA, BETA, SE)
  
  # Save lifted files to a sub-folder
  return(dat.lifted) 
}

# Save liftover summary stats file
liftOver(path_to_file, liftFrom, liftTo, file_lifted, chain) %>% fwrite(paste0(dirname(path_to_file), "/", file_lifted), sep="\t", col.names=T)


##EOF


