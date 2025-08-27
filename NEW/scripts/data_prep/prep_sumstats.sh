#!/bin/bash

#$ -l h_vmem=50G
#$ -l h_rt=00:30:00
#$ -o reports/


#$ -j y
#$ -cwd



macro=$1



source /broad/software/scripts/useuse
use R-4.1



## Prepare GWAS summary statistics 

echo preparing PGS summary statistics for $macro ...

R --no-save <<EOF
library(data.table) ; library(tidyverse)
macro=fread("../data/raw/Diet.GWAS.Merino.NatHumBehav.2022/GWAS.SumStats.BMI/meta_${macro}_CHARGE_UKBB_BMI1_noindels.txt.gz") %>%
separate(MarkerName, into=c("chr" ,"pos", "ref", "alt"), sep = ":") %>% 
mutate(var_id=paste0("chr", paste(chr, pos, ref, alt, sep=":"))) %>%
mutate(across(c(Allele1, Allele2), ~toupper(.))) %>%
mutate(ea=ifelse(Effect>0,Allele1,Allele2), nea=ifelse(Effect>0,Allele2,Allele1)) %>%
mutate(effect=abs(Effect)) %>%
select(var_id, chr, pos, ref, alt, ea, nea, effect, se=StdErr)

cat("Done preparing summary stats for ${macro}:/n")
head(macro)

macro %>% fwrite("../data/processed/pgs/merino_meta_chargeukbb_${macro}_adjBMI_forLift.txt", sep="\t")
EOF


#EOF

