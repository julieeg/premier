#!/bin/bash

#$ -l h_vmem=10G
#$ -l h_rt=1:00:00
#$ -o reports/

#$ -j y
#$ -cwd



macro=$1




mgbb_bgen_dir=/humgen/florezlab/users/rmandla/MGBB-genotypes/mega-gsa_a-d/bgen_v2
pgs_dir=../data/processed/pgs
opt=../../opt



source /broad/software/scripts/useuse

reuse -q Anaconda3
source activate $opt/bgen





## Extract bgen with all variants to calculate pgs

#awk -v CHR=$chr '{ if (NR>1 && $2 == CHR) print $1 }' ../data/processed/merino_sumstats_carb.txt > ${pgs_dir}/${macro}_pgs_chr${chr}.snplist

for i in {1..22} ; do
bgenix \
-g $mgbb_bgen_dir/MGBB.TOPMED.65K.20221012.RM.chr${i}.bgen \
-incl-rsids ../data/raw/merino_allsnps_snplist_hg38.txt \
> ../data/temp/merino_allsnps_snplist_chr${i}.bgen
done

cat-bgen -g ../data/temp/merino_allsnps_snplist_chr*.bgen -og ../data/temp/merino_allsnps_snplist.bgen -clobber 



## Calculate PGS using plink

$opt/plink2 \
--bgen ../data/temp/merino_allsnps_snplist.bgen ref-first \
--sample $mgbb_bgen_dir/MGBB.samples \
--q-score-range ../data/raw/q.ranges.txt merino_allsnps_macro_snpvals \
--score ../data/raw/merino_allsnps_macro_sumstats_forPGS.txt 1 3 4 header center \
--memory 5000 \
--out ../data/processed/pgs/macro_allsnps_pgs


#--keep ../data/raw/MGBB.premier.samples \
#--read-freq ${pgs_dir}/carb_pgs_chr22.afreq \


#EOF

