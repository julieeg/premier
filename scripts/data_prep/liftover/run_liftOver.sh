#!/bin/bash

#$ -l h_vmem=50G
#$ -l h_rt=0:15:00
#$ -o reports/

#$ -j y
#$ -cwd


macro=$1

file_to_lift=../data/processed/pgs/merino_meta_chargeukbb_${macro}_adjBMI_forLift.txt
liftFrom=hg19
liftTo=hg38


source /broad/software/scripts/useuse
use R-4.1


Rscript --vanilla ../scripts/data_prep/liftOver.R $file_to_lift $liftFrom $liftTo

#EOF