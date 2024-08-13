#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=GWAS
#################
#a file for job output, you can check job progress for each sample
#SBATCH --output=out/gwas.%j.out
#################
# a file for errors from the job for each sample
#SBATCH --error=out/gwas.%j.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH -t 6:00:00
#################
#quality of service; think of it as job priority
#SBATCH -p RM-shared
#################
#number of nodes- for RM-shared -N 1
##SBATCH -N 1
##SBATCH --ntasks-per-node 24
#################
##SBATCH --mem=128G
#################
#assigns primary group
#SBATCH -A deb200006p
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=END
##################
set -x

#
BEAGF="/ocean/projects/deb200006p/rachbay/HETH/data/HETH.contemp.male.beagle.gz"
FAI="/ocean/projects/deb200006p/adamsne2/HETH/reference/GCF_009819885.2_bCatUst1.pri.v2_genomic.fna.NoU.fai"
OUTDIR="/ocean/projects/deb200006p/rachbay/HETH/gwas"

# PCs and library covariate - this is the one to use!
/jet/home/adamsne2/programs/angsd/angsd \
-beagle $BEAGF \
-P 10 \
-fai $FAI -doAsso 2 -doMaf 4 -Pvalue 1 -minMaf 0.1 \
-sampleFile $OUTDIR/HETH.contemp.male.sample \
-whichPhe Tarsus \
-whichCov Lib,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
-out $OUTDIR/HETH.contemp.male.PCcovar


#Tarsus as covariate for Bill and Wing
/jet/home/adamsne2/programs/angsd/angsd \
-beagle $BEAGF \
-P 10 \
-fai $FAI -doAsso 2 -doMaf 4 -Pvalue 1 -minMaf 0.1\
-sampleFile $OUTDIR/HETH.contemp.male.tarsuscov.sample \
-whichPhe Bill,Wing \
-whichCov Lib,Tarsus,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
-out $OUTDIR/HETH.contemp.male.tarsPCcovar
