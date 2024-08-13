#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=CAND
#################
#a file for job output, you can check job progress for each sample
#SBATCH --output=out/cand.%j.out
#################
# a file for errors from the job for each sample
#SBATCH --error=out/cand.%j.err
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
##SBATCH --ntasks-per-node 2
#################
##SBATCH --mem=4G
#################
#assigns primary group
#SBATCH -A deb200006p
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=END
#################
##################
set -x

#
BEAG="/ocean/projects/deb200006p/rachbay/HETH/data/HETH.219.beagle.gz"
FAI="/ocean/projects/deb200006p/adamsne2/HETH/reference/GCF_009819885.2_bCatUst1.pri.v2_genomic.fna.NoU.fai"
OUTDIR="/ocean/projects/deb200006p/rachbay/HETH/gwas"
CANDIR="/ocean/projects/deb200006p/rachbay/HETH/cands/cutoff5"

PRE=$1
#PRE=HETH.contemp.male.PCcovar.lrt1

##All sig snps
zcat $OUTDIR/$PRE.gz | awk '($8>0 && $8<0.00001)' | awk '{print $1"_"$2}' > $CANDIR/$PRE.sites

##Filter beagle for significant snps only
zcat $BEAG | head -1 > $CANDIR/$PRE.beagle
zcat $BEAG | \
awk 'NR==FNR{_[$1];next}$1 in _' $CANDIR/$PRE.sites - \
>> $CANDIR/$PRE.beagle
gzip $CANDIR/$PRE.beagle

##Run beagle to get dosage and imputed genotypes (note - this has to be the old beagle 3)
BEAGLE="/jet/home/adamsne2/programs/beagle.jar"

java -Xmx4000m -jar $BEAGLE \
like=$CANDIR/$PRE.beagle.gz \
out="$CANDIR"/impute
