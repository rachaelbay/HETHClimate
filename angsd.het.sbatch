#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=het
#################
#a file for job output, you can check job progress
#SBATCH --output=out/het.%j.out
#################
# a file for errors from the job
#SBATCH --error=out/het.%j.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH -t 48:00:00
#################
#resource
#SBATCH -p RM-shared
#################
##SBATCH -N 1
#SBATCH --ntasks 23
#################
#SBATCH --mem=44G
#################
#assigns primary group
#SBATCH -A deb200006p
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=END
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=rachaelbay@gmail.com
#################
#now run normal batch commands
##################
#export MALLOC_PER_THREAD=1
#echo commands to stdout

set -x

ANGSD=/jet/home/adamsne2/programs/angsd/angsd
REALSFS=/jet/home/adamsne2/programs/angsd/misc/realSFS

REFERENCE="/ocean/projects/deb200006p/adamsne2/HETH/reference/GCF_009819885.2_bCatUst1.pri.v2_genomic.fna"
BAMDIR="/ocean/projects/deb200006p/adamsne2/HETH/bwa_map/toSWTH/"
OUTDIR="/ocean/projects/deb200006p/rachbay/HETH/het/"


ls $BAMDIR/*sw.mrkdup.bam > bamfiles.txt

mapfile -t FILENAMES < HETH.219.bamlist
b=${FILENAMES[$SLURM_ARRAY_TASK_ID]}

sample="$(basename $b .sw.mrkdup.bam)"

$ANGSD -i $b -anc $REFERENCE -ref $REFERENCE -P 10 \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
	-trim 0 -skipTriallelic 1 -baq 1 -C 50 \
	-minMapQ 20 -minQ 20 -doMajorMinor 4 -rf ../data/chr.list \
	-dosaf 1 -GL 1 -out "$OUTDIR"/"$sample"

$REALSFS $OUTDIR/$sample.saf.idx -fold 1 -P 10 > $OUTDIR/"$sample".ml
