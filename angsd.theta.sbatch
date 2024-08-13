#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=theta
#################
#a file for job output, you can check job progress
#SBATCH --output=out/theta.%j.out
#################
# a file for errors from the job
#SBATCH --error=out/theta.%j.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH -t 48:00:00
#################
#resource
#SBATCH -p RM
#################
##SBATCH -N 1
#SBATCH --ntasks 144
#################
#SBATCH --mem=240G
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
NGSTOOLS="/jet/home/adamsne2/programs/ngsTools"

REFERENCE="/ocean/projects/deb200006p/adamsne2/HETH/reference/GCF_009819885.2_bCatUst1.pri.v2_genomic.fna"
BAMDIR="/ocean/projects/deb200006p/adamsne2/HETH/bwa_map/toSWTH/"
OUTDIR="/ocean/projects/deb200006p/rachbay/HETH/thetas/"

YEAR=$1
BAMLIST=${OUTDIR}/bamlist.$YEAR.txt

grep -f ${OUTDIR}/samples.$YEAR.txt HETH.219.bamlist > $BAMLIST

###Generate SFS
$ANGSD -bam $BAMLIST -P 24 -ref $REFERENCE -anc $REFERENCE \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
	-trim 0 -skipTriallelic 1 -baq 1 -C 50 \
	-minMapQ 20 -minQ 20  -minInd 5 -doCounts 1 -doMajorMinor 5 \
	-GL 1 -doSaf 1 -doMaf 1 -setMinDepth 1 \
	-out ${OUTDIR}/$YEAR

$REALSFS ${OUTDIR}/$YEAR.saf.idx -P 6 -fold 1 > ${OUTDIR}/$YEAR.sfs

$NGSTOOLS/angsd/angsd -bam $BAMLIST -P 48 -ref $REFERENCE -anc $REFERENCE \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
        -trim 0 -skipTriallelic 1 -baq 1 -C 50 \
        -minMapQ 20 -minQ 20  -minInd 5 -doMajorMinor 5 \
        -GL 1 -doSaf 1 -doCounts 1 -setMinDepth 3 \
        -pest ${OUTDIR}/$YEAR.sfs -doThetas 1 -out ${OUTDIR}/$YEAR.thet


$NGSTOOLS/angsd/misc/thetaStat do_stat -out ${OUTDIR}/$YEAR.thet.thetas.idx

$NGSTOOLS/angsd/misc/thetaStat do_stat ${OUTDIR}/$YEAR.thet.thetas.idx -win 10000 -step 10000 -outnames ${OUTDIR}/$YEAR.thetas
