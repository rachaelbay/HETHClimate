#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=ANNOT
#################
#a file for job output, you can check job progress for each sample
#SBATCH --output=out/annot.%j.out
#################
# a file for errors from the job for each sample
#SBATCH --error=out/annot.%j.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH -t 72:00:00
#################
#quality of service; think of it as job priority
#SBATCH -p RM-shared
#################
#number of nodes- for RM-shared -N 1
##SBATCH -N 1
#SBATCH --ntasks-per-node 6
#################
#SBATCH --mem=10G
#################
#assigns primary group
#SBATCH -A deb200006p
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=END
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=rachaelbay@gmail.com
##################
##################
##################
set -x

#
module load anaconda3
conda activate plink
module load vcftools

FAI="/ocean/projects/deb200006p/adamsne2/HETH/reference/GCF_009819885.2_bCatUst1.pri.v2_genomic.fna.NoU.fai"
DIR="/ocean/projects/deb200006p/rachbay/HETH/annot"
VCF="/ocean/projects/deb200006p/rachbay/HETH/data/HETH.219.nosex.vcf"

###Note: vcf and gff need to be edited so that chr and sample names do not include special characters
cat ../data/HETH.219.nosex.vcf | \
sed 's/.sw.mrkdup.bam//g' | \
sed 's/\/ocean\/projects\/deb200006p\/adamsne2\/HETH\/bwa_map\/toSWTH\///g' | \
sed -Ee 's/NC\_(\S*)\.1/chr\1/g' | \
sed -Ee 's/NW\_(\S*)\.1/scaff\1/g' | \
sed -Ee 's/NC\_(\S*)\.2/chr\1/g' > $DIR/HETH.219.fixednames.vcf

zcat $DIR/GCF_009819885.2_bCatUst1.pri.v2_genomic.gff.gz | \
sed -Ee 's/NC\_(\S*)\.1/chr\1/g' | \
sed -Ee 's/NW\_(\S*)\.1/scaff\1/g' | \
sed -Ee 's/NC\_(\S*)\.2/chr\1/g' > $DIR/fixedchrom.gff

vcftools --vcf $DIR/HETH.219.fixednames.vcf --012 --out $DIR/HETH.219.fixednames

cd $DIR

python3 LD-annot/LD-annot0.4.py HETH.219.fixednames.vcf fixedchrom.gff HETH.219.fixednames.012.pos gene 0.9 HETH_0.9
