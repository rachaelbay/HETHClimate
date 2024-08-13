#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=PRUNE
#################
#a file for job output, you can check job progress for each sample
#SBATCH --output=out/prune.%j.out
#################
# a file for errors from the job for each sample
#SBATCH --error=out/prune.%j.err
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
#assigns primary group to CBossu's Xsede allocation, Rachaels is deb200006p (old=bi5613p), mine is mcb200015p (old=mcz3afp)
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


FAI="/ocean/projects/deb200006p/adamsne2/HETH/reference/GCF_009819885.2_bCatUst1.pri.v2_genomic.fna.NoU.fai"
DIR="/ocean/projects/deb200006p/rachbay/HETH/data"
############
## All 219
#############
BEAG="/ocean/projects/deb200006p/adamsne2/HETH/analyses/toSWTH/HETH.SWTH.219.noSex.noUn.beagle.gz"

 Turn beagle into plink
/jet/home/adamsne2/programs/angsd/angsd \
-beagle $BEAG \
-P 10 \
-fai $FAI -doMaf 4 -doPlink 2 -doGeno 1 \
-out $DIR/HETH.219


# Filter by ld
plink --tped $DIR/HETH.219.tped \
--tfam $DIR/HETH.219.tfam \
--indep-pairwise 25 10 0.7 \
--allow-extra-chr \
--out $DIR/HETH.219.R2.pruned

# Filter beagle file
zcat $DIR/HETH.219.beagle.gz | head -1 > $DIR/HETH.219.pruned.beagle
zcat $DIR/HETH.219.beagle.gz | \
awk 'NR==FNR{_[$1];next}$1 in _' $DIR/HETH.219.R2.pruned.prune.out - \
>> $DIR/HETH.219.pruned.beagle
gzip $DIR/HETH.219.pruned.beagle

###############
## Contemporary only
##############
BEAG="/ocean/projects/deb200006p/rachbay/HETH/data/HETH_contemp.2005.noUn.beagle.gz"

# Turn beagle into plink
/jet/home/adamsne2/programs/angsd/angsd \
-beagle $BEAG \
-P 10 \
-fai $FAI -doMaf 4 -doPlink 2 -doGeno 1 \
-out $DIR/HETH.contemp


# Filter by ld
plink --tped $DIR/HETH.contemp.tped \
--tfam $DIR/HETH.contemp.tfam \
--indep-pairwise 25 10 0.7 \
--allow-extra-chr \
--out $DIR/HETH.contemp.R2.pruned

# Filter beagle file
zcat $DIR/HETH_contemp.2005.noUn.beagle.gz | head -1 > $DIR/HETH.contemp.pruned.beagle
zcat $DIR/HETH_contemp.2005.noUn.beagle.gz | \
awk 'NR==FNR{_[$1];next}$1 in _' $DIR/HETH.contemp.R2.pruned.prune.out - \
>> $DIR/HETH.contemp.pruned.beagle
gzip $DIR/HETH.contemp.pruned.beagle

##Make separate files for males and females (for both pruned and full)
PROG="/ocean/projects/deb200006p/rachbay/programs"

zcat $DIR/HETH.contemp.pruned.beagle.gz | java -jar $PROG/filtercolumns.jar -1 $DIR/FemaleList.contemp.txt | gzip > $DIR/HETH.contemp.pruned.male.beagle.gz
zcat $DIR/HETH.contemp.pruned.beagle.gz | java -jar $PROG/filtercolumns.jar -1 $DIR/MaleList.contemp.txt | gzip > $DIR/HETH.contemp.pruned.female.beagle.gz
zcat $DIR/HETH_contemp.2005.noUn.beagle.gz | java -jar $PROG/filtercolumns.jar -1 $DIR/FemaleList.contemp.txt | gzip > $DIR/HETH.contemp.male.beagle.gz
zcat $DIR/HETH_contemp.2005.noUn.beagle.gz | java -jar $PROG/filtercolumns.jar -1 $DIR/MaleList.contemp.txt | gzip > $DIR/HETH.contemp.female.beagle.gz

##Make separate files for sex for full dataset
zcat $DIR/HETH.219.pruned.beagle.gz | java -jar $PROG/filtercolumns.jar -1 $DIR/FemaleList.all.txt | gzip > $DIR/HETH.all.pruned.male.beagle.gz
zcat $DIR/HETH.219.beagle.gz | java -jar $PROG/filtercolumns.jar -1 $DIR/FemaleList.all.txt | gzip > $DIR/HETH.all.male.beagle.gz
zcat $DIR/HETH.219.pruned.beagle.gz | java -jar $PROG/filtercolumns.jar -1 $DIR/MaleList.all.txt | gzip > $DIR/HETH.all.pruned.female.beagle.gz
zcat $DIR/HETH.219.beagle.gz | java -jar $PROG/filtercolumns.jar -1 $DIR/MaleList.all.txt | gzip > $DIR/HETH.all.female.beagle.gz
