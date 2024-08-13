#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=PCA
#################
#a file for job output, you can check job progress for each sample
#SBATCH --output=out/pca.%j.out
#################
# a file for errors from the job for each sample
#SBATCH --error=out/pca.%j.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH -t 12:00:00
#################
#quality of service; think of it as job priority
#SBATCH -p RM-shared
#################
#number of nodes- for RM-shared -N 1
##SBATCH -N 1
#SBATCH --ntasks-per-node 12
#################
#SBATCH --mem=23G
#################
#assigns primary group
#SBATCH -A deb200006p
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=END
#################

##################
set -x


module load anaconda3
conda activate pcangsd

#PCA on entire pruned dataset
BEAG="/ocean/projects/deb200006p/rachbay/HETH/data/HETH.219.pruned.beagle.gz"
OUTDIR="/ocean/projects/deb200006p/rachbay/HETH/pca"
pcangsd -b $BEAG -o $OUTDIR/HETH.219.pruned --threads 6

#pruned males and females separately (contemporary)
BEAGM="/ocean/projects/deb200006p/rachbay/HETH/data/HETH.contemp.pruned.male.beagle.gz"
pcangsd -b $BEAGM -o $OUTDIR/HETH.contemp.pruned.male --threads 6

BEAGF="/ocean/projects/deb200006p/rachbay/HETH/data/HETH.contemp.pruned.female.beagle.gz"
pcangsd -b $BEAGF -o $OUTDIR/HETH.contemp.pruned.female --threads 6

BEAGC="/ocean/projects/deb200006p/rachbay/HETH/data/HETH.contemp.pruned.beagle.gz"
pcangsd -b $BEAGC -o $OUTDIR/HETH.contemp.pruned --threads 6

#pruned males and females separately (all)
BEAGAM="/ocean/projects/deb200006p/rachbay/HETH/data/HETH.all.pruned.male.beagle.gz"
pcangsd -b $BEAGAM -o $OUTDIR/HETH.all.pruned.male --threads 6

BEAGAF="/ocean/projects/deb200006p/rachbay/HETH/data/HETH.all.pruned.female.beagle.gz"
pcangsd -b $BEAGAF -o $OUTDIR/HETH.all.pruned.female --threads 6
