#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=ABBA-BABA_chrom
#SBATCH --nodes=1 --ntasks=2
#SBATCH --partition nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-34
#SBATCH --mail-user=arrice@ttu.edu
#SBATCH --mail-type=ALL

workdir=/lustre/scratch/arrice/greenbul_introgression/ABBA_BABA2/Dstats_by_chrom_Dsuite
source ~/anaconda3/etc/profile.d/conda.sh
conda activate phylostats_env

chrom_array=$( head -n${SLURM_ARRAY_TASK_ID} scaffolds.txt | tail -n1 )

~/Dsuite/Build/Dsuite Dtrios -n ${chrom_array} -t species_tree_1.nwk ${chrom_array}__ABBA_BABA.recode.vcf.gz \
Greenbul_Sets)_1.txt -o ${workdir}/Greenbul_Sets
