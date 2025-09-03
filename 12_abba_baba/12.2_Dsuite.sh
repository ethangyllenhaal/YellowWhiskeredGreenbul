#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=ABBA-BABA_chrom
#SBATCH --nodes=1 --ntasks=2
#SBATCH --partition nocona
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-34

workdir1=/lustre/scratch/arrice/greenbul_introgression/ABBA_BABA1/
workdir2=/lustre/scratch/arrice/greenbul_introgression/ABBA_BABA2/
source ~/anaconda3/etc/profile.d/conda.sh
conda activate phylostats_env

chrom_array1=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir1}/scaffolds.txt | tail -n1 )
chrom_array2=$( head -n${SLURM_ARRAY_TASK_ID} ${workdir2}/scaffolds.txt | tail -n1 )

~/Dsuite/Build/Dsuite Dtrios -n ${chrom_array1} -t ${workdir1}/species_tree_1.nwk ${chrom_array1}__ABBA_BABA.recode.vcf.gz \
${workdir1}/Greenbul_Sets_1.txt -o ${workdir1}/Dstats_by_chrom_Dsuite/Greenbul_Sets

~/Dsuite/Build/Dsuite Dtrios -n ${chrom_array2} -t ${workdir2}/species_tree_1.nwk ${chrom_array2}__ABBA_BABA.recode.vcf.gz \
${workdir2}/Greenbul_Sets_1.txt -o ${workdir2}/Dstats_by_chrom_Dsuite/Greenbul_Sets
