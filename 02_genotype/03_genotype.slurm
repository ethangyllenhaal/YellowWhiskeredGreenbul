#!/bin/bash

#SBATCH --chdir=./
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --job-name=filter_greenbul
#SBATCH --output=_geno_out
#SBATCH --error=_geno_error
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=egyllenh@ttu.edu
#SBATCH --mem-per-cpu=4G
#SBATCH --time=48:00:00
#SBATCH --array=1-51

# bcftools is the important part here
source activate bcftools-env

threads=8
reads=/lustre/scratch/egyllenh/redo_greenbul/clean_reads/

# define main working directory
workdir=/lustre/scratch/egyllenh/redo_greenbul

# base name of fastq files, intermediate files, and output files
basename_array=$( head -n${SLURM_ARRAY_TASK_ID} sample_list_filter.txt | tail -n1 )
# define the reference genome
ref=/lustre/scratch/egyllenh/redo_greenbul/eula_ragtag_correct_chr0/ragtag.scaffold.fasta

# run bcftools to genotype
bcftools mpileup --skip-indels -C 0 -d 200 --min-MQ 10 --threads ${threads} \
	 -f ${ref} ${workdir}/01_bam_files/${basename_array}_final.bam | \
    bcftools call -m --threads ${threads} -o ${workdir}/02_vcf/${basename_array}.vcf

# bgzip
bgzip ${workdir}/02_vcf/${basename_array}.vcf

#tabix
tabix ${workdir}/02_vcf/${basename_array}.vcf.gz

# filter individual vcf files
bcftools view -i 'MIN(DP)>5' ${workdir}/02_vcf/${basename_array}.vcf.gz > \
	 ${workdir}/03_vcf/${basename_array}.vcf

# bgzip
bgzip ${workdir}/03_vcf/${basename_array}.vcf

#tabix
tabix ${workdir}/03_vcf/${basename_array}.vcf.gz
