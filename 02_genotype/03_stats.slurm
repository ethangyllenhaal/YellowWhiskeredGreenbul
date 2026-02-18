#!/bin/bash

#SBATCH --chdir=./
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --job-name=stat_greenbul
#SBATCH --output=_stat_out
#SBATCH --error=_stat_error
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=egyllenh@ttu.edu
#SBATCH --mem-per-cpu=4G
#SBATCH --time=48:00:00
#SBATCH --array=1-52

# samtools is the relevant bit here
source activate bcftools-env

# define main working directory
workdir=/lustre/scratch/egyllenh/redo_greenbul

# base name of fastq files, intermediate files, and output files
basename_array=$( head -n${SLURM_ARRAY_TASK_ID} sample_list_full.txt | tail -n1 )

# alignment stats
echo ${basename_array} > ${workdir}/stats/${basename_array}.stats

# samtools depth sum of aligned sites
echo "samtools depth sum of aligned sites" >> ${workdir}/stats/${basename_array}.stats
samtools depth ${workdir}/01_bam_files/${basename_array}_final.bam  |  awk '{sum+=$3} END { print "Sum = ",sum}' >> ${workdir}/stats/${basename_array}.stats

echo "samtools mean depth" >> ${workdir}/stats/${basename_array}.stats
samtools depth ${workdir}/01_bam_files/${basename_array}_final.bam |  awk '{sum+=$3} END { print "Average = ",sum/NR}' >> ${workdir}/stats/${basename_array}.stats

# proportion dupes
echo "proportion duplicates" >> ${workdir}/stats/${basename_array}.stats
head -n8 ${workdir}/01_bam_files/${basename_array}_markdups_metric_file.txt | tail -n1 | cut -f9 >> ${workdir}/stats/${basename_array}.stats

# number of genotyped sites passing minimum depth filter
echo "sites genotyped" >> ${workdir}/stats/${basename_array}.stats
gzip -cd ${workdir}/03_vcf/${basename_array}.vcf.gz | grep -v "^#" | wc -l >> ${workdir}/stats/${basename_array}.stats

# number of reads in read 1 then read 2 files
echo "number of reads (read 1 then read 2)" >> ${workdir}/stats/${basename_array}.stats
echo $(zcat clean_reads/${basename_array}_R1.fastq.gz|wc -l)/4|bc >> ${workdir}/stats/${basename_array}.stats
echo $(zcat clean_reads/${basename_array}_R2.fastq.gz|wc -l)/4|bc >> ${workdir}/stats/${basename_array}.stats
