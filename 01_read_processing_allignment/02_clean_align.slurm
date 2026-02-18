#!/bin/bash

#SBATCH --chdir=./
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --partition=nocona
#SBATCH --time=48:00:00
#SBATCH --job-name=filter_greenbul
#SBATCH --output=_align_out
#SBATCH --error=_align_error
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=egyllenh@ttu.edu
#SBATCH --mem-per-cpu=4G
#SBATCH --time=48:00:00
#SBATCH --array=1-52

source activate bcftools-env

threads=8
reads=/lustre/scratch/egyllenh/redo_greenbul/clean_reads/

# define main working directory
workdir=/lustre/scratch/egyllenh/redo_greenbul

# base name of fastq files, intermediate files, and output files
basename_array=$( head -n${SLURM_ARRAY_TASK_ID} sample_list_full.txt | tail -n1 )

# define the reference genome
ref=/lustre/scratch/egyllenh/redo_greenbul/eula_ragtag_correct_chr0/ragtag.scaffold.fasta

# clean with bbduk
bbduk.sh \
    in1=${workdir}/raw_reads/${basename_array}__R1.fastq.gz \
    in2=${workdir}/raw_reads/${basename_array}__R2.fastq.gz \
    out1=${workdir}/clean_reads/${basename_array}_R1.fastq.gz \
    out2=${workdir}/clean_reads/${basename_array}_R2.fastq.gz \
    minlen=50 ftl=10 qtrim=rl trimq=10 ktrim=r k=25 mink=7 \
    ref=/lustre/work/jmanthey/bbmap/resources/adapters.fa hdist=1 tbo tpe

# run bwa mem, adding read groups
bwa mem -t ${threads} ${ref} \
    -R "@RG\tID:{}\tPL:ILLUMINA\tLB:${basename_array}\tSM:${basename_array}" \
    ${reads}/${basename_array}_R1.fastq.gz \
    ${reads}/${basename_array}_R2.fastq.gz | \
    samtools view -b -f 2 -@ ${threads} - > \
    ${workdir}/01_bam_files/${basename_array}.bam

# clean up the bam file
picard CleanSam \
       I=${workdir}/01_bam_files/${basename_array}.bam \
       O=${workdir}/01_bam_files/${basename_array}_cleaned.bam

# remove the raw bam
rm ${workdir}/01_bam_files/${basename_array}.bam

# sort the cleaned bam file
picard SortSam \
       I=${workdir}/01_bam_files/${basename_array}_cleaned.bam \
       O=${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam SORT_ORDER=coordinate

# remove the cleaned bam file
rm ${workdir}/01_bam_files/${basename_array}_cleaned.bam

# remove duplicates to sorted, cleaned, and read grouped bam file (creates final bam file)
picard MarkDuplicates REMOVE_DUPLICATES=true \
     MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100 \
     M=${workdir}/01_bam_files/${basename_array}_markdups_metric_file.txt \
     I=${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam \
     O=${workdir}/01_bam_files/${basename_array}_final.bam

# remove sorted, cleaned, and read grouped bam file
rm ${workdir}/01_bam_files/${basename_array}_cleaned_sorted.bam

# index the final bam file
samtools index ${workdir}/01_bam_files/${basename_array}_final.bam
