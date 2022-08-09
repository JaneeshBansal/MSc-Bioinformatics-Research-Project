#!/bin/sh
#$ -cwd              # Set the working directory for the job to the current directory
#$ -j y              # Combine the standard output and the standard error stream
#$ -pe smp 16         # Request 16 cores
#$ -l h_rt=72:0:0    # Request 72 hours runtime
#$ -l h_vmem=6G      # Request 6GB RAM
#$ -m bea            # Send an email once job starts, finishes and aborts

module load bwa
module load anaconda3


conda activate pairtools

bwa mem -5SP -T0 -t ${NSLOTS} input/bwa/hg38.fa /data/home/bt211032/scratch/research_project/2022-04-08-mappy_seq_monk/input/SRR15503419_1.fastq /data/home/bt211032/scratch/research__project/2022-04-08-mappy_seq_monk/input/SRR15503419_2.fastq | \
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in 2 --nproc-out 2 --chroms-path input/hg38.genome | \
pairtools sort --nproc 4 --tmpdir=/data/home/bt211032/scratch/research_project/2022-04-13-bwa_pairtools/tmp | \
pairtools dedup --nproc-in 2 --nproc-out 2 --mark-dups | \
pairtools split --nproc-in 2 --nproc-out 2 --output-sam whole_mapped.bam

conda deactivate
