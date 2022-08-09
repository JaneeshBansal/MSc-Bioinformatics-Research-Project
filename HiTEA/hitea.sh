#!/bin/sh
#$ -cwd
#$ -V
#$ -pe smp 4
#$ -l h_rt=72:0:0
#$ -l h_vmem=4G
#$ -m bea


module load perl
module load R
module load bedtools
module load samtools
module load parallel
module load use.own
module load hitea
module load pairtools
module load bwa

# input bam file, output directory, output prefix, genome built to be used, digestive enzyme used, remap unmapped clipped reads
hitea -i /data/home/bt211032/scratch/research_project/2022-04-13-bwa_pairtools/mapped_merged.bam -w /data/home/bt211032/scratch/research_project/2022-04-07-hitea_whole_sorted/results_3 -o whole -g hg38 -e 'DpnII' -r 'T'
