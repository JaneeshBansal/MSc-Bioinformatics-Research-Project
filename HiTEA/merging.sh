#!/bin/sh
#$ -cwd              # Set the working directory for the job to the current directory
#$ -j y              # Combine the standard output and the standard error stream
#$ -pe smp 4        # Request 16 cores
#$ -l h_rt=72:0:0    # Request 72 hours runtime
#$ -l h_vmem=4G      # Request 6GB RAM
#$ -m bea            # Send an email once job starts, finishes and aborts

module load samtools

samtools merge mapped_merged.bam *_mapped.bam
