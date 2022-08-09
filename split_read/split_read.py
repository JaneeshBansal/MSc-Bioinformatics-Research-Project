print("START OF CODE", flush=True)

####################################### importing the dependencies ###################################

# import the required dependencies 
import mappy as mp
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
import numpy as np
import csv
import bisect
import time
from operator import itemgetter
import re
print("imported the required dependencies", flush=True)

############################################# mapping to L1 ##################################################

# get the L1 consensus sequence 
reference = mp.Aligner("/data/home/bt211032/scratch/research_project/2022-04-08-mappy_seq_monkk/input/L1_consensus.fa", preset = "sr") 
# reference = mp.Aligner("input/L1_consensus.fa") 

if not reference:
    raise Exception("ERROR: failed to load/build index")
print("got the L1 reference", flush=True)

# extracting the human reference 
human_ref = mp.Aligner("/data/home/bt211032/scratch/research_project/2022-04-08-mappy_seq_monkk/input/hg38.fa", preset = "sr")
# human_ref = mp.Aligner("input/hg38.fa", preset = "sr")

if not human_ref:
    raise Exception("ERROR: failed to load/build index")
print("got the human reference", flush=True)

# reading the fastq reads 
# read_1 = SeqIO.parse("/data/home/bt211032/scratch/research_project/2022-04-08-mappy_seq_monk/input/SRR15503419_1.fastq", "fastq")
# read_2 = SeqIO.parse("/data/home/bt211032/scratch/research_project/2022-04-08-mappy_seq_monk/input/SRR15503419_2.fastq", "fastq")

# read_1 = SeqIO.parse("input/sub_sim-R1.fastq", "fastq")
# read_2 = SeqIO.parse("input/sub_sim-R2.fastq", "fastq")

read_1 = SeqIO.parse("/data/home/bt211032/scratch/research_project/2022-05-09-simulated_data/input_py/sim-R1.fastq", "fastq")
read_2 = SeqIO.parse("/data/home/bt211032/scratch/research_project/2022-05-09-simulated_data/input_py/sim-R2.fastq", "fastq")

t0 = time.time()

split_read_R1 = open("results_super_lenient/split_read_R1.txt", "a")
split_read_R2 = open("results_super_lenient/split_read_R2.txt", "a")

# getting a paired read from the fastq files 
for R1, R2 in zip(read_1, read_2):

#####

# RESTRICTION SITE TRIMMING

    # define the restriction site of interest
    restri_site = "GATCGATC"
    # minimum number of bases required before trimming at 3' end
    min_cut = 5

    # R1
    # if the restriction site is found in the R1 sequence 
    if restri_site in R1.seq:
        R1_seq = str(R1.seq)
        # split the string according to the site
        split_seq = R1_seq.split("GATCGATC", 1)
        # keep the bases before the site 
        R1_seq = str(split_seq[0])
    else:
        R1_seq = str(R1.seq)
        # loop for site with minimum number of bases in 3' end 
        for n in range(len(restri_site)-min_cut, -1, -1):
            # if the site is present in the sequence 
            if restri_site[0:min_cut+n] == R1_seq[-(min_cut+n):]:
                # delete the site at the 3' end
                R1_seq = re.sub(restri_site[0:min_cut+n], "", str(R1_seq))
                break
            else:
                # else leave the string as is
                R1_seq = str(R1_seq)

    # R2
    # if the restriction site is found in the R1 sequence 
    if restri_site in R2.seq:
        R2_seq = str(R2.seq)
        # split the string according to the site
        split_seq = R2_seq.split("GATCGATC", 1)
        # keep the bases before the site 
        R2_seq = str(split_seq[0])
    else:
        R2_seq = str(R2.seq)
        # loop for site with minimum number of bases in 3' end 
        for n in range(len(restri_site)-min_cut, -1, -1):
            # if the site is present in the sequence 
            if restri_site[0:min_cut+n] == R2_seq[-(min_cut+n):]:
                # delete the site at the 3' end
                R2_seq = re.sub(restri_site[0:min_cut+n], "", str(R2_seq))
                break
            else:
                # else leave the string as is
                R2_seq = str(R2_seq)

######
    
    # if the length of the resulting sequences are more than 40b
    if len(R1_seq) >= 40 and len(R2_seq) >= 40:

        # mpaping the read to L1
        hits_1 = list(reference.map(R1_seq))
        hits_2 = list(reference.map(R2_seq))

        # checking if only read 2 mapped to L1
        if len(hits_1) == 0 and len(hits_2) > 0:    

            # extracting the hit when mapping the read to L1
            for hit in reference.map(R2_seq):

                # calculate the length of alignment to the L1 consensus
                L1_len = hit.q_en - hit.q_st
                # calculate sequence length that aligned to hg38
                hg38_len = len(R2_seq) - hit.q_en

                # subset the sequence for only portion mapping to L1 consensus
                L1_seq = R2_seq[hit.q_st : hit.q_en]
                # subset sequence for only portion mapping to hg38
                hg38_seq = R2_seq[hit.q_en :]

                L1_start = hit.q_st

                # map the portion mapping to hg38 to hg38 again and get all the hits
                for hit in human_ref.map(hg38_seq):
                    # store mapping to chr, start and end position 
                    hg38_chr = hit.ctg
                    hg38_start = hit.r_st
                    hg38_end = hit.r_en

                    # save details and write to file 
                    R2_info = [str(len(R2_seq)), str(L1_len), str(hg38_len), str(L1_start), str(hg38_chr), str(hg38_start), str(hg38_end)]
                    split_read_R2.write("\t".join(R2_info))
                    split_read_R2.write("\n")




        # checking if only read 1 mapped to L1
        elif len(hits_1) > 0 and len(hits_2) == 0:

            # extracting the hit when mapping the read to L1
            for hit in reference.map(R1_seq):

                # calculate the length of alignment to the L1 consensus
                L1_len = hit.q_en - hit.q_st
                # calculate sequence length that aligned to hg38
                hg38_len = len(R1_seq) - hit.q_en

                # subset the sequence for only portion mapping to L1 consensus
                L1_seq = R1_seq[hit.q_st : hit.q_en]
                # subset sequence for only portion mapping to hg38
                hg38_seq = R1_seq[hit.q_en :]

                L1_start = hit.q_st

                # map the portion mapping to hg38 to hg38 again and get all the hits
                for hit in human_ref.map(hg38_seq):
                    # store mapping to chr, start and end position 
                    hg38_chr = hit.ctg
                    hg38_start = hit.r_st
                    hg38_end = hit.r_en

                    # save details and write to file 
                    R1_info = [str(len(R1_seq)), str(L1_len), str(hg38_len), str(L1_start), str(hg38_chr), str(hg38_start), str(hg38_end)]
                    split_read_R1.write("\t".join(R1_info))
                    split_read_R1.write("\n")




        # checking no reads mapped to L1
        elif len(hits_1) == 0 and len(hits_2) == 0:
            continue

        # checking both reads mapped to L1
        elif len(hits_1) > 0 and len(hits_2) > 0:
            continue

    else:
        continue 

print("finished split read", flush=True)

split_read_R1.close()

split_read_R2.close()

print("closed files", flush=True)

print("END OF CODE", flush=True)
