###############################################################################################################
############################################### CHR REGIONS ###################################################
###############################################################################################################

print("starting the code...", flush=True)
 
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
print("imported the required dependencies", flush=True)

############################################# mapping to L1 ##################################################

# get the L1 consensus sequence 
reference = mp.Aligner("/data/home/bt211032/scratch/research_project/2022-04-05-bedfiles/results/chr_regions.fa") 
# reference = mp.Aligner("input/chr_regions.fa") 

if not reference:
    raise Exception("ERROR: failed to load/build index")
print("got the reference", flush=True)

# reading the fastq reads 
read_1 = SeqIO.parse("/data/home/bt211032/scratch/research_project/2022-04-08-mappy_seq_monk/input/SRR15503419_1.fastq", "fastq")
read_2 = SeqIO.parse("/data/home/bt211032/scratch/research_project/2022-04-08-mappy_seq_monk/input/SRR15503419_2.fastq", "fastq")

# read_1 = SeqIO.parse("input/sub2_SRR15503419_1.fastq", "fastq")
# read_2 = SeqIO.parse("input/sub2_SRR15503419_2.fastq", "fastq")

# producing empty lists
position_mapped = []
seq_unmapped = []
ID = []
chr_map = []

t0 = time.time()

# getting a paired read from the fastq files 
for R1, R2 in zip(read_1, read_2):
    # mpaping the read to L1
    hits_1 = list(reference.map(R1.seq))
    hits_2 = list(reference.map(R2.seq))

    tmp_info_ls = []

    # checking if only read 2 mapped to L1
    if len(hits_1) == 0 and len(hits_2) > 0:
        # mapping the read to the reference 
        for hit in reference.map(R2.seq):
            # continue if the length fo alignment is >= 95% and the mapq is 60
            if ((hit.blen)/(len(R2.seq))) >= 0.95 and hit.mapq >= 60: 
                # mapq, R2 L1 pos map, R1 unmapped seq, R2 ID, R2 chr mapped 
                tmp_info = (hit.mapq, hit.r_st, R1.seq, R2.id, hit.ctg)
                # append the information to the master list 
                tmp_info_ls.append(tmp_info)
            # if thresholds are not met pass this read 
            else:
                pass

        # if the list is empty pass 
        if tmp_info_ls == []:
            pass
        # if there are reads which met the threshold 
        else:
            # get the information for the read has the highest mapq
            info = max(tmp_info_ls, key=itemgetter(0))
            # appending the position where R2 mapped to L1 
            position_mapped.append(info[1])
            # appending the unmapped sequence of read 1
            seq_unmapped.append(info[2])
            # append R2 ID
            ID.append(info[3])
            # append on what chr R2 maps 
            chr_map.append(info[4])


    # checking if only read 1 mapped to L1
    elif len(hits_1) > 0 and len(hits_2) == 0:
        # mapping the sequence to the reference 
        for hit in reference.map(R1.seq):
            # continue if the length fo alignment is >= 95% and the mapq is 60
            if ((hit.blen)/(len(R1.seq))) >= 0.95 and hit.mapq >= 60:
                # mapq, R1 L1 pos map, R2 unmapped seq, R1 id, R1 chr map
                tmp_info = (hit.mapq, hit.r_st, R2.seq, R1.id, hit.ctg)
                # append the information to the master list 
                tmp_info_ls.append(tmp_info)
            # if the threshold requirements are not met pass this read 
            else:
                pass

        # if the list if empty pass 
        if tmp_info_ls == []:
            pass
        # if there are reads which matched the threshold 
        else:
            # max mapq taken 
            info = max(tmp_info_ls, key=itemgetter(0))
            # appending the position where R1 mapped to L1
            position_mapped.append(info[1])
            # appending the unmapped sequence of read 2
            seq_unmapped.append(info[2])
            # R1 id
            ID.append(info[3])
            # R1 chr mapped
            chr_map.append(info[4])

    # checking no reads mapped to L1
    elif len(hits_1) == 0 and len(hits_2) == 0:
        continue

    # checking both reads mapped to L1
    elif len(hits_1) > 0 and len(hits_2) > 0:
        continue

print("have the unmapped reads", flush=True)

t1 = time.time()
print("Time elapsed: ", t1 - t0) # CPU seconds elapsed (floating point)

# joining the normalised position and the sequence unmapped 
pos_seq = zip(ID, chr_map, position_mapped, seq_unmapped)


################################# mapping the unmapped sequence to the human genome #################################

checking_things = []

t0 = time.time()

# extracting the human reference 
human_ref = mp.Aligner("/data/home/bt211032/scratch/research_project/2022-04-08-mappy_seq_monk/input/hg38.fa", preset = "sr")
# human_ref = mp.Aligner("input/hg38.fa", preset = "sr")
if not human_ref:
    raise Exception("ERROR: failed to load/build index")
print("got the human reference", flush=True)

# opening an text file to add the outputs 
output_txt = open("results_mapq_60/chr_mapped_60.txt", "a")

# extracting a read from the list of tuples 
for read in pos_seq:
    # extracting the unmapped sequence
    seq = read[3]
    # empty list to store info for mapq>30
    tmp_info_ls = []

    # mapping the unmapped sequence to the human reference
    for hit in human_ref.map(seq):

        # subsetting if the mapping quality is greater than 30 
        if hit.mapq >= 30:
            # mapq, ID, chr region, chr_start, start hg38, end hg38 and chr 
            tmp_info = (hit.mapq, read[0], read[1], read[2], hit.r_st, hit.r_en, hit.ctg)
            # append the tuple to the list 
            tmp_info_ls.append(tmp_info)

        # mapping qualities lower than 30 are passed 
        else:
            pass

    # when mapq <30 there are no results in list 
    if tmp_info_ls == []:
        pass

    else:
        # getting the max mapq from the list of tuples
        info = max(tmp_info_ls, key=itemgetter(0))

        # producing a list of the required information
        hit_info = [str(info[0]), str(info[1]), str(info[2]), str(info[3]), str(info[4]), str(info[5]), str(info[6])]
        # writing the info in the text file and tab separating it
        output_txt.write("\t".join(hit_info))
        # writing a new line 
        output_txt.write("\n")


# closing the text file 
output_txt.close()

t1 = time.time()
print("Time elapsed: ", t1 - t0) # CPU seconds elapsed (floating point)

print("mapped the unmapped reads and made the text file!", flush=True)

#########################################################################################################################
