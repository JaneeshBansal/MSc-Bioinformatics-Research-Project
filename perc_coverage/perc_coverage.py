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
reference = mp.Aligner("/data/home/bt211032/scratch/research_project/2022-04-08-mappy_seq_monkk/input/L1_consensus.fa") 
# reference = mp.Aligner("input/L1_consensus.fa") 

if not reference:
    raise Exception("ERROR: failed to load/build index")
print("got the reference", flush=True)

# reading the fastq reads 
# read_1 = SeqIO.parse("/data/home/bt211032/scratch/research_project/2022-04-08-mappy_seq_monk/input/SRR15503419_1.fastq", "fastq")
# read_2 = SeqIO.parse("/data/home/bt211032/scratch/research_project/2022-04-08-mappy_seq_monk/input/SRR15503419_2.fastq", "fastq")

# read_1 = SeqIO.parse("input/sub2_SRR15503419_1.fastq", "fastq")
# read_2 = SeqIO.parse("input/sub2_SRR15503419_2.fastq", "fastq")

read_1 = SeqIO.parse("/data/home/bt211032/scratch/research_project/2022-05-09-simulated_data/input_py/sim-R1.fastq", "fastq")
read_2 = SeqIO.parse("/data/home/bt211032/scratch/research_project/2022-05-09-simulated_data/input_py/sim-R2.fastq", "fastq")

# producing empty lists
position_mapped = []
seq_unmapped = []

L1_pos = [] ###

t0 = time.time()

# getting a paired read from the fastq files 
for R1, R2 in zip(read_1, read_2):
    # mpaping the read to L1
    hits_1 = list(reference.map(R1.seq))
    hits_2 = list(reference.map(R2.seq))

    # checking if only read 2 mapped to L1
    if len(hits_1) == 0 and len(hits_2) > 0:
        # extracting the hit when mapping the read to L1
        for hit in reference.map(R2.seq):

            # subsetting values for lengths of alignment greater than 95%
            if ((hit.blen)/(len(R2.seq))) >= 0.95: 
                # appending the position where R2 mapped to L1 
                position_mapped.append(hit.r_st)
                # appending the unmapped sequence of read 1
                seq_unmapped.append(R1.seq)

                L1_pos.append((R2.id, hit.r_st))

            # if the length of alignment is less than 95% pass this hit
            else:
                pass

    # checking if only read 1 mapped to L1
    elif len(hits_1) > 0 and len(hits_2) == 0:
        # extracting the hit when mapping the read to L1
        for hit in reference.map(R1.seq):

            # subsetting values for lengths of alignment greater than 95%
            if ((hit.blen)/(len(R1.seq))) >= 0.95:
                # appending the position where R1 mapped to L1
                position_mapped.append(hit.r_st)
                # appending the unmapped sequence of read 2
                seq_unmapped.append(R2.seq)

                L1_pos.append((R1.id, hit.r_st))

            # if the length of alignment is less than 95% pass this hit
            else:
                pass

    # checking no reads mapped to L1
    elif len(hits_1) == 0 and len(hits_2) == 0:
        continue

    # checking both reads mapped to L1
    elif len(hits_1) > 0 and len(hits_2) > 0:
        continue

print("have the unmapped reads", flush=True)

t1 = time.time()
print("Time elapsed: ", t1 - t0) # CPU seconds elapsed (floating point)

# opening the L1 fasta file 
L1 = open("/data/home/bt211032/scratch/research_project/2022-04-08-mappy_seq_monkk/input/L1_consensus.fa", "r")
# L1 = open("input/L1_consensus.fa", "r")

# get the length of the sequence in the fasta file  
for name, seq in SimpleFastaParser(L1):
    seqLen = len(seq)
print("sequence length is ...", seqLen)

# closing the file 
L1.close()

# joining the normalised position and the sequence unmapped 
pos_seq = zip(position_mapped, seq_unmapped)

print("zipped L1 start positions mapped", flush=True)


################################# mapping the unmapped sequence to the human genome #################################

t0 = time.time()

# produce empty lists 
L1_start = []
human_start = []
human_chr = []

# extracting the human reference 
human_ref = mp.Aligner("/data/home/bt211032/scratch/research_project/2022-04-08-mappy_seq_monkk/input/hg38.fa", preset = "sr")
# human_ref = mp.Aligner("input/hg38.fa", preset = "sr")

if not human_ref:
    raise Exception("ERROR: failed to load/build index")
print("got the human reference", flush=True)

# opening an text file to add the outputs 
output_txt = open("results/mapped.txt", "a")

# extracting a read from the list of tuples 
for read in pos_seq:
    # extracting the unmapped sequence
    seq = read[1]
    # empty list to store info for mapq>30
    tmp_info_ls = []

    # mapping the unmapped sequence to the human reference
    for hit in human_ref.map(seq):

        # subsetting if the mapping quality is greater than 30 
        if hit.mapq >= 30:
            # mapq, L1_start, start hg38, end hg38 and chr 
            tmp_info = (hit.mapq, read[0], hit.r_st, hit.r_en, hit.ctg)
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

        # appending the L1 mapping start position
        L1_start.append(info[1])
        # appending the start position of mapping to the hg38
        human_start.append(info[2])
        # appending what chromosome it maps to 
        human_chr.append(info[4])
        # producing a list of the required information
        hit_info = [str(info[1]), str(info[4]), str(info[2]), str(info[3])]
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

# make a list of tuple values with the lists produced 
L1_human_start = list(map(lambda x,y,z:(x,y,z), L1_start, human_chr, human_start))

# adding the L1 position and hg38 position with the chr as the key in the dictionary 
mapping_dict = {}
# extracting from the tuple
for L1_start, human_chr, human_start in L1_human_start:
    # adding to a chr key if already present 
    if human_chr in mapping_dict:
        mapping_dict[human_chr].append((L1_start, human_start))
    # adding a new key to the dictionary if the chr is not already present and adding the positions
    else:
        mapping_dict[human_chr] = [(L1_start, human_start)]

print("got a dictionary for the L1 start and hg38 start position mapped", flush=True)

# key is set to sort by the second element (first hg38 position mapped)
for chromosome, values in mapping_dict.items():
    mapping_dict[chromosome] = sorted(values, key = lambda x: x[1])

print("sorted the dictionary according to hg38 start position mapped", flush=True)

################################## SW ####################################

t0 = time.time()

L1_human_dict = {}
for chromosome, values in mapping_dict.items():
    L1_start, human_start = zip(*mapping_dict[chromosome])
    L1_human_dict[chromosome] = [list(human_start), list(L1_start)]

#### user can select ####
window = 20000
step = 1000
#########################

# empty dictionary to store the variance values 
all_coverage_dict = {}

# extracting information from the dictionary 
for chromosome, values in L1_human_dict.items():
    # indexing the hg38 first start position 
    start = values[0][0]
    # reassigning to the lower bound
    lower_bound = start
    # addition of the window to make the uppper bound 
    upper_bound = lower_bound + window

    # continue the loop until the upper bound exceeds the last hg38 position recorded
    while upper_bound < values[0][-1]:
        # finding the index for the lower bound 
        lower_bound_i = bisect.bisect_left(values[0], lower_bound)
        # finding the index for the upper bound 
        upper_bound_i = bisect.bisect_right(values[0], upper_bound, lo=lower_bound_i)
        # extracting values within this window for the L1 start values 
        L1_start_ls = values[1][lower_bound_i:upper_bound_i]

        if len(L1_start_ls) >= 1:

            L1_start_ls = sorted(L1_start_ls)

            coverage_dict = {"<1":0, "1-2":0, "2-3":0, "3-4":0, "4-5":0, "5<":0}

            coverage_dict_keys = list(coverage_dict.keys())

            ### 
            mini_window = 1000
            mini_step = 1000
            ###

            mini_start = 0

            mini_lower_bound = mini_start
            mini_upper_bound = mini_window

            while mini_upper_bound < seqLen:
                
                for n in range(0, len(coverage_dict_keys)):
                    
                    if n < (len(coverage_dict_keys)-1):

                        mini_lower_bound_i = bisect.bisect_left(L1_start_ls, mini_lower_bound)
                        mini_upper_bound_i = bisect.bisect_right(L1_start_ls, mini_upper_bound, lo=mini_lower_bound_i)
                        L1_start_subset = L1_start_ls[mini_lower_bound_i:mini_upper_bound_i]
                    
                    else:
                        mini_lower_bound_i = bisect.bisect_left(L1_start_ls, mini_lower_bound)
                        L1_start_subset = L1_start_ls[mini_lower_bound_i:]

                    coverage_dict[coverage_dict_keys[n]] = len(L1_start_subset)

                    mini_upper_bound += mini_step
                    mini_lower_bound += mini_step
                    n += 1

            # percentage coverage
            tot_covered = sum(value != 0 for value in coverage_dict.values())
            tot = len(coverage_dict_keys)
            perc_coverage = (tot_covered/tot)*100
            
            if chromosome in all_coverage_dict:
                all_coverage_dict[chromosome].append((perc_coverage, lower_bound, upper_bound, L1_start_ls, len(L1_start_ls)))
            # addition of information for a new chr record 
            else:
                all_coverage_dict[chromosome] = [(perc_coverage, lower_bound, upper_bound, L1_start_ls, len(L1_start_ls))]

        # addition of a step to the windows
        lower_bound += step
        upper_bound += step

# print(var_dict)
        
t1 = time.time()

print("Time elapsed: ", t1 - t0) # CPU seconds elapsed (floating point)
        
print("got a dictionary for L1 mapping percentage coverage", flush=True)

####################################################################################
            
# writing the tsv file 
with open("results/window.tsv", "wt") as out_file:
    # producing the writer which separates values by tabs 
    tsv_writer = csv.writer(out_file, delimiter='\t')
    # writing the header
    tsv_writer.writerow(['chr', 'perc_coverage', 'start_hg38', 'end_hg38', "L1_start", "num_L1_start"])
    # extracting the information from the variance dictionary 
    for key, value in all_coverage_dict.items():
        # looping over the multiple results for a chromosome 
        for n in range(len(value)):
            # writing the chromosome, variance, start hg38 and end hg38 position mapped 
            tsv_writer.writerow((key, value[n][0], value[n][1], value[n][2], value[n][3], value[n][4]))

# closing the output file 
out_file.close()

print("produced the tsv file", flush=True)

