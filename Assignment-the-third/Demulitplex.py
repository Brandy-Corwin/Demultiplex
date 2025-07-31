#!/usr/bin/env python

import argparse
import gzip
import bioinfo

# Define arguments that the user needs to input
def get_args():
    parser = argparse.ArgumentParser(description="This script demultiplexes given gzipped fastq files and finds the number of matched, unknown, and swapped reads.")
    parser.add_argument("-f1", "--file1", help="The path to the R1 input gzipped fastq file (str)", required=True, type=str)
    parser.add_argument("-f2", "--file2", help="The path to the R2 input gzipped fastq file (str)", required=True, type=str)
    parser.add_argument("-f3", "--file3", help="The path to the R3 input gzipped fastq file (str)", required=True, type=str)
    parser.add_argument("-f4", "--file4", help="The path to the R4 input gzipped fastq file (str)", required=True, type=str)
    parser.add_argument("-i", "--indexes", help="The path to the indexes file (str)", required=True, type=str)
    parser.add_argument("-q", "--qscore", help="The qscore cutoff for low quality indexes (int)", required=True, type=int)
    return parser.parse_args()
	
args = get_args()

# Create reverse complement function
def reverse_complement(seq: str) -> str:
    '''Takes a sequence and converts it to the reverse complement'''
    rev_comp = ""
    for index in range(len(seq)-1, -1, -1):
        if seq[index] == "A":
            letter = "T"
        elif seq[index] == "T":
            letter = "A"
        elif seq[index] == "C":
            letter = "G"
        elif seq[index] == "G":
            letter = "C"
        elif seq[index] == "N":
            letter = "N"
        rev_comp = rev_comp + letter
    return rev_comp
assert reverse_complement("AAGT") == "ACTT", "Wrong reverse compliment for ACTT"

# Open all read files
with (open(args.indexes, "r") as index_file, gzip.open(args.file1, "r") as R1, gzip.open(args.file2, "r") as R2,
      gzip.open(args.file3, "r") as R3, gzip.open(args.file4, "r") as R4):
    
    # Open all write files and create set of indexes
    files_dict = {}
    index_set = set()
    for line in index_file:
        if "index sequence" in line:
            continue
        index = line.strip('\n').split('\t')[-1]
        index_set.add(index)
        files_dict[index] = (open("fastqs/" + index + "_R1.fastq", "w"), open("fastqs/" + index + "_R2.fastq", "w"))
    files_dict["unknown"] = (open("fastqs/unknown_R1.fastq", "w"), open("fastqs/unknown_R2.fastq", "w"))
    files_dict["swapped"] = (open("fastqs/swapped_R1.fastq", "w"), open("fastqs/swapped_R2.fastq", "w"))

    # Initiate counts, swaps dictionary, and matched dictionary
    num_matched = 0
    num_swapped = 0
    num_unknown = 0
    swaps = {}
    matched = {}

    # Loop through input fastq files by record
    while True:
        
        R1_header = R1.readline().decode('utf-8')

        # Stop loop if at end of file
        if R1_header == '':
            break

        # record for R1 file
        R1_header = R1_header.strip('\n')
        R1_seq = R1.readline().decode('utf-8').strip('\n')
        R1_plus = R1.readline().decode('utf-8').strip('\n')
        R1_qscore = R1.readline().decode('utf-8').strip('\n')

        # record for R2 file
        R2_header = R2.readline().decode('utf-8').strip('\n')
        R2_index = R2.readline().decode('utf-8').strip('\n')
        R2_plus = R2.readline().decode('utf-8').strip('\n')
        R2_qscore = R2.readline().decode('utf-8').strip('\n')

        # record for R3 file
        R3_header = R3.readline().decode('utf-8').strip('\n')
        R3_index = R3.readline().decode('utf-8').strip('\n')
        R3_plus = R3.readline().decode('utf-8').strip('\n')
        R3_qscore = R3.readline().decode('utf-8').strip('\n')
        
        # record for R4 file
        R4_header = R4.readline().decode('utf-8').strip('\n')
        R4_seq = R4.readline().decode('utf-8').strip('\n')
        R4_plus = R4.readline().decode('utf-8').strip('\n')
        R4_qscore = R4.readline().decode('utf-8').strip('\n')

        # get average quality score for the two indexes
        R2_avg_qscore = bioinfo.qual_score(R2_qscore)
        R3_avg_qscore = bioinfo.qual_score(R3_qscore)
        
        # get the reverse complement of the R3 index
        R3_RC_index = reverse_complement(R3_index)

        if R2_index not in index_set or R3_RC_index not in index_set or R2_avg_qscore < args.qscore or R3_avg_qscore < args.qscore:
            files_dict["unknown"][0].write(f"{R1_header} {R2_index}-{R3_RC_index}\n")
            files_dict["unknown"][0].write(f"{R1_seq}\n")
            files_dict["unknown"][0].write(f"{R1_plus}\n")
            files_dict["unknown"][0].write(f"{R1_qscore}\n")

            files_dict["unknown"][1].write(f"{R4_header} {R2_index}-{R3_RC_index}\n")
            files_dict["unknown"][1].write(f"{R4_seq}\n")
            files_dict["unknown"][1].write(f"{R4_plus}\n")
            files_dict["unknown"][1].write(f"{R4_qscore}\n")

            num_unknown += 1

        elif R2_index != R3_RC_index:
            files_dict["swapped"][0].write(f"{R1_header} {R2_index}-{R3_RC_index}\n")
            files_dict["swapped"][0].write(f"{R1_seq}\n")
            files_dict["swapped"][0].write(f"{R1_plus}\n")
            files_dict["swapped"][0].write(f"{R1_qscore}\n")

            files_dict["swapped"][1].write(f"{R4_header} {R2_index}-{R3_RC_index}\n")
            files_dict["swapped"][1].write(f"{R4_seq}\n")
            files_dict["swapped"][1].write(f"{R4_plus}\n")
            files_dict["swapped"][1].write(f"{R4_qscore}\n")

            num_swapped +=1
            if (R2_index, R3_RC_index) not in swaps:
                swaps[(R2_index, R3_RC_index)] = 1
            elif (R2_index, R3_RC_index) in swaps:
                swaps[(R2_index, R3_RC_index)] += 1

        elif R2_index == R3_RC_index:
            files_dict[R2_index][0].write(f"{R1_header} {R2_index}-{R3_RC_index}\n")
            files_dict[R2_index][0].write(f"{R1_seq}\n")
            files_dict[R2_index][0].write(f"{R1_plus}\n")
            files_dict[R2_index][0].write(f"{R1_qscore}\n")

            files_dict[R2_index][1].write(f"{R4_header} {R2_index}-{R3_RC_index}\n")
            files_dict[R2_index][1].write(f"{R4_seq}\n")
            files_dict[R2_index][1].write(f"{R4_plus}\n")
            files_dict[R2_index][1].write(f"{R4_qscore}\n")

            num_matched += 1
            if R2_index not in matched:
                matched[R2_index] = 1
            elif R2_index in matched:
                matched[R2_index] += 1

    # close all write files
    for index in files_dict:
        files_dict[index][0].close()
        files_dict[index][1].close()
    files_dict["swapped"][0].close()
    files_dict["swapped"][1].close()
    files_dict["unknown"][0].close()
    files_dict["unknown"][1].close()

print(f"There were {num_matched} matched read-pairs ({num_matched/(num_matched+num_swapped+num_unknown)*100}% of total reads)")
print(f"There were {num_swapped} read-pairs that had index-hopping")
print(f"There were {num_unknown} read-pairs that either had low quality indexes or unknown indexes")
print("\n")

print(f"Index1,Index2 Swap Pairs\tCount")
for swap in swaps:
    print(f"{swap}\t{swaps[swap]}")
print("\n")

print(f"Index Matched Pairs\tCount")
for match in matched:
    print(f"{match}\t{matched[match]}")