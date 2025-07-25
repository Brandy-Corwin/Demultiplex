#!/usr/bin/env python

import argparse
import gzip
import bioinfo
import matplotlib.pyplot as plt

# Define arguments that the user needs to input
def get_args():
    parser = argparse.ArgumentParser(description="This script generates a per base distribution of quality scores at each position for all reads.")
    parser.add_argument("-f", "--file", help="The path to the input fastq file (str)", required=True, type=str)
    parser.add_argument("-l", "--length", help="The length of the reads (int)", required=True, type=int)
    return parser.parse_args()
	
args = get_args()

def init_list(lst: list, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    for i in range(args.length):
        if value == []:
            lst.append([])
        else:
            lst.append(value)
    return lst

# Create and initialize a list for the quality scores
my_list: list = []
my_list = init_list(my_list)

def populate_list(file: str) -> tuple[list, int]:
    """Sums the quality score at the same postion of each record in the given file.
    Returns list of summed quality scores and total number of lines in the file"""
    qc_sums: list = []
    qc_sums = init_list(qc_sums,args.length)
    num_lines: int = 0
    with gzip.open(file, "r") as fh:
        for line in fh:
            num_lines+=1
            line = line.decode('utf-8').strip('\n')
            if num_lines%4 == 0:
                for index,record in enumerate(line):
                    qc_sums[index] += bioinfo.convert_phred(record)
    return qc_sums,num_lines

# Populate the list
my_list, num_lines = populate_list(args.file)

# Calculate the mean
for index in range(len(my_list)):
    my_list[index] = my_list[index]/(num_lines/4)

# Get read number from file name
read = args.file.split("/")[-1].split(".")[0].split('_')[3]

# Make histogram
fig, ax = plt.subplots(figsize=(8, 4))
ax.bar(range(args.length), my_list)
ax.set_xlabel("Base Pair")
ax.set_ylabel("Average Phred Quality Score")
ax.set_title("Average Quality Score for Each Base Pair Position")
plt.savefig(f"{read}_hist.png")
print(f"created {read}_hist.png")