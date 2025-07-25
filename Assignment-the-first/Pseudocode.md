## Problem
This is my strategy for demultiplexing fastq files and reporting index hopping

## Output
The final output should have 2 * number of properly matched read pair fastq file, as well as 2 hopped fastq files and 2 unknown/low quality index fastq files.  
It should also print out/report the number of read-pairs that match, had hopping, and that had unknown indexes.

## Test files
Made four test fastq files with 2 matching indexes (one index has two sequences), 2 swapped, 1 unknown, and 1 low quality  
Made expected output files of test files

## Algorithm outline
- create argparse arguments to input all fastq names
- open indexes.txt as read and create a dictionary of index_seq and names
- open all four fastq files as read
    - loop through by record accross all four files
    - calculate index quality score
    - check if either index is not in index_dict.keys or has a quality score below decided threshold
        - if yes, have filename variable be "unknown" and increment unknown counter
    - then check if the indexes match
        - if yes, have filename variable be name form index dictionary and increment matched counter
    - else have filename variable be "hopped" and increment hopped counter
        make and populate dictionary of swaps {(index1_name, index2_name): count}
    - wirte seq1 to R1 of filename variable file and seq2 to R2 filename variable file and add indexes to header
- once looped through all input files, print to terminal the counter for unknown, hopped, and mapped

## Functions
using some functions from bioinfo.py

```python
def reverse_complement(letter: str) -> str:
    '''Takes an index and converts it to the reverse complement'''
    return rev_comp
Input: "AAGT"
Expected output: "ACTT"
```