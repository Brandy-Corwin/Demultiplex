# Start Demultiplex Assignment (07/24/2025)
## Assignment the first

### Data exploration
bash commands used:  
*Commands for R1 are shown, but were used on all 4 files*  
```bash
cd /projects/bgmp/shared/2017_sequencing
zcat 1294_S1_L008_R1_001.fastq.gz | head
zcat 1294_S1_L008_R1_001.fastq.gz | head -2 | tail -1 | wc
```

### Part 2
- made Pseudocode.md to write out pseudocode  
- made test fastq files called test_R<1,2,3,4>.fastq  
- Made expected output files of test files
- Made test_indexes.txt
- Wrote out pseudocode to demultiplex the fastq files
- Wrote outline of functions I will use

### Distribution of quality scores
- made quality_dist.py
- made make_hists.sh to run quality_dist.py for each fastq file
- Currently running

- used ```zcat 1294_S1_L008_R2_001.fastq.gz 1294_S1_L008_R3_001.fastq.gz | sed -n '2~4p' | grep -E "N" | wc -l``` to get the number of indexes that had base N  
    answer was 7304664

**(07/25/2025)**
- Finished running histograms - look good
    - picked index cutoff as 34 and bio read pairs as 36

## Assignment the third
### Writing the code
**(07/31/2025)**
- Created Demultiplex.py
    - imported modules and wrote out argparse
    - opened all read and write files
    - using a while true loop to go through the four input fastqs by record
        - all fastqs are same length so if one file is done, the rest are, so only have to check if one file is finished
    - created reverse_complement function and checked it, it works
    - initialized counts and a dictionary to keep track of the amount each index swap
    - wrote if statements to check if read-pairs are either unknown, swapped, or matched
    - Wrote print statements to decribe output
    - added heatmap for swaps and bar graph for matched pairs

- Ran correctly on my test data!
    - used diff to check output files with expected output files

- Created run_demultiplex.sh to run my python script on sbatch

- Running on full data: data good, but xticks on plots are messed up
    - Fixed it

### Done!