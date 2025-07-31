#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=bcorwin@uoregon.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=demultiplex
#SBATCH --output=demultiplex_logs/demultiplex_%j.out
#SBATCH --error=demultiplex_logs/demultiplex_%j.err

mkdir fastqs

file1="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
file2="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
file3="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
file4="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"
indexes="/projects/bgmp/shared/2017_sequencing/indexes.txt"

./Demulitplex.py -f1 $file1 -f2 $file2 -f3 $file3 -f4 $file4 -i $indexes -q 34 > output_summary.txt

exit