#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH --mail-user=bcorwin@uoregon.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=distribution_hists
#SBATCH --output=distribution_hists_logs/distribution_hists%j.out
#SBATCH --error=distribution_hists_logs/distribution_hists%j.err

path1="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
path2="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
path3="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
path4="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"

./quality_dist.py -f $path1 -l 101
./quality_dist.py -f $path2 -l 101
./quality_dist.py -f $path3 -l 101
./quality_dist.py -f $path4 -l 101

exit