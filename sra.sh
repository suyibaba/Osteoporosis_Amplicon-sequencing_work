#!/bin/bash
#SBATCH --job-name=sam_work
#SBATCH --mail-type=ALL
#SBATCH --mail-user=o.akinsuyi@ufl.edu
#SBATCH --ntasks=1
#SBATCH --mem=50gb
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=20
#SBATCH --output=serial_test_%j.log
#SBATCH --account=microbiology-dept
#SBATCH --qos=microbiology-dept-b

module load sra/3.0.8

# Define your SRA accession number
SRA_ACCESSION=" SRR22425876"

# Create a directory to store downloaded data
mkdir -p SRA_data
cd SRA_data



# Download the SRA file using prefetch
prefetch $SRA_ACCESSION

# Convert the SRA file to FASTQ format using fastq-dump
fastq-dump --split-files --outdir ./ $SRA_ACCESSION



