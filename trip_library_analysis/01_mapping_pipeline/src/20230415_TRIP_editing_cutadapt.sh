#!/bin/bash
# Move up one level from the src directory to access the fastqs/editing/ directory
parent_dir="../"

# Update the log file path
datum=$(date +"%Y%m%d")
touch $parent_dir$datum"_trip_editing_cutadaptlog.txt"

# Update the loop to look for files in ../fastqs/editing/
for filename in $parent_dir"fastq/editing/"*R1_001.fastq.gz; do
    # Remove the directory path from the filename to create the shortname
    shortname=$(basename $filename | sed 's/_R1_001.fastq.gz//')
    # Update the input file path and output file path for the cutadapt command
    cutadapt -j 0 -g AAGGGCCGGCCACAa --discard-untrimmed -o $parent_dir"fastq/editing/Output/"$shortname"_barcode_5trim.fastq.gz" $filename >> $parent_dir$datum"_trip_editing_cutadaptlog.txt"
done
