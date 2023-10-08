#!/bin/bash
datum=$(date +"%Y%m%d")
# set up log file
touch $datum"_cutadaptlog.txt"

fastq_folder="fastqs/"
output_folder="Output/"

# define number of characters to remove from filename
forremove=14
keep=17

# loop through all R1 files and trim to barcode and target sequence
for filename in ${fastq_folder}*R1_001.fastq.gz; do
		# Extract the desired part of the filename
		shortname="${filename:${forremove}:${keep}}"
        cutadapt -j 0 -g aagtggtagagtag --discard-untrimmed -o "${output_folder}${shortname}_barcode_5trim.fastq.gz" $filename >> $datum"_cutadaptlog.txt"
        cutadapt -j 0 -g GTCATAGGAGTCAC --discard-untrimmed -o "${output_folder}${shortname}_target_5trim.fastq.gz" $filename >> $datum"_cutadaptlog.txt"
done

# loop through all R2 files and trim to protospacer and RTT sequence
for filename in ${fastq_folder}*R2_001.fastq.gz; do
		shortname="${filename:${forremove}:${keep}}"
        cutadapt -j 0 -g gaaacaccG --discard-untrimmed -o "${output_folder}${shortname}_proto_5trim.fastq.gz" $filename >> $datum"_cutadaptlog.txt"
        cutadapt -j 0 -a GTTTCAGA -m 19 -M 19 --discard-untrimmed -o "${output_folder}${shortname}_Protospacer.fastq.gz" "${output_folder}${shortname}_proto_5trim.fastq.gz" >> $datum"_cutadaptlog.txt"
        cutadapt -j 0 -g AGTCGGTGC --discard-untrimmed -o "${output_folder}${shortname}_RTTtrim.fastq.gz" $filename >> $datum"_cutadaptlog.txt"
done
