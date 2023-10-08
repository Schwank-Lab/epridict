#!/bin/bash

### Shell script to run the TRIP library analysis pipeline. Written by Nicolas Mathis (@nicopmat). ###

# Declare the overview CSV file name, containing information about the samples
overview_filename="20230410_library_overview_file.csv"

# pipeline suffix for running pipeline with different settings
pipsuffix="Library_finalmapping"

# Set cutadapt parameters
adapter_rev="cgtcaattttacgcagactatctttctagggttaa"
adapter_for="gtacgtcacaatatgattatctttctagggttaa"
adapter_forbc="agttggtaccgatca"  # crop before barcode on ForBC_R2 read
adapter_forbc_reverse="ttaaccctagaaagataat"  # crop forBC reads if they go into piggybac repeats
adapter_forbc_shortread="GGGGGGGGGGGGGGGX"  # remove reads where read is shorter than NGS read length (leading to 3' GGGGG)
adapter_forbc_reporterread="agtcacaaggcctgcaggtcaca"  # remove reads where reporter sequence is in amplicon which is further away from genomic sequence than 150bp
barcode_length=16
forbc_mapping_length=40  # bases to map for forBC alignment; only 100% matches will be used for mapping

# Set bowtie2 parameters
bowtie2_options="--threads 32 --no-unal"

genome_version="hg38"  # define genome version (hg19, hg38 or chm13v2)

index_location="/mnt/TRIP/DSB_TRIP_protocol/index/${genome_version}"

alignment_dir="results/alignment"  # relative to shell script location

# Set filtering parameters
mapq_limit_for_rev=30 # only use for/rev mapping if mapq score is > 30
mapq_limit_forbc=30 # only use forBC mapping if mapq score is > 30 (could be less stringent than for/rev because forBC will only map if for/rev are in close proximity)
drop_percentage=0.33  # if barcode occurrence is less than this value of the previous row, remove all rows below (to remove noise)
min_occurrence=10  # drop alignments in For and Rev if their occurrence is less than this value (remove noise)

percentage_alignment=50  # percentage of occurrence of a specific alignment for each barcode; barcodes where the most occuring match is less than this value will be removed from analysis
percentage_second_alignment=10 # max. percentage of second most common alignment; if the second-most alignment is more than this percentage, remove barcode from analysis
min_mapping_count=10 # minimum number of mappings for each barcode. Lower value leads to more mappings, but also more noise

# Set log file names
datum=$(date +"%Y%m%d")
log_rev="${datum}_${genome_version}_${pipsuffix}_rev_log.txt"
log_for="${datum}_${genome_version}_${pipsuffix}_for_log.txt"
log_forbc_barcode="${datum}_${genome_version}_${pipsuffix}_forBC_barcode_log.txt"
log_forbc_R1="${datum}_${genome_version}_${pipsuffix}_forBC_R1_log.txt"

# Create log files
touch "log/${log_rev}"
touch "log/${log_for}"
touch "log/${log_forbc_barcode}"
touch "log/${log_forbc_R1}"

# Read the overview CSV file and skip the header line
input_data=$(tail -n +2 "input/$overview_filename")


# Loop through files and run cutadapt with specified parameters
# Process Rev files
echo "Processing Rev files..."
echo -e "$input_data" | while IFS=',' read -r name prefix pcr_type filename_r1 filename_r2 folder
do
  if [ "$pcr_type" == "Rev" ]; then
    shortname="$name"
    file="$filename_r2"
    fastq_location="$folder"
    echo "Processing $shortname Rev R2 file: $fastq_location/$file"

	cutadapt -j 0 -g "$adapter_rev" -o "${fastq_location}/trimmed/${shortname}_${pipsuffix}_5trim.fastq.gz" "${fastq_location}/${file}" >> "log/${log_rev}"
	bowtie2 -x "$index_location" -U "${fastq_location}/trimmed/${shortname}_${pipsuffix}_5trim.fastq.gz" -S "${alignment_dir}/${shortname}_${genome_version}_${pipsuffix}_tagmentation_output.sam" "$bowtie2_options" >> "log/${log_rev}" 2>&1
	(echo "read_name,chromosome,start_location,CIGAR,flag,MAPQ" && samtools view -h "${alignment_dir}/${shortname}_${genome_version}_${pipsuffix}_tagmentation_output.sam" | awk '{if($1 !~ /@(SQ|HD|PG)/) print $1 "," $3 "," $4 "," $6 "," $2 "," $5}') > "${alignment_dir}/${shortname}_${genome_version}_${pipsuffix}_tagmentation_output.csv"
	python3 src/20230226_Alignment_Clustering.py "${shortname}_${genome_version}_${pipsuffix}_tagmentation_output.csv" "${shortname}_Rev_${genome_version}_${pipsuffix}_alignment_summary.csv" --min_occurrence $min_occurrence --min_mapq $mapq_limit_for_rev
	fi
done

# Process For files
echo "Processing For files..."
echo -e "$input_data" | while IFS=',' read -r name prefix pcr_type filename_r1 filename_r2 folder
do
  if [ "$pcr_type" == "For" ]; then
    shortname="$name"
    file="$filename_r2"
    fastq_location="$folder"
    echo "Processing $shortname For R2 file: $fastq_location/$file"
    
	cutadapt -j 0 -g "$adapter_for" -o "${fastq_location}/trimmed/${shortname}_${pipsuffix}_5trim.fastq.gz" "${fastq_location}/${file}" >> "log/${log_for}"
	bowtie2 -x "$index_location" -U "${fastq_location}/trimmed/${shortname}_${pipsuffix}_5trim.fastq.gz" -S "${alignment_dir}/${shortname}_${genome_version}_${pipsuffix}_tagmentation_output.sam" "$bowtie2_options" >> "log/${log_for}" 2>&1
	(echo "read_name,chromosome,start_location,CIGAR,flag,MAPQ" && samtools view -h "${alignment_dir}/${shortname}_${genome_version}_${pipsuffix}_tagmentation_output.sam" | awk '{if($1 !~ /@(SQ|HD|PG)/) print $1 "," $3 "," $4 "," $6 "," $2 "," $5}') > "${alignment_dir}/${shortname}_${genome_version}_${pipsuffix}_tagmentation_output.csv"
	python3 src/20230226_Alignment_Clustering.py "${shortname}_${genome_version}_${pipsuffix}_tagmentation_output.csv" "${shortname}_For_${genome_version}_${pipsuffix}_alignment_summary.csv" --min_occurrence $min_occurrence --min_mapq $mapq_limit_for_rev
	fi
done

# Process ForBC files
echo "Processing ForBC files..."
echo -e "$input_data" | while IFS=',' read -r name prefix pcr_type filename_r1 filename_r2 folder
do
  if [ "$pcr_type" == "ForBC" ]; then
    shortname="$name"
    file_r1="$filename_r1"
    file_r2="$filename_r2"
    fastq_location="$folder"
    
    echo "Processing $shortname ForBC R2 file: $fastq_location/$file_r2"
	cutadapt -j 0 -g "$adapter_forbc" -l "$barcode_length" -o "${fastq_location}/trimmed/${shortname}_ForBC_${pipsuffix}_barcode_trim.fastq.gz" "${fastq_location}/${file_r2}" >> "log/${log_forbc_barcode}"


	echo "Processing $shortname ForBC R1 file: $fastq_location/$file_r1"
	# remove reads that are too short so they have GGGGGG at the end of the read
	cutadapt -j 0 -a "$adapter_forbc_shortread" --discard-trimmed -o "${fastq_location}/trimmed/temp_${shortname}_${pipsuffix}_R1_trim.fastq.gz" "${fastq_location}/${file_r1}" >> "log/${log_forbc_R1}"
	# remove reads that are too short and have only the reporter itself in the mapping read
	cutadapt -j 0 -a "$adapter_forbc_reporterread" --discard-trimmed -o "${fastq_location}/trimmed/temp2_${shortname}_${pipsuffix}_R1_trim.fastq.gz" "${fastq_location}/trimmed/temp_${shortname}_${pipsuffix}_R1_trim.fastq.gz" >> "log/${log_forbc_R1}"
	# trim reads that have the start of the reporter at this position to prevent mapping of the reporter sequence
	cutadapt -j 0 -a "$adapter_forbc_reverse" -o "${fastq_location}/trimmed/${shortname}_${pipsuffix}_R1_trim.fastq.gz" "${fastq_location}/trimmed/temp2_${shortname}_${pipsuffix}_R1_trim.fastq.gz" >> "log/${log_forbc_R1}"
	
	bowtie2 -x "$index_location" -U "${fastq_location}/trimmed/${shortname}_${pipsuffix}_R1_trim.fastq.gz" -S "${alignment_dir}/${shortname}_${genome_version}_${pipsuffix}_tagmentation_output.sam" "$bowtie2_options" >> "log/${log_forbc_R1}" 2>&1
	(echo "read_name,chromosome,start_location,CIGAR,flag,MAPQ" && samtools view -h "${alignment_dir}/${shortname}_${genome_version}_${pipsuffix}_tagmentation_output.sam" | awk '{if($1 !~ /@(SQ|HD|PG)/) print $1 "," $3 "," $4 "," $6 "," $2 "," $5}') > "${alignment_dir}/${shortname}_ForBC_${genome_version}_${pipsuffix}_tagmentation_output.csv"
	
	
	# Run alignment-checker on the ForBC alignments and check whether alignment can be found within 1000bp of the For and Rev location
	python src/20230226_Barcode_Alignment_Checker_For_Rev.py "${shortname}_For_${genome_version}_${pipsuffix}_alignment_summary.csv" "${shortname}_Rev_${genome_version}_${pipsuffix}_alignment_summary.csv" "${shortname}_ForBC_${genome_version}_${pipsuffix}_tagmentation_output.csv" "${shortname}_ForBC_${pipsuffix}_barcode_trim.fastq.gz" "${shortname}_ForBC_${genome_version}_${pipsuffix}_tagmentation_with_checks_output.csv" "$fastq_location" "$alignment_dir" --tolerance 1000 --min_mapq $mapq_limit_forbc 
	
	# Assign barcode to alignments if majority mapping is not ambiguous
	python src/20230219_Barcode_Assignment_to_Alignments.py "${shortname}_ForBC_${genome_version}_${pipsuffix}_tagmentation_with_checks_output.csv" "${shortname}_${genome_version}_${pipsuffix}_mappings_with_barcodes.csv" "$alignment_dir" --percentage_alignment "$percentage_alignment" --percentage_second_alignment "$percentage_second_alignment" --min_count "$min_mapping_count"
  fi
done

python src/20230219_alignment_to_editing_mapping.py --editingpath "results/editingPCR/" --editingfile "20230106_TRIP_PE_corrected_editing_values.csv" --finalname "TRIP_PE_editing_with_mapping_${genome_version}_${pipsuffix}_fullpipeline.csv" --mappingpath "results/alignment/" --mappingfile "NM-TRIP_Library_${genome_version}_${pipsuffix}_mappings_with_barcodes.csv"
