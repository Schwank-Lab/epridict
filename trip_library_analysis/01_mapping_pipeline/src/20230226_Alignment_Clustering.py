import argparse
import pandas as pd

# Define argument parser
parser = argparse.ArgumentParser(description='Count occurrence of full matches in bowtie alignments')
parser.add_argument('input', type=str, help='Input CSV file of alignments')
parser.add_argument('output', type=str, help='Output CSV file of full match occurrence')
#parser.add_argument('--drop_percentage', type=float, default=0.33, help='Drop percentage for detecting sudden drop in occurrence')
parser.add_argument('--min_occurrence', type=int, default=50, help='Minimum occurrence number to keep')
parser.add_argument('--min_mapq', type=int, default=30, help='Minimum MAPQ score to keep alignment')
args = parser.parse_args()

# Read the CSV file
df = pd.read_csv("results/alignment/"+args.input)


# Filter for full alignments and get the start locations
full_alignments = df[df['CIGAR'].str.contains('^(\d+)M$', regex=True)].copy()

# only keep alignments that have better MAPQ score than defined minimum
full_alignments = full_alignments[full_alignments['MAPQ'] > args.min_mapq].copy()


# Set correct start_position (some reads are mapped in reverse complement to the genome):
full_alignments['start_location'] = full_alignments.apply(lambda x: x.start_location + int(x.CIGAR[:-1]) if x.flag == 16 else x.start_location ,axis=1)

# Only keep start location and chromosome in dataframe
full_alignments = full_alignments[['start_location', 'chromosome']]


# Count the occurrence of each chromosome and start location combination
occurrence = full_alignments.groupby(['chromosome', 'start_location']).size().reset_index(name='occurrence')

# Sort the occurrence dataframe by decreasing order of occurrence
occurrence = occurrence.sort_values('occurrence', ascending=False)


### only uncomment this part if there you want to filter for a drop_percentage (e.g. for clones with only few barcodes)
# Find the first index where there is a sudden drop in occurrence
#for i in range(1, len(occurrence)):
#    if occurrence['occurrence'].iloc[i] < (args.drop_percentage * occurrence['occurrence'].iloc[i-1]):
#        drop_index = i
#        break
#else:
#    drop_index = len(occurrence)

# Filter the dataframe to keep only rows above the drop index
#filtered_df = occurrence.iloc[:drop_index]
filtered_df = occurrence

# Write the full match occurrence to a new CSV file
#filtered_df.to_csv("bowtie_alignment/"+args.output, index=False)


### cluster close occurrences together
# Create a copy of filtered_df and reset the index
merged_df = filtered_df.copy().reset_index(drop=True)

# Loop through all rows from top to bottom
for i, row in merged_df.iterrows():
    # Check if the occurrence value is less than the minimum specified, if so break the loop
    if row['occurrence'] < args.min_occurrence:
        break
    # Check if there are other rows with the same chromosome within +/-10 of start_location
    close_rows = merged_df[(merged_df['chromosome'] == row['chromosome']) &
                           (abs(merged_df['start_location'] - row['start_location']) <= 10) &
                           (merged_df.index > i)]
    if len(close_rows) > 0:
        # Add the occurrence of the close rows to the current row and remove the close rows
        merged_occurrence = row['occurrence'] + close_rows['occurrence'].sum()
        merged_df.loc[i, 'occurrence'] = merged_occurrence
        merged_df = merged_df.drop(close_rows.index)

# Sort the merged_df by decreasing order of occurrence
merged_df = merged_df.sort_values('occurrence', ascending=False)

# drop alignments with less than the minimum specified occurrences
merged_df = merged_df[merged_df['occurrence'] >= args.min_occurrence]

merged_df.to_csv("results/alignment/"+args.output, index=False)
