import argparse
import pandas as pd
import re
from tqdm import tqdm
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator

# Define command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("for_file", help="Input file for For alignment clustering")
parser.add_argument("rev_file", help="Input file for Rev alignment clustering")
parser.add_argument("forbc_file", help="Input file for ForBC tagmentation output")
parser.add_argument("barcode_file", help="Input file for trimmed barcodes")
parser.add_argument("output_file", help="Output file for tagged ForBC reads")
parser.add_argument("fastq_location", help="Directory for fastq files")
parser.add_argument("alignment_dir", help="Directory for input/output files")
parser.add_argument("--tolerance", type=int, default=1000, help="Tolerance for matching start locations")
parser.add_argument('--min_mapq', type=int, default=15, help='Minimum MAPQ score to keep alignment')
args = parser.parse_args()

# Read the three CSV files
df_for = pd.read_csv(f"{args.alignment_dir}/{args.for_file}")
df_rev = pd.read_csv(f"{args.alignment_dir}/{args.rev_file}")
df_forbc = pd.read_csv(f"{args.alignment_dir}/{args.forbc_file}")

# Read the barcode file and create a dictionary with read name as key and sequence as value
barcode_dict = {}
with gzip.open(args.fastq_location+"/trimmed/"+args.barcode_file, "rt") as f:
    for title, seq, qual in tqdm(FastqGeneralIterator(f), desc="Loading barcodes"):
        read_name = title.split(" ")[0]
        barcode_dict[read_name] = seq

# df_for = pd.read_csv("NM-TRIP_Library_For_alignment_summary.csv")
# df_rev = pd.read_csv("NM-TRIP_Library_Rev_alignment_summary.csv")
# df_forbc = pd.read_csv("NM-TRIP_Library_ForBC_tagmentation_output.csv")

#df_forbc = df_forbc[:200000]

# Define function to check for matches
def check_for_matches(row):
    # If CIGAR contains anything other than "M", return None
    if not re.match(r"^\d+M$", row["CIGAR"]):
        return None
    # if mapping quality is lower than this threshold, do not map read
    if row['MAPQ'] < args.min_mapq:
        return None
    # Otherwise, check if there are any matches within 1000 bp and on the same chromosome
    df_for_matches = df_for[(df_for["chromosome"] == row["chromosome"]) & 
                            (abs(df_for["start_location"] - row["start_location"]) <= 1000)]
    df_rev_matches = df_rev[(df_rev["chromosome"] == row["chromosome"]) & 
                            (abs(df_rev["start_location"] - row["start_location"]) <= 1000)]
    # If there are no matches, return None
    if df_for_matches.empty or df_rev_matches.empty:
        return None
    # If there are matches, add start positions to the row and return it
    result = row.copy()
    for i, match in enumerate(df_for_matches["start_location"]):
        result[f"for_start_position_{i+1}"] = match
    for i, match in enumerate(df_rev_matches["start_location"]):
        result[f"rev_start_position_{i+1}"] = match
    return result

# Apply the function to the forbc dataframe with tqdm progress bar
result_rows = []
for index, row in tqdm(df_forbc.iterrows(), total=df_forbc.shape[0]):
    row_result = check_for_matches(row)
    if row_result is not None:
        row_result['barcode'] = barcode_dict[row.read_name]  # add barcode
        result_rows.append(row_result)

df_forbc_temp = pd.DataFrame(result_rows)


# Drop rows with both for_start_position_2 and rev_start_position_2 values
try:  # only do this if more than one for_start_position_2 and rev_start_position_2 values exist (might not be the case for clone analysis)
    df_forbc_final = df_forbc_temp.drop(df_forbc_temp[(df_forbc_temp["for_start_position_2"].notnull()) & (df_forbc_temp["rev_start_position_2"].notnull())].index)
except KeyError:
    df_forbc_final = df_forbc_temp

# Save the resulting dataframe
df_forbc_final.to_csv(f"{args.alignment_dir}/{args.output_file}", index=False)