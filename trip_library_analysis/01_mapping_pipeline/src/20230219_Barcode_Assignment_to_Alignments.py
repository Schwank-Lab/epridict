# -*- coding: utf-8 -*-
"""
Created on Sat Feb 18 18:47:02 2023

@author: nicol
"""

from collections import Counter
import pandas as pd
from tqdm import tqdm
import argparse
from Bio.Seq import Seq

# Define argparse
parser = argparse.ArgumentParser(description="Create a new dataframe based on an original dataframe with most common mapping of barcodes")
parser.add_argument("original_dataframe", help="The original dataframe to create the new dataframe from")
parser.add_argument("new_dataframe", help="The name of the new dataframe to create")
parser.add_argument("alignment_dir", help="Directory for input/output files")
parser.add_argument("--percentage_alignment", help="The min percentage threshold for the most-common mapping", type=float, default=0.0)
parser.add_argument("--percentage_second_alignment", help="The max. percentage of the second-most common mapping", type=float, default=10.0)
parser.add_argument("--min_count", help="The minimum count threshold to filter the new dataframe", type=int, default=0)

# Parse the arguments
args = parser.parse_args()

# Read the original dataframe
df = pd.read_csv(f"{args.alignment_dir}/{args.original_dataframe}")

# # Read the original dataframe
# df = pd.read_csv("NM-TRIP_Library_ForBC_tagmentation_with_checks_output.csv")

# Create a list of unique barcodes
unique_barcodes = df["barcode"].unique()

# Create an empty list to store rows of the new dataframe
new_rows = []

# Loop through each barcode in the original dataframe
for barcode in tqdm(unique_barcodes, desc="Processing barcodes"):
    # Filter the dataframe for the current barcode
    barcode_df = df[df["barcode"] == barcode]
    
    # Create a Counter of the number of occurrences of each chromosome/for_start_position/rev_start_position combination
    counter = Counter(zip(barcode_df["chromosome"], barcode_df["for_start_position_1"], barcode_df["rev_start_position_1"]))
    total_occurrences = sum(counter.values())
    
    # Create a new row for the current barcode
    new_row = {"barcode": barcode}
    
    # Loop through the most common chromosome/for_start_position/rev_start_position combinations and add them to the new row
    for i, (key, count) in enumerate(counter.most_common()):
        if i == 10:  # only check the 10 most common mappings since rest is mostly noise
            break
        new_row[f"{i+1}_chromosome"] = key[0]
        new_row[f"{i+1}_for_start_position"] = key[1]
        new_row[f"{i+1}_rev_start_position"] = key[2]
        new_row[f"{i+1}_distance_for_rev"] = abs(key[1]-key[2])
        new_row[f"{i+1}_center_position"] = (key[1] + key[2])/2
        new_row[f"{i+1}_count"] = count
        
        # Add the percentage column for this mapping
        new_row[f"{i+1}_percentage"] = count / total_occurrences * 100
    
    # Add the new row to the list of rows for the new dataframe
    new_rows.append(new_row)

# Create the new dataframe from the list of rows
new_df = pd.DataFrame(new_rows)

# Filter the new dataframe based on the percentage threshold
new_df = new_df[new_df["1_percentage"] >= args.percentage_alignment]
new_df = new_df[new_df["1_count"] >= args.min_count]
new_df = new_df[new_df["2_percentage"] <= args.percentage_second_alignment]


# reverse_complement the barcodes so they fit the output from the "indel_PCR"
new_df['barcode'] = new_df['barcode'].apply(lambda x: str(Seq(x).reverse_complement()))

#Save the new dataframe to a CSV file
new_df.to_csv(f"{args.alignment_dir}/{args.new_dataframe}", index=False)