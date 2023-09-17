#!/bin/bash
# Shell script to download ENCODE bigWig files for ePRIDICT


# ENCODE base URL
BASE_URL="https://www.encodeproject.org"

# Get the model type from the command line argument; default is 'slim'
MODEL_TYPE=${1:-slim}

# Confirm the download based on the model type
if [ "$MODEL_TYPE" == "full" ]; then
  read -p "This will download 455 ENCODE datasets with a total of 624 GB to the folder 'bigwig'. Are you sure to continue (yes/no)? " CONFIRM
  TARGET_COLUMN="full_model"
elif [ "$MODEL_TYPE" == "slim" ]; then
  read -p "This will download 6 ENCODE datasets with a total of 5.3 GB to the folder 'bigwig'. Are you sure to continue (yes/no)? " CONFIRM
  TARGET_COLUMN="slim_model"
else
  echo "Invalid model type. Please specify either 'slim' or 'full'."
  exit 1
fi

# Exit if the user does not confirm
if [ "$CONFIRM" != "yes" ]; then
  echo "Operation cancelled."
  exit 0
fi

# Use Python for CSV parsing and row filtering
python3 - <<END
import csv
import os

with open('misc/ePRIDICT_ENCODE_Selected_Datasets_Overview.csv', 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        if row["$TARGET_COLUMN"] == "x":
            full_url = f"$BASE_URL{row['Download URL']}"
            accession = row['Accession']
            os.system(f"wget -O bigwig/{accession}.bigWig {full_url}")
END