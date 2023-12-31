#!/bin/bash
# Shell script to download ENCODE bigWig files for ePRIDICT


# ENCODE base URL
BASE_URL="https://www.encodeproject.org"

# Get the model type from the command line argument; default is 'light'
MODEL_TYPE=${1:-light}

# Confirm the download based on the model type
if [ "$MODEL_TYPE" == "full" ]; then
  read -p "This will download 455 ENCODE datasets with a total of 624 GB to the folder 'bigwig'. Are you sure to continue (yes/no)? " CONFIRM
  TARGET_COLUMN="full_model"
elif [ "$MODEL_TYPE" == "light" ]; then
  read -p "This will download 6 ENCODE datasets with a total of 5.3 GB to the folder 'bigwig'. Are you sure to continue (yes/no)? " CONFIRM
  TARGET_COLUMN="light_model"
else
  echo "Invalid model type. Please specify either 'light' or 'full'."
  exit 1
fi

# Exit if the user does not confirm
case "$CONFIRM" in
  [yY] | [yY][eE][sS])
    # proceed with the operation
    ;;
  [nN] | [nN][oO])
    echo "Operation cancelled."
    exit 0
    ;;
  *)
    echo "Invalid answer. Operation cancelled."
    exit 1
    ;;
esac

# Use Python for CSV parsing and row filtering
python3 - <<END
import csv
import os

with open('misc/ePRIDICT_ENCODE_Selected_Datasets_Overview.csv', 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    count = 0
    commands = []
    for row in reader:
        if row["$TARGET_COLUMN"] == "x":
            full_url = f"$BASE_URL{row['Download URL']}"
            accession = row['Accession']
            commands.append(f"wget -O bigwig/{accession}.bigWig {full_url}")
            count += 1
            if count % 6 == 0:
                os.system(" & ".join(commands) + " & wait")
                commands = []
    if commands:
        os.system(" & ".join(commands) + " & wait")
END