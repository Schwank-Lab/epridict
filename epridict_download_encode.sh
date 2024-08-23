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

# Function to show loading animation
loading_animation() {
    local pid=$1
    local delay=0.5
    local spinstr='|/-\'
    while [ "$(ps a | awk '{print $1}' | grep $pid)" ]; do
        local temp=${spinstr#?}
        printf " [%c]  " "$spinstr"
        local spinstr=$temp${spinstr%"$temp"}
        sleep $delay
        printf "\b\b\b\b\b\b"
    done
    printf "    \b\b\b\b"
}

echo "Downloading files ..."

# Use Python to handle CSV parsing and downloads
python3 - <<END &
import csv
import os
import subprocess
import concurrent.futures

def download_file(url, output, accession):
    cmd = f"curl -L -s -o '{output}' '{url}'"
    result = subprocess.run(cmd, shell=True)
    if result.returncode == 0:
        print(f"Successfully downloaded {accession}")
    else:
        print(f"Failed to download {accession}")

with open('misc/ePRIDICT_ENCODE_Selected_Datasets_Overview.csv', 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    download_tasks = []
    for row in reader:
        if row["$TARGET_COLUMN"] == "x":
            full_url = f"$BASE_URL{row['Download URL']}"
            accession = row['Accession']
            output_file = f"bigwig/{accession}.bigWig"
            download_tasks.append((full_url, output_file, accession))

    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
        executor.map(lambda x: download_file(*x), download_tasks)

END

PID=$!
loading_animation $PID

wait $PID

echo "All downloads completed."