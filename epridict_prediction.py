### epigenetic PRIme editing preDICTion (ePRIDICT) from genomic location (hg38) trained in K562 ###

import os
import pandas as pd
from joblib import Parallel, delayed
import numpy as np
import xgboost as xgb
import pyBigWig
import argparse
import bisect
from tqdm import tqdm

### Settings ###

model = "slim"  # "slim" or "full"

bigwig_folder = "./bigwig"  # Folder containing bigWig files; files must be downloaded first (see README.md)

input_folder = "./input"  # Folder containing input files for batch analysis; formatting see README.md
output_folder = "./predictions"  # Folder to save output files

### End of settings ###

genomewide_prediction_slim_df = pd.read_csv('misc/genome_wide_ePRIDICT_slim_predictions_K562.csv')

def find_closest_percentile(df: pd.DataFrame, value_to_find: float) -> float:
    """
    Find the closest percentile to a given value in a sorted dataframe.
    
    Parameters:
    - df (pd.DataFrame): The dataframe containing genomic data. Assumes it is sorted by 'prediction'.
    - value_to_find (float): The value for which to find the closest percentile.
    
    Returns:
    - float: The percentile of the closest value.
    """
    # Perform binary search to find the index of the closest value
    idx = bisect.bisect_left(df['prediction'].values, value_to_find)
    
    # Determine the closest value based on the index
    if idx == 0:
        closest_value = df.iloc[0]
    elif idx == len(df):
        closest_value = df.iloc[-1]
    else:
        before = df.iloc[idx - 1]
        after = df.iloc[idx]
        if after['prediction'] - value_to_find < value_to_find - before['prediction']:
            closest_value = after
        else:
            closest_value = before
    
    # Retrieve the percentile of the closest value
    closest_percentile = closest_value['percentile']
    
    return closest_percentile

def create_directory(folder_name, directory="current"):
    """create directory/folder (if it does not exist) and returns the path of the directory
       Args:
           folder_name: string representing the name of the folder to be created
       Keyword Arguments:
           directory: string representing the directory where to create the folder
                      if `current` then the folder will be created in the current directory
    """
    if directory == "current":
        path_current_dir = os.path.dirname(__file__)  # __file__ refers to utilities.py
    else:
        path_current_dir = directory
    path_new_dir = os.path.join(path_current_dir, folder_name)
    if not os.path.exists(path_new_dir):
        os.makedirs(path_new_dir)
    return(path_new_dir)

# Function to parse encode name and bin size from column name
def parse_column_name(column_name):
    parts = column_name.split("_")
    encode_name = parts[1] + ".bigWig"
    bin_size = int(parts[-1])
    return encode_name, bin_size

# Define function to process each bigWig file
def process_bigwig_file(bigwig_file, chromosome, position, bin_size):
    # Open the bigWig file
    bw = pyBigWig.open(os.path.join("./bigwig", bigwig_file))  # Look for the file in the "bigwig" subfolder

    # Calculate start and end positions
    start = position - bin_size
    end = position + bin_size

    # Get values from bigWig file for the genomic location
    try:
        # Attempt to get values from the BigWig file
        values = bw.values(chromosome, start, end)
    except RuntimeError:
        # If an error occurs, fetch the valid range for the chromosome
        chrom_info = bw.chroms()
        if chromosome in chrom_info:
            valid_end = chrom_info[chromosome]
            valid_start = 1  # Generally, chromosome positions start from 1
            raise RuntimeError(f"Position error! {position} is not a valid position for the prediction. Position must be within {valid_start+5000} to {valid_end-5000} for {chromosome} (hg38).")
    
    # Compute average value
    average_value = sum(values) / len(values)

    bw.close()

    return average_value

# Function to predict efficiency based on genomic location and model file
def predict_efficiency_single(chromosome, position, output_dir, model_file, model_name, preloaded_model=None):

    # Load xgboost model only if not preloaded
    if preloaded_model is None:
        model = xgb.Booster()
        model.load_model(model_file)
        single = True
    else:
        model = preloaded_model
        single = False

    # Convert position to int
    position = int(position)

    chromosomelist = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7","chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14","chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21","chr22", "chrX", "chrY"]
    if not chromosome in chromosomelist:
        raise ValueError(f"{chromosome} is not a valid chromosome argument. Chromosome must be one of the following: 'chr1', 'chr2', ... ,'chrX', 'chrY'")
    # Define a helper function for parallel processing
    def process_column(column_name):
        # Parse encode name and bin size
        encode_name, bin_size = parse_column_name(column_name)

        # Compute average value
        average_value = process_bigwig_file(encode_name, chromosome, position, bin_size)

        return average_value

    # Use joblib to apply the helper function in parallel
    results = Parallel(n_jobs=-1)(delayed(process_column)(col) for col in column_names_dict[model_name])

    # Create a dataframe from the results
    features_df = pd.DataFrame(dict(zip(column_names_dict[model_name], results)), index=[0])

    # Check if any feature is NaN or all features less than 0.5
    if features_df.isnull().values.any() or np.all(features_df.values < 0.5):
        features_df = features_df.transpose()
        # rename column "0" to "value"
        features_df = features_df.rename(columns={0: "feature_value"})
        if single:
            features_df.to_csv(os.path.join(output_dir, f"{chromosome}_{position}_failed.csv"))
            print("One or more features are missing or less than 0.5 at this genomic location. Prediction is therefore not possible.")
        return None, None, features_df, "One or more features are missing or less than 0.5 at this genomic location. Prediction is therefore not possible."

    # Use the model to predict efficiency
    # Convert the DataFrame into DMatrix
    data_dmatrix = xgb.DMatrix(data=features_df)

    # Use the model to predict efficiency
    prediction = model.predict(data_dmatrix)

    percentile = find_closest_percentile(genomewide_prediction_slim_df, prediction[0])

    features_df[f"ePRIDICT_prediction_{model_name}"] = prediction[0]
    features_df[f"percentile_slim_ePRIDICT_genomewide_K562"] = percentile
    features_df = features_df.transpose()
    # rename column "0" to "value"
    features_df = features_df.rename(columns={0: "feature_value"})

    if single:
        # Save the prediction to a csv file
        features_df.to_csv(os.path.join(output_dir, f"{chromosome}_{position}.csv"))

        print(f"ePRIDICT score ({model_name} model): ", round(prediction[0], 2))
        print(f"Percentile in context of {len(genomewide_prediction_slim_df)} sampled genomic locations in K562: {round(percentile,2)}%")
        print()
        print(f"Output stored in {output_dir} as {chromosome}_{position}.csv")
        print()
        print()

    return prediction[0], percentile, features_df, ""

def predict_efficiency_df(batchdf, out_dir, model_file, model_name):
    print("...")

    feature_dfs = []

    # Check if required columns exist in the DataFrame
    if "chromosome" not in batchdf.columns or "position_hg38" not in batchdf.columns:
        raise ValueError("The DataFrame must contain 'chromosome' and 'position_hg38' columns.")

    # Initialize lists to store predictions and percentiles
    predictions = []
    percentiles = []
    notelist = []

    # Load xgboost model once to avoid reloading for each prediction
    model = xgb.Booster()
    model.load_model(model_file)

    for idx, row in tqdm(batchdf.iterrows(), total=len(batchdf), desc="Batch prediction progress"):
        chromosome = row['chromosome']
        position = row['position_hg38']
        notes = ""

        try:
            prediction, percentile, feature_df_row, notes = predict_efficiency_single(chromosome, position, out_dir, model_file, model_name, preloaded_model=model)
            feature_df_row = feature_df_row.transpose()
            feature_df_row = feature_df_row.drop(columns=[col for col in ['ePRIDICT_prediction_slim', 'percentile_slim_ePRIDICT_genomewide_K562'] if col in feature_df_row.columns])

        except Exception as e:
            prediction = None
            percentile = None
            feature_df_row = pd.DataFrame(data=np.nan, index=[0], columns=column_names_dict[model_name])
            notes = str(e)

        # Append to lists
        predictions.append(prediction)
        percentiles.append(percentile)
        notelist.append(notes)
        feature_dfs.append(feature_df_row)

    # Add the lists as new columns to the DataFrame
    batchdf[f"ePRIDICT score ({model_name} model)"] = predictions
    batchdf[f"ePRIDICT percentile ({model_name} model)"] = percentiles
    batchdf[f"Error notes"] = notelist

    # Concatenate the feature DataFrames and add to batchdf
    feature_df_combined = pd.concat(feature_dfs, ignore_index=True)
    batchdf = pd.concat([batchdf, feature_df_combined], axis=1)

    batchdf.to_csv(os.path.join(out_dir, f"{out_fname}.csv"), index=False)
    print("Done!")
    return batchdf


# Read the full model column names from the text file
def read_model_columns(file_path):
    with open(file_path, 'r') as f:
        return [line.strip() for line in f.readlines()]

# Hardcoded slim model columns
# slim_model_columns = ['Dnase-seq_ENCFF972GVB_average_value_100',
#  'Dnase-seq_ENCFF972GVB_average_value_1000',
#  'Dnase-seq_ENCFF972GVB_average_value_2000',
#  'Dnase-seq_ENCFF972GVB_average_value_5000',
#  'HDAC2_ENCFF954LGE_average_value_100',
#  'HDAC2_ENCFF954LGE_average_value_1000',
#  'HDAC2_ENCFF954LGE_average_value_2000',
#  'HDAC2_ENCFF954LGE_average_value_5000',
#  'H3K4me1_ENCFF834SEY_average_value_100',
#  'H3K4me1_ENCFF834SEY_average_value_1000',
#  'H3K4me1_ENCFF834SEY_average_value_2000',
#  'H3K4me1_ENCFF834SEY_average_value_5000',
#  'H3K9me3_ENCFF601JGK_average_value_100',
#  'H3K9me3_ENCFF601JGK_average_value_1000',
#  'H3K9me3_ENCFF601JGK_average_value_2000',
#  'H3K9me3_ENCFF601JGK_average_value_5000',
#  'H3K27me3_ENCFF139KZL_average_value_100',
#  'H3K27me3_ENCFF139KZL_average_value_1000',
#  'H3K27me3_ENCFF139KZL_average_value_2000',
#  'H3K27me3_ENCFF139KZL_average_value_5000',
#  'H3K4me2_ENCFF959YJV_average_value_100',
#  'H3K4me2_ENCFF959YJV_average_value_1000',
#  'H3K4me2_ENCFF959YJV_average_value_2000',
#  'H3K4me2_ENCFF959YJV_average_value_5000']

# Read full model columns from the text file
full_model_columns = read_model_columns('misc/ePRIDICT_full_model_column_names.txt')

slim_model_columns = read_model_columns('misc/ePRIDICT_slim_model_column_names.txt')

# Update the column_names_dict with the new full model columns
column_names_dict = {'slim': slim_model_columns, 'full': full_model_columns}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Running ePRIDICT to predict prime editing efficiencies in K562 based on chromatin context of a genomic location.")

    subparser = parser.add_subparsers(dest='command')
    manual_m = subparser.add_parser('manual')
    batch_m  = subparser.add_parser('batch')

    manual_m.add_argument("--chromosome", type=str, help="Chromosome of location to predict. ('chr1', 'chr2', ... , 'chrX'", required=True)
    manual_m.add_argument("--position_hg38", type=str, help="Genomic position of location to predict within selected chromosome. ('11939204',...)", required=True)
    manual_m.add_argument("--use_full_model", action='store_true', help="Use full model with 455 ENCODE datasets. Default is to use slim model")

    batch_m.add_argument("input_fname", type=str, help="Input filename - name of csv file that has three columns {identifier, chromosome ('chr1', 'chr2', ...), position (within chromosome)}. See batch_template.csv in the ./input folder ")
    # batch_m.add_argument("--input-dir", type=str, default=input_folder, help="Input directory where the input csv file is found on disk")
    # batch_m.add_argument("--output-dir", type=str, default=output_folder, help="Output directory where results are dumped on disk")
    batch_m.add_argument("--output-fname", type=str, help="Output filename for the resulting dataframe. If not specified, the name of the input file will be used")
    batch_m.add_argument("--use_full_model", action='store_true', help="Use full model with 455 ENCODE datasets. Default is to use slim model")
  
    args = parser.parse_args()

    if args.command == 'manual':
        print('Running in manual mode...')

        out_dir = create_directory(output_folder, os.getcwd())
        
        try:
            chromosome = args.chromosome.lower()
            position = int(args.position_hg38)
        except:
            print("Check your input. Chromosome must be a string (e.g. 'chr1') and position must be an integer (e.g. 11939204)")
            exit()
                   
        if args.use_full_model:
            model="full"
        else:
            model="slim"

        if model == "slim":
            model_file = "models/epridict_slim_xgboost_model.json"
        elif model == "full":
            model_file = "models/epridict_full_xgboost_model.json"
        
        predict_efficiency_single(chromosome, position, out_dir, model_file, model)

    elif args.command == 'batch':
        print('Running in batch mode...')

        inp_dir = create_directory(input_folder, os.getcwd())
        print('input directory:', inp_dir)

        inp_fname = args.input_fname

        out_dir = create_directory(output_folder, os.getcwd())
        print('output directory:', out_dir)

        if args.output_fname:
            out_fname = args.output_fname
        else:
            out_fname = args.input_fname.split('.')[0]+"_output"
        
        if args.use_full_model:
            model="full"
        else:
            model="slim"

        if model == "slim":
            model_file = "models/epridict_slim_xgboost_model.json"
        elif model == "full":
            model_file = "models/epridict_full_xgboost_model.json"

        batchdf = pd.read_csv(os.path.join(inp_dir, inp_fname))

        predict_efficiency_df(batchdf, out_dir, model_file, model)

    else:
        print('Please specify how to run ePRIDICT ("manual" or "batch") as argument after the script name.')