# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 16:20:38 2023
Map locations to editing dataframe

@author: nicol
"""

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Process editing and mapping files')
parser.add_argument('--editingpath', type=str, help='path to the editing file')
parser.add_argument('--editingfile', type=str, help='editing file name')
parser.add_argument('--finalname', type=str, help='final editing file name')
parser.add_argument('--mappingpath', type=str, help='path to the mapping file')
parser.add_argument('--mappingfile', type=str, help='mapping file name')
args = parser.parse_args()

editingdf = pd.read_csv(args.editingpath + args.editingfile)
editingdf = editingdf.set_index('barcode')

mappingdf = pd.read_csv(args.mappingpath + args.mappingfile)
mappingdf = mappingdf.set_index('barcode')

# Filter out barcodes where distance is not 4bp between for and rev location
mappingdf = mappingdf[mappingdf["1_distance_for_rev"] == 4]

for ind, row in mappingdf.iterrows():
    chromosome = row["1_chromosome"]
    center_position = row["1_center_position"]
    if ind in editingdf.index:
        editingdf.loc[ind,'chromosome'] = chromosome
        editingdf.loc[ind,'center_position'] = center_position

final_editingdf = editingdf.dropna(subset=['chromosome'])

# Get editing columns:
cols = [col for col in final_editingdf.columns if col.endswith("_corrected_percentage")]

cols.extend(['chromosome', 'center_position'])
final_editingdf = final_editingdf[cols]

final_editingdf.to_csv(args.editingpath + args.finalname)

