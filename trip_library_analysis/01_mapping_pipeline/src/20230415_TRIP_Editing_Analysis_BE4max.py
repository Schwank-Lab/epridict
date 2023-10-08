# -*- coding: utf-8 -*-
"""
Analysis of TRIP editing.

@author: nimath
"""
import pandas as pd
from Bio import SeqIO
import gzip
import concurrent.futures
import itertools
from collections import Counter

overview_path = '../input/'
barcode_location_within_read = 0
barcode_length = 16

target_location_within_read = 16
target_length = 20


# BE4max editing sequences within protospacer
BE4max_original_seq = 'TGATCGGTACCAACTCCAGC'
# create lists "C_positions" that contains the positions of all C's in BE4max_original_seq
C_positions = [4, 9, 10, 13, 15, 16, 19]
        
combinations = []

# Generate all combinations of length 1 to 7
for i in range(1, 8):
    comb = list(itertools.combinations(C_positions, i))
    combinations.extend(comb)

combinations = [list(x) for x in combinations]  # convert tuples to list to avoid having comma in column names for single indices
BE4max_modified_list = []
for combi in combinations:
        # Make a copy of the original string
    new_string = BE4max_original_seq[:]
    # Iterate through the list of indices in reverse order
    for index in reversed(combi):
        # Replace the character at the index with 'T'
        new_string = new_string[:index] + 'T' + new_string[index+1:]
    # Append the modified string to the list
    BE4max_modified_list.append(new_string)

distance_to_window = 0 # distance from read start to the start of quantification window



def importfiles(input_folder, filename):
    barcodedict = {}
    with gzip.open("../"+input_folder+filename, "rt") as fasta_file:
        count = 0
        counttemp = 0
        print('Check editing...')
        for seq_record in SeqIO.parse(fasta_file, 'fastq'):
            barcode = seq_record.seq[barcode_location_within_read:barcode_length]
            target = seq_record.seq[target_location_within_read:target_location_within_read+target_length]

            if target == BE4max_original_seq:
                edit = 'wt'
            else:
                edit = 'modified'
                # print(target)
                # print(BE4max_original_seq)
                # print()
                if target in BE4max_modified_list:  # if target is a (or multiple) C to T editing, then label it as such instead of only modified:
                    ind = BE4max_modified_list.index(target)
                    edit = str(combinations[ind])  # needs to be converted to a string since script will break otherwise downstream
            #print(target)
            #print(BE4max_original_seq)
            #print(edit)
            #print()
            if barcode in barcodedict:
                    barcodedict[barcode].append(edit)
            else:
                barcodedict[barcode] = [edit]

            counttemp+=1
            count+=1
            if counttemp == 100000:
                print(count)
                counttemp = 0
                
            
        
    editingpercentagedict = {}
    for barcode in barcodedict:
        editingpercentagedict[barcode] = Counter(barcodedict[barcode])

    for key, counter in editingpercentagedict.items():
        total = sum(counter.values())  # Total number of counts for this key
        for element, count in counter.items():
            counter[element] = count / total

    editingpercentagedf = pd.DataFrame.from_dict(editingpercentagedict,orient='index')
    editingpercentagedf = editingpercentagedf.fillna(0)
    for index, row in editingpercentagedf.iterrows():
        barcode = index
        editingpercentagedf.at[index, 'barcodecount'] = len(barcodedict[barcode])

    # sort editingpercentagedf by barcodecount (descending)
    editingpercentagedf = editingpercentagedf.sort_values(by=['barcodecount'], ascending=False)
    print('Editing evaluation done')
    return editingpercentagedf

# import file overview containing the editor, replicate, input folder, and filename
filedf = pd.read_csv(overview_path+'20230415_TRIP_editing_files_overview.csv')

# only keep rows with "PE" in the column "analysis"
filedf = filedf[filedf['analysis'] == 'BE4max']

def process_row(row):
    shortname = row.shortname
    editor = row.editor
    replicate = str(row.replicate)
    analysis = row.analysis
    input_folder = row.input_folder
    filename = row.filename

    barcodedf = importfiles(input_folder, filename)
    #print('Length of barcodedf: ', len(barcodedf))
    #print('test1')
    savepath = '../results/editing/analysis/'
    savename = savepath+'20231415_'+shortname+'_'+editor+'_'+replicate+'_'+analysis+'_editing_per_barcode.csv'
    #print('test2')
    #print(savename)
    barcodedf.to_csv(savename)

# Create a ThreadPoolExecutor to run the loop in parallel
with concurrent.futures.ThreadPoolExecutor() as executor:
    # Use the executor to run process_row() for each row in the DataFrame
    executor.map(process_row, filedf.itertuples(index=False))

   