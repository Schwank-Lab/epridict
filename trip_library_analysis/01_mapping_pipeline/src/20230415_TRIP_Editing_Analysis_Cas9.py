# -*- coding: utf-8 -*-
"""
Analysis of TRIP editing.

@author: nimath
"""
import pandas as pd
from Bio import SeqIO
import gzip
import concurrent.futures

overview_path = '../input/'
barcode_location_within_read = 0
barcode_length = 16

# use 4 bp around nick as quantification window for Cas9 editing
target_location_within_read = 29
target_length = 8

def importfiles(input_folder, filename):
    barcodedict = {}
    with gzip.open("../"+input_folder+filename, "rt") as fasta_file:
        count = 0
        counttemp = 0
        Cas9_original_seq = 'CTCCAGCA'
        print('Check editing...')
        for seq_record in SeqIO.parse(fasta_file, 'fastq'):
            barcode = seq_record.seq[barcode_location_within_read:barcode_length]
            target = seq_record.seq[target_location_within_read:target_location_within_read+target_length]
            if target == Cas9_original_seq:
                edit = 'wt'
            else:
                edit = 'modified'

            #print(edit)
            #print(target)
            #print(Cas9_original_seq)
            #print(barcode)
            #print()
            
            # add edit to the list of edits for this barcode (barcodedict) and create a new list if it does not exist yet
            if barcode in barcodedict:
                barcodedict[barcode].append(edit)
            else:
                barcodedict[barcode] = [edit]
            counttemp+=1
            count+=1
            if counttemp == 100000:
                print(count)
                counttemp = 0
    # for each key in barcodedict, count the number of 'wt', 'edited' and 'unintended_editing' and add this to new dictionary "barcodedict_evaluation"
    # also add a "totalreads" key to the dictionary which contains the total number of reads for this barcode
    barcodedict_evaluation = {}
    for key in barcodedict:
        barcodedict_evaluation[key] = {'wt': barcodedict[key].count('wt'), 'modified': barcodedict[key].count('modified'), 'totalreads': len(barcodedict[key])}

    # create a dataframe from the dictionary
    barcodedf = pd.DataFrame.from_dict(barcodedict_evaluation,orient='index')

    # calculate the Cas9rcentage of edited and unintended editing reads for each barcode
    barcodedf['modified_percentage'] = barcodedf['modified'] / barcodedf['totalreads']
    barcodedf['wt_percentage'] = barcodedf['wt'] / barcodedf['totalreads']

    # sort barcodedf by barcodecount (descending)
    barcodedf = barcodedf.sort_values(by=['totalreads'], ascending=False)

    print('Editing evaluation done')
    return barcodedf



# import file overview containing the editor, replicate, input folder, and filename
filedf = pd.read_csv(overview_path+'20230415_TRIP_editing_files_overview.csv')

# only keep rows with "Cas9" in the column "analysis"
filedf = filedf[filedf['analysis'] == 'Cas9']

def process_row(row):
    shortname = row.shortname
    editor = row.editor
    replicate = str(row.replicate)
    analysis = row.analysis
    input_folder = row.input_folder
    filename = row.filename

    barcodedf = importfiles(input_folder, filename)

    savepath = '../results/editing/analysis/'
    barcodedf.to_csv(savepath+'20231415_'+shortname+'_'+editor+'_'+replicate+'_'+analysis+'_editing_per_barcode.csv')

# Create a ThreadPoolExecutor to run the loop in parallel
with concurrent.futures.ThreadPoolExecutor() as executor:
    # Use the executor to run process_row() for each row in the DataFrame
    executor.map(process_row, filedf.itertuples(index=False))