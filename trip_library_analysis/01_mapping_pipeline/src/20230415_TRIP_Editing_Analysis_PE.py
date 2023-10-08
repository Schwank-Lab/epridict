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

target_location_within_read = 80
target_length = 18

def importfiles(input_folder, filename):
    barcodedict = {}
    with gzip.open("../"+input_folder+filename, "rt") as fasta_file:
        count = 0
        counttemp = 0
        PE_original_seq = 'TGAAGTGGTAAGAGGTCG'
        PE_edited_seq = 'TGAAGTCGTAAGAGGTCG'
        print('Check editing...')
        for seq_record in SeqIO.parse(fasta_file, 'fastq'):
            barcode = seq_record.seq[barcode_location_within_read:barcode_length]
            target = seq_record.seq[target_location_within_read:target_location_within_read+target_length]
            if target == PE_original_seq:
                edit = 'wt'
            elif target == PE_edited_seq:
                edit = 'edited'
            else:
                edit = 'unintended_editing'

            #print(edit)
            #print(target)
            #print(PE_edited_seq)
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
        barcodedict_evaluation[key] = {'wt': barcodedict[key].count('wt'), 'edited': barcodedict[key].count('edited'), 'unintended_editing': barcodedict[key].count('unintended_editing'), 'totalreads': len(barcodedict[key])}

    # create a dataframe from the dictionary
    barcodedf = pd.DataFrame.from_dict(barcodedict_evaluation,orient='index')

    # calculate the percentage of edited and unintended editing reads for each barcode
    barcodedf['edited_percentage'] = barcodedf['edited'] / barcodedf['totalreads']
    barcodedf['unintended_editing_percentage'] = barcodedf['unintended_editing'] / barcodedf['totalreads']
    barcodedf['wt_percentage'] = barcodedf['wt'] / barcodedf['totalreads']

    # sort barcodedf by barcodecount (descending)
    barcodedf = barcodedf.sort_values(by=['totalreads'], ascending=False)

    print('Editing evaluation done')
    return barcodedf



# import file overview containing the editor, replicate, input folder, and filename
filedf = pd.read_csv(overview_path+'20230415_TRIP_editing_files_overview.csv')

# only keep rows with "PE" in the column "analysis"
filedf = filedf[filedf['analysis'] == 'PE']

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