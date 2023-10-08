# -*- coding: utf-8 -*-
"""
Analysis of self-targeting library.

@author: nimath
"""
import pandas as pd
from Bio import SeqIO
import gzip
from os import listdir
import numpy as np
import os
import time
from Bio.Seq import Seq
import argparse
from itertools import combinations

# Get the current working directory
cwd = os.path.join(os.getcwd(), '')


def lookup(prototemplate, exttemplate, barcodetemplate):  # generate lookup dictionaries for all three relevant regions (protospacer, RTstart, endoftarget)
    protolookup = {}
    for i, z in enumerate(prototemplate):
            protolookup.setdefault(z, []).append(i)
    
    rtlookup = {}
    for i, z in enumerate(exttemplate):
            rtlookup.setdefault(z[:6], []).append(i)  # only consider first 15 bases of RT
    
    barcodelookup = {}
    for i, z in enumerate(barcodetemplate):
            barcodelookup.setdefault(z, []).append(i)
      
    
    return protolookup, rtlookup, barcodelookup


def importfiles(filename, protolookup, rtlookup, barcodelookup):
    diseasedict = {}
    filename = shortname+'_Protospacer'+filtered+'.fastq.gz'
    with gzip.open(path+filename, "rt") as fasta_file:
        count = 0
        counttemp = 0
        print('Start Protospacer lookup loop...')
        for seq_record in SeqIO.parse(fasta_file, 'fastq'):
            protospacer = seq_record.seq
            protomatch = protolookup.get("G"+str(protospacer))
            identifier = seq_record.id
            if identifier in diseasedict:
                diseasedict[identifier]["protomatch"] = protomatch
            else:
                diseasedict[identifier] = {'protomatch':protomatch}
            count+=1
            counttemp+=1
            if counttemp == 1000000:
                print(count)
                counttemp = 0
    print('Protospacer done')
    
    
    filename = shortname+'_target_5trim'+filtered+'.fastq.gz'
    with gzip.open(path+filename, "rt") as fasta_file:
        count = 0
        counttemp = 0
        print('Start barcode lookup loop...')
        for seq_record in SeqIO.parse(fasta_file, 'fastq'):
            barcode = seq_record.seq[5:13]
            # print(barcode)
            barcodematch = barcodelookup.get(barcode)
            # print(barcodematch)
            # print()
            identifier = seq_record.id
            if identifier in diseasedict:
                diseasedict[identifier]["barcodematch"] = barcodematch
            # no else since we are only interested in reads which have protospacer match
            # else:
            #     diseasedict[identifier] = {'barcode':str(barcode)}
            count+=1
            counttemp+=1
            if counttemp == 1000000:
                print(count)
                counttemp = 0
    print('Barcode done')
      
    
    filename = shortname+'_RTTtrim'+filtered+'.fastq.gz'
    with gzip.open(path+filename, "rt") as fasta_file:
        count = 0
        counttemp = 0
        print('Start extension lookup loop...')
        for seq_record in SeqIO.parse(fasta_file, 'fastq'):
            rtt = seq_record.seq
            exttemplate = str(rtt[:6])
            exttemplatematch = rtlookup.get(exttemplate)

            # print(exttemplate)
            # print(exttemplatematch)
            # print()
            identifier = seq_record.id
            if identifier in diseasedict:
                diseasedict[identifier]["exttemplatematch"] = exttemplatematch

            count+=1
            counttemp+=1
            if counttemp == 1000000:
                print(count)
                counttemp = 0

    print('extension done')
    
    return diseasedict


def mergeDict(dict1, dict2):
        ''' Merge dictionaries and keep values of common keys in list'''
        dict3 = {**dict1, **dict2}
        for key, value in dict3.items():
            if key in dict1 and key in dict2:
                dict3[key] = [value , dict1[key]]
        return dict3


def get_intermediate_strings(str1, str2):
    '''Get all intermediate strings between two strings with one or more differences. Used to find all possible intermediate strings between a mutated and a wildtype sequence, to quantify incomplete editing.'''
    intermediate_strings = []
    changed_positions = []
    diff_indexes = [i for i in range(len(str1)) if str1[i] != str2[i]]
    
    for i in range(len(diff_indexes)):
        for comb in combinations(diff_indexes, i + 1):
            intermediate = list(str1)
            changes = []
            for index in comb:
                intermediate[index] = str2[index]
                changes.append(index)
            intermediate = "".join(intermediate)
            if intermediate != str1 and intermediate != str2:
                intermediate_strings.append(intermediate)
                changed_positions.append(changes)
                
    return intermediate_strings, changed_positions


scaffold = 'GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC'
templatedforiginal = pd.read_csv('20220922_30k_library_withbystander_with_barcodes_and_primers.csv')
templatedforiginal['barcode_rev'] = templatedforiginal['barcode'].apply(lambda x: str(Seq(x).reverse_complement()))
templatedforiginal['extension'] = templatedforiginal['RTrevcomp'] + templatedforiginal['PBSrevcomp13bp']
templatedforiginal['extension_15bp'] = templatedforiginal['extension'].apply(lambda x: x[:15])


# apply the get_intermediate_strings function to the rows of the dataframe which have "MultiBpReplacement" in "Mutation_Type" and add the resulting first part of the tuple to the column "wide_mutated_target_intermediate" and the second to "deepeditposition_intermediate"

templatedforiginal['wide_mutated_target_intermediate'] = templatedforiginal.apply(lambda x: get_intermediate_strings(x['wide_initial_target'], x['wide_mutated_target'])[0] if x['Mutation_Type'] == 'MultibpReplacement' else None, axis=1)
templatedforiginal['deepeditposition_intermediate'] = templatedforiginal.apply(lambda x: ['_'.join(str(item) for item in inner_list) for inner_list in get_intermediate_strings(x['wide_initial_target'], x['wide_mutated_target'])[1]] if x['Mutation_Type'] == 'MultibpReplacement' else None, axis=1)

templatedf = templatedforiginal.copy()

path = './Output/'

# get filename from argument
parser = argparse.ArgumentParser(description='Process a given file.')
parser.add_argument('filename', type=str, help='The input file to process.')
args = parser.parse_args()
filelist = [args.filename]


filtered = ''

for filename in filelist:  # loop through all SAMPLES in a directory with "R1" in the name (listfiles1 function)
    print(filename)
    shortname = filename[0:-21]
    prototemplate = templatedf['Protospacer-Sequence'].values.flatten()
    exttemplate = templatedf['extension_15bp'].values.flatten()
    barcodetemplate = templatedf['barcode_rev'].values.flatten()
    

    protolookup, rtlookup, barcodelookup = lookup(prototemplate,exttemplate, barcodetemplate)
    del prototemplate, exttemplate, barcodetemplate
    diseasedict = importfiles(filename,protolookup, rtlookup, barcodelookup)
    diseasedf = pd.DataFrame.from_dict(diseasedict,orient='index')  # make dataframe from dict
    del diseasedict

    diseasedf.dropna(subset = ['protomatch', 'exttemplatematch', 'barcodematch'], inplace=True) # do not drop elements without barcode match, since barcode match only includes barcodes from ambiguous sequences
    diseasedf['match'] = [list(set(a).intersection(set(b))) for a, b in zip(diseasedf.protomatch, diseasedf.barcodematch)]
    diseasedf['matchwithextension'] = [list(set(a).intersection(set(b), set(c))) for a, b, c in zip(diseasedf.protomatch, diseasedf.exttemplatematch, diseasedf.barcodematch)]
    diseasedf['recombination'] = diseasedf.apply(lambda x: True if len(x.matchwithextension) == 0 else False,axis=1)

    final_diseasedf = diseasedf[diseasedf['match'].map(lambda d: len(d)) > 0][['match','matchwithextension']]  # all reads which have at least one match; those with multiple matches will be assigned by barcodematch

print("calculating...")
templatedf["editedcount"] = 0
templatedf["uneditedcount"] = 0
templatedf["indelcount"] = 0
templatedf["nickindelcount"] = 0
templatedf["intermediatecount"] = 0
templatedf["intermediatepositions"] = np.empty((len(templatedf), 0)).tolist()

templatedf["beforeflapindelcount"] = 0
templatedf["totalreads"] = 0
templatedf['barcodes'] = np.empty((len(templatedf), 0)).tolist()
filename = shortname+'_target_5trim'+filtered+'.fastq.gz'
templatenumpy = templatedf.to_numpy()

# c = 0
readcounter = 0
readcounttemp = 0
start = time.time()

# define column position for used column in numpy array:
uneditedcountnr = templatedf.columns.get_loc("uneditedcount")
editedcountnr = templatedf.columns.get_loc("editedcount")
indelcountnr = templatedf.columns.get_loc("indelcount")
nickindelcountnr = templatedf.columns.get_loc("nickindelcount")
intermediatecountnr = templatedf.columns.get_loc("intermediatecount")
intermediatepositionsnr = templatedf.columns.get_loc("intermediatepositions")
beforeflapindelcountnr = templatedf.columns.get_loc("beforeflapindelcount")
wide_initial_targetnr = templatedf.columns.get_loc("wide_initial_target")
wide_mutated_targetnr = templatedf.columns.get_loc("wide_mutated_target")
wide_mutated_target_intermediate = templatedf.columns.get_loc("wide_mutated_target_intermediate")
deepeditposition_intermediatenr = templatedf.columns.get_loc("deepeditposition_intermediate")
editingpositionnr = templatedf.columns.get_loc("Editing_Position")
RTTlengthpositionnr = templatedf.columns.get_loc("RTlength")
correction_typepositionnr = templatedf.columns.get_loc("Correction_Type")
mutation_type_positionnr = templatedf.columns.get_loc("Mutation_Type")
correction_lengthpositionnr = templatedf.columns.get_loc("Correction_Length")
namenr = templatedf.columns.get_loc("Name")

basesbeforenick = 2
baseafterflap = 5
distanceto_wide_start = 8 # number of bases in target read that come before wide_initial_target, 
                          #  note that the first 5 bases are not matching wide_target, but they are not used for analysis here.

with gzip.open(path+filename, "rt") as fasta_file:
    for seq_record in SeqIO.parse(fasta_file, 'fastq'):
        # only analyze high quality reads
        if len(seq_record.letter_annotations["phred_quality"]) > 0:
            if sum(seq_record.letter_annotations["phred_quality"])/len(seq_record.letter_annotations["phred_quality"]) > 27.5:
                targetend = str(seq_record.seq)
                
                identifier = seq_record.id
                readcounter+=1
                readcounttemp+=1
                
                if readcounttemp == 100000:
                    print(readcounter)
                    readcounttemp = 0
                if not identifier in final_diseasedf.index:
                    # if identifier in multiple_diseasedf.index:
                        
                    continue
                variantindex = diseasedf.at[identifier,'match'][0]
                
                RTTlength = int(templatenumpy[variantindex,RTTlengthpositionnr])
                correction_type = templatenumpy[variantindex,correction_typepositionnr]
                mutation_type = templatenumpy[variantindex,mutation_type_positionnr]  # can distinguish between 1bp and multibp replacements
                correction_length = int(templatenumpy[variantindex,correction_lengthpositionnr])
                
                if correction_type == 'Replacement':
                    sequenceWT = targetend[distanceto_wide_start+21-basesbeforenick:distanceto_wide_start+21+RTTlength+baseafterflap]
                    sequenceWT_beforeflap = targetend[distanceto_wide_start+21-basesbeforenick:distanceto_wide_start+21+RTTlength-2]  # full sequence until 2bp before flap
                    sequenceWT_nick = targetend[distanceto_wide_start+21-basesbeforenick:distanceto_wide_start+21+2]
                    
                    sequenceMUT = targetend[distanceto_wide_start+21-basesbeforenick:distanceto_wide_start+21+RTTlength+baseafterflap]
                    sequenceMUT_beforeflap = targetend[distanceto_wide_start+21-basesbeforenick:distanceto_wide_start+21+RTTlength-2]  # full sequence until 2bp before flap
                    sequenceMUT_nick = targetend[distanceto_wide_start+21-basesbeforenick:distanceto_wide_start+21+2]
                    
                    controlWT = templatenumpy[variantindex,wide_initial_targetnr][21-basesbeforenick:21+RTTlength+baseafterflap]
                    controlWT_beforeflap = templatenumpy[variantindex,wide_initial_targetnr][21-basesbeforenick:21+RTTlength-2]
                    controlWT_nick = templatenumpy[variantindex,wide_initial_targetnr][21-basesbeforenick:21+2]
                    
                    controlMUT = templatenumpy[variantindex,wide_mutated_targetnr][21-basesbeforenick:21+RTTlength+baseafterflap]
                    controlMUT_beforeflap = templatenumpy[variantindex,wide_mutated_targetnr][21-basesbeforenick:21+RTTlength-2]
                    controlMUT_nick = templatenumpy[variantindex,wide_mutated_targetnr][21-basesbeforenick:21+2]
                    
        
                    if sequenceWT == controlWT:
                        templatenumpy[variantindex,uneditedcountnr] += 1 # uneditedcount is column 42
                    elif sequenceMUT == controlMUT:
                        templatenumpy[variantindex,editedcountnr] += 1
                    else:
                        # print(sequenceWT)
                        # print(controlWT)
                        templatenumpy[variantindex,indelcountnr] += 1

                        if mutation_type == 'MultibpReplacement':
                            intermediate_editlist = templatenumpy[variantindex,wide_mutated_target_intermediate]
                            deepeditposition_intermediatelist = templatenumpy[variantindex,deepeditposition_intermediatenr]
                            # check if sequenceWT or sequenceMUT occur in intermediate_editlist and store index in variable
                            if sequenceWT in intermediate_editlist:
                                intermediateindex = intermediate_editlist.index(sequenceWT)
                                templatenumpy[variantindex,intermediatecountnr] += 1
                                templatenumpy[variantindex,intermediatepositionsnr].append(deepeditposition_intermediatelist[intermediateindex])
                            elif sequenceMUT in intermediate_editlist:
                                intermediateindex = intermediate_editlist.index(sequenceMUT)
                                templatenumpy[variantindex,intermediatecountnr] += 1
                                templatenumpy[variantindex,intermediatepositionsnr].append(deepeditposition_intermediatelist[intermediateindex])

                        
                        if (sequenceWT_nick != controlWT_nick) and (sequenceMUT_nick != controlMUT_nick):  # check if 4bp window around nick has unintended edits
                            templatenumpy[variantindex,nickindelcountnr] += 1
                        
                        if (sequenceWT_beforeflap != controlWT_beforeflap) and (sequenceMUT_beforeflap != controlMUT_beforeflap):  # check if 4bp window around nick has unintended edits
                            templatenumpy[variantindex,beforeflapindelcountnr] += 1
                            
                elif correction_type == 'Deletion':
                    sequenceWT = targetend[distanceto_wide_start+21-basesbeforenick:distanceto_wide_start+21+RTTlength+correction_length+baseafterflap]
                    sequenceWT_beforeflap = targetend[distanceto_wide_start+21-basesbeforenick:distanceto_wide_start+21+RTTlength+correction_length-2]  # full sequence until 2bp before flap
                    sequenceWT_nick = targetend[distanceto_wide_start+21-basesbeforenick:distanceto_wide_start+21+2]
                    
                    sequenceMUT = targetend[distanceto_wide_start+21-basesbeforenick:distanceto_wide_start+21+RTTlength+baseafterflap]
                    sequenceMUT_beforeflap = targetend[distanceto_wide_start+21-basesbeforenick:distanceto_wide_start+21+RTTlength-2]  # full sequence until 2bp before flap
                    sequenceMUT_nick = targetend[distanceto_wide_start+21-basesbeforenick:distanceto_wide_start+21+2]
                    
                    controlWT = templatenumpy[variantindex,wide_initial_targetnr][21-basesbeforenick:21+RTTlength+correction_length+baseafterflap]
                    controlWT_beforeflap = templatenumpy[variantindex,wide_initial_targetnr][21-basesbeforenick:21+RTTlength+correction_length-2]
                    controlWT_nick = templatenumpy[variantindex,wide_initial_targetnr][21-basesbeforenick:21+2]
                    
                    controlMUT = templatenumpy[variantindex,wide_mutated_targetnr][21-basesbeforenick:21+RTTlength+baseafterflap]
                    controlMUT_beforeflap = templatenumpy[variantindex,wide_mutated_targetnr][21-basesbeforenick:21+RTTlength-2]
                    controlMUT_nick = templatenumpy[variantindex,wide_mutated_targetnr][21-basesbeforenick:21+2]
        
                    if sequenceWT == controlWT:
                        templatenumpy[variantindex,uneditedcountnr] += 1 # uneditedcount is column 42
                    elif sequenceMUT == controlMUT:
                        templatenumpy[variantindex,editedcountnr] += 1
                    else:
                        templatenumpy[variantindex,indelcountnr] += 1
                        if (sequenceWT_nick != controlWT_nick) and (sequenceMUT_nick != controlMUT_nick):  # check if 4bp window around nick has unintended edits
                            templatenumpy[variantindex,nickindelcountnr] += 1
                        
                        if (sequenceWT_beforeflap != controlWT_beforeflap) and (sequenceMUT_beforeflap != controlMUT_beforeflap):  # check if 4bp window around nick has unintended edits
                            templatenumpy[variantindex,beforeflapindelcountnr] += 1
                            
                elif correction_type == 'Insertion':
                    
                    sequenceWT = targetend[distanceto_wide_start+21-basesbeforenick:distanceto_wide_start+21+RTTlength-correction_length+baseafterflap]
                    sequenceWT_beforeflap = targetend[distanceto_wide_start+21-basesbeforenick:distanceto_wide_start+21+RTTlength-correction_length-2]  # full sequence until 2bp before flap
                    sequenceWT_nick = targetend[distanceto_wide_start+21-basesbeforenick:distanceto_wide_start+21+2]
                    
                    sequenceMUT = targetend[distanceto_wide_start+21-basesbeforenick:distanceto_wide_start+21+RTTlength+baseafterflap]
                    sequenceMUT_beforeflap = targetend[distanceto_wide_start+21-basesbeforenick:distanceto_wide_start+21+RTTlength-2]  # full sequence until 2bp before flap
                    sequenceMUT_nick = targetend[distanceto_wide_start+21-basesbeforenick:distanceto_wide_start+21+2]
                    
                    controlWT = templatenumpy[variantindex,wide_initial_targetnr][21-basesbeforenick:21+RTTlength-correction_length+baseafterflap]
                    controlWT_beforeflap = templatenumpy[variantindex,wide_initial_targetnr][21-basesbeforenick:21+RTTlength-correction_length-2]
                    controlWT_nick = templatenumpy[variantindex,wide_initial_targetnr][21-basesbeforenick:21+2]
                    
                    controlMUT = templatenumpy[variantindex,wide_mutated_targetnr][21-basesbeforenick:21+RTTlength+baseafterflap]
                    controlMUT_beforeflap = templatenumpy[variantindex,wide_mutated_targetnr][21-basesbeforenick:21+RTTlength-2]
                    controlMUT_nick = templatenumpy[variantindex,wide_mutated_targetnr][21-basesbeforenick:21+2]
                    
        
                    if sequenceWT == controlWT:
                        templatenumpy[variantindex,uneditedcountnr] += 1 # uneditedcount is column 42
                    elif sequenceMUT == controlMUT:
                        templatenumpy[variantindex,editedcountnr] += 1
                    else:
                        templatenumpy[variantindex,indelcountnr] += 1
                        if (sequenceWT_nick != controlWT_nick) and (sequenceMUT_nick != controlMUT_nick):  # check if 4bp window around nick has unintended edits
                            templatenumpy[variantindex,nickindelcountnr] += 1
                        
                        if (sequenceWT_beforeflap != controlWT_beforeflap) and (sequenceMUT_beforeflap != controlMUT_beforeflap):  # check if 4bp window around nick has unintended edits
                            templatenumpy[variantindex,beforeflapindelcountnr] += 1
            
end = time.time()
print('Time for loop:',end-start)

#make again dataframe out of numpy array for easier saving as csv:
templatedf = pd.DataFrame(data = templatenumpy, 
                  index = templatedf.index.tolist(), 
                  columns = templatedf.columns.tolist())
  

# del final_diseasedf # remove final_diseasedf from memory
totalreads = templatedf["uneditedcount"] + templatedf["editedcount"] + templatedf["indelcount"]
templatedf['totalreads'] = totalreads
templatedf['totalreads'].replace(0, np.nan, inplace=True)
percentageedited = (templatedf["editedcount"]/templatedf['totalreads'])*100
templatedf['percentageedited'] = percentageedited
percentageunedited = (templatedf["uneditedcount"]/templatedf['totalreads'])*100
templatedf['percentageunedited'] = percentageunedited
percentageindels = (templatedf["indelcount"]/templatedf['totalreads'])*100
templatedf['percentageindel'] = percentageindels
# templatedf['barcodenr'] = templatedf.apply(lambda row: len(set(row.barcodes)), axis=1)

# fix spelling errors in dataframe columns:
templatedf['RTToverhanglength'] = templatedf['RToverhanglength']
templatedf['RTTlength'] = templatedf['RTlength']
templatedf['RTTseqoverhangrevcomp'] = templatedf['RTseqoverhangrevcomp']
templatedf['RTTrevcomp'] = templatedf['RTrevcomp']
templatedf['RTTmt'] = templatedf['RTmt']
templatedf['RTToverhangmt'] = templatedf['RToverhangmt']
templatedf['RTToverhangmatches'] = templatedf['RToverhangmatches']
templatedf['Spacer'] = templatedf['Protospacer-Sequence']
templatedf['PBSrevcomp'] = templatedf['PBSrevcomp13bp']
templatedf['spacermt'] = templatedf['protospacermt']
templatedf['spacermt'] = templatedf['protospacermt']



cols = ['Unnamed: 0.1',
 'Name',
 'position',
 'Original_Sequence',
 'Edited-Sequences',
 'Target-Strand',
 'Mutation_Type',
 'Correction_Type',
 'Correction_Length',
 'Editing_Position',
 'PBSlength',
 'RTToverhanglength',
 'RTTlength',
 'EditedAllele',
 'OriginalAllele',
 'Spacer',
 'PBSrevcomp',
 'RTTseqoverhangrevcomp',
 'RTTrevcomp',
 'pegRNA',
 'Editor_Variant',
 'spacermt',
 'extensionmt',
 'RTTmt',
 'RTToverhangmt',
 'PBSmt',
 'original_base_mt',
 'edited_base_mt',
 'original_base_mt_nan',
 'edited_base_mt_nan',
 'RTToverhangmatches',
 'PolyT_Proto_Extension',
 'wide_initial_target',
 'wide_mutated_target',
 'wide_mutated_target_intermediate',
 'protospacerlocation_only_initial',
 'PBSlocation',
 'RT_initial_location',
 'RT_mutated_location',
 'deepeditposition',
 'deepeditposition_intermediate',
 'barcode','barcode_rev',
 'specificprimerpart',
 'extension',
 'extension_15bp',
 'editedcount',
 'uneditedcount',
 'indelcount',
 'intermeidatecount',
 'intermediatepositions',
 'nickindelcount',
 'beforeflapindelcount',
 'totalreads',
 'percentageedited',
 'percentageunedited',
 'percentageindel']

templatedf = templatedf[cols]

final_diseasedf['match'] = final_diseasedf['match'].apply(lambda x: x[0])

print(len(final_diseasedf['match'].unique()))
print('Total reads:',len(diseasedf))
print('Successful assignment:',len(final_diseasedf))
print('Recombined:',diseasedf.recombination.sum())
print('Percentage recombined:',round((diseasedf.recombination.sum()/(diseasedf.recombination.sum()+len(final_diseasedf)))*100,2))

templatedf.to_csv('Analysis/20230426_'+shortname+'_analysisdf_focused.csv')