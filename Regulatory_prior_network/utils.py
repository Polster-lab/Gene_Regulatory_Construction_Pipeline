from tqdm import tqdm
import matplotlib.pyplot as plt 
import numpy as np  
import os
from Bio import SeqIO         
import pandas as pd          
import multiprocessing        
from functools import partial


def reduceSequence(sequence):
    seq=sequence[:2000]
    return seq



def FIMO(tfi, output_dir, geneNames, tfNames, pqval, pathos):
    os.chdir(pathos)
    os.chdir(output_dir)
    regMatQval = np.zeros((1,len(geneNames))) 
    regMatQ = pd.DataFrame(data = regMatQval, columns = geneNames, index = [tfNames[tfi]])
    
    bashCommand = '/Users/anwer/meme/bin/fimo --thresh 1e-5 --no-qvalue --verbosity 1 --max-stored-scores 100000000 --oc '+tfNames[tfi]+' --motif ' + tfNames[tfi] + ' /Users/anwer/Desktop/Gene_Regulatory_Network_ROSMAP/code/Regulatory_Prior_Network_Data/output/ConvPWMs/allTFs.meme  /Users/anwer/Desktop/Gene_Regulatory_Network_ROSMAP/code/Regulatory_Prior_Network_Data/output/hg38_sequence_Tss1000_Tss1000.fasta'
    res = os.system(bashCommand)
    if res != 0:
        print('could not find fimo')
    
    else:
        
        os.chdir(pathos)
        os.chdir(output_dir + tfNames[tfi])
        
    try:
        tf = pd.read_csv('fimo.tsv',sep='\t', comment='#')
        tf.iloc[:,3] = 1000 - tf.iloc[:,3]
        for geneName in geneNames:
            jset = np.where(np.in1d(tf.iloc[:,2],geneName))
            if jset[0].size != 0:
                regMatQ[geneName].iloc[0] = tf.iloc[jset[0],7].min()
                
    except: 
        print('empty fimo result')
    
    return regMatQ