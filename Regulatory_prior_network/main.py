from tqdm import tqdm
import matplotlib.pyplot as plt 
import numpy as np  
import os
from Bio import SeqIO         
import pandas as pd          
import multiprocessing        
from functools import partial
from utils import FIMO
def main():
    pathos=os.getcwd()
    os.chdir(pathos)
    os.chdir('./data/convPWMs')
    tfNames=[]

    with open("tfNames.txt", "r") as f:
        for line in f:
            tfNames.append(str(line.strip()))


    os.chdir(pathos)
    geneNames  = []
    fasta_sequences = SeqIO.parse(open('./data/hg38_sequence_Tss1000_Tss1000.fasta'),'fasta')
    for fasta in tqdm(fasta_sequences):
        name, sequence = fasta.id, str(fasta.seq)
        geneNames.append(name)

    nTFs = len(tfNames)
    pool = multiprocessing.Pool(8)
    res  = pool.map(partial(FIMO,output_dir = './data/fimo_output/',tfNames=tfNames,geneNames=geneNames,pqval=0,pathos = pathos), range(nTFs))
    res  = pd.concat(res)
    res.to_csv('regMatPval1e3.csv')


if __name__ == '__main__':
    main()
