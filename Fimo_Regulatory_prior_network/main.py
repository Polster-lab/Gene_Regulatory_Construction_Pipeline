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
    os.chdir('../Regulatory_prior_Network_Data/input_for_fimo/convPWMs')
    tfNames=[]

    with open("tfNames.txt", "r") as f:
        for line in f:
            tfNames.append(str(line.strip()))


    os.chdir(pathos)
    geneNames  = []
    seq_path = '/Users/anwer/Desktop/Gene_Regulatory_Network_ROSMAP/code/Regulatory_prior_Network/Regulatory_prior_Network_Data/input_for_fimo/t2t_TSS1000_TSS1000_cat.fasta'
    fasta_sequences = SeqIO.parse(open('/Users/anwer/Desktop/Gene_Regulatory_Network_ROSMAP/code/Regulatory_prior_Network/Regulatory_prior_Network_Data/input_for_fimo/t2t_TSS1000_TSS1000_cat.fasta'),'fasta')
    for fasta in tqdm(fasta_sequences):
        name, sequence = fasta.id, str(fasta.seq)
        geneNames.append(name)

    nTFs = len(tfNames)
    pool = multiprocessing.Pool(6)
    res  = pool.map(partial(FIMO,output_dir = '../Regulatory_prior_Network_Data/fimo_output/',tfNames=tfNames,geneNames=geneNames,pqval=0,pathos = pathos,seq_path=seq_path), range(nTFs))
    res  = pd.concat(res)
    res.to_csv('../Regulatory_prior_Network_Data/fimo_output/regMatpval0.000005.csv')


if __name__ == '__main__':
    main()
