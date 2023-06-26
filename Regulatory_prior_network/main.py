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
    os.chdir('/Users/anwer/Desktop/Gene_Regulatory_Network/code/Regulatory_Prior_Network_Data/output/convPWMs')
    tfNames=[]

    with open("tfNames.txt", "r") as f:
        for line in f:
        tfNames.append(str(line.strip()))


    os.chdir(pathos)
    geneNames  = []
    fasta_sequences = SeqIO.parse(open('/Users/anwer/Desktop/Gene_Regulatory_Network/code/Regulatory_Prior_Network_Data/output/hg38_sequence_Tss1000_Tss1000.fasta'),'fasta')
    for fasta in tqdm(fasta_sequences):
        name, sequence = fasta.id, str(fasta.seq)
        geneNames.append(name)

    nTFs = len(tfNames)


if __name__ == '__main__':
    main()
