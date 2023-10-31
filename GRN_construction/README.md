# Preprocessing Data for Panda Tool

## Steps
1. Create t2t_TSS1000_TSS1000_cat.fasta from cat_gene_info.csv and cat_gene_sequence.fasta (building_a_regulation_prior_network_T2T.ipynb)
2. Create convPWMs, allTFs.meme and tfNames.txt from PWMs and Human_TF_MotifList_v_1.01.csv (building_a_regulation_prior_network_T2T.ipynb)
3. Run python main.py to get regMatpval0.000005.csv (main.py)
4. Run section 3 (few cells) and 4 to get TF_genes.txt (building_a_regulation_prior_network_T2T.ipynb)

## Input Files

- ROSMAP_all_counts_matrix.txt
- ROSMAP_biospecimen_metadata.csv
- ROSMAP_clinical.csv
- ROSMAP_assay_rnaSeq_metadata.csv
- TF_Genes.txt (Complete TF_genes)
- TF_TF_PPI.csv (Complete TF_TF_PPI)


## Output
- clinical_data_preprocessed.csv 
- count_matrix_mrn.csv
- TF_Genes.txt (Filtered TF_Genes)
- TF_PPI.csv (Filtered TF_PPI)



# Construct GRN

- git clone https://github.com/netZoo/netZooPy.git
- Place GRN_construction.py in netZooPy/netZooPy/ directory
- replace netZooPy/netZooPy/lioness/lioness.py with modified lioness.pu
- Execute python run_lioness.py

Warning: Change directory path 

