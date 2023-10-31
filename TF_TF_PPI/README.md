# Constructing TF-Genes network (motif data) for T2T consortium 

## Steps
1. Create t2t_TSS1000_TSS1000_cat.fasta from cat_gene_info.csv and cat_gene_sequence.fasta (building_a_regulation_prior_network_T2T.ipynb)
2. Create convPWMs, allTFs.meme and tfNames.txt from PWMs and Human_TF_MotifList_v_1.01.csv (building_a_regulation_prior_network_T2T.ipynb)
3. Run python main.py to get regMatpval0.000005.csv (main.py)
4. Run section 3 (few cells) and 4 to get TF_genes.txt (building_a_regulation_prior_network_T2T.ipynb)

## Input Files

- cat_gene_sequence.fasta
- cat_gene_info.csv
- ROSMAP_all_counts_matrix.txt
- Directory: PWMs
- Human_TF_MotifList_v_1.01.csv


## Output
- t2t_TSS1000_TSS1000_cat.fasta : sequence[:1250] #since start is tss-1000, then we take the 2000bp upstream 
- allTFs.meme
- tfNames.txt
- regMatpval0.000005.csv
- **TF_Genes.txt** : TF-Motif data for GRN
- TF_symbol_to_ensembleid.csv
- Directory:  convPWMs

