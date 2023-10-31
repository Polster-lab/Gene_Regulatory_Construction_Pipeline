from netZooPy.panda.panda import Panda
from netZooPy.lioness import Lioness
import pickle


def pipeline():
    expression_data = '../GRN_input/count_matrix_mrn.csv'
    motif_data = '../GRN_input/TF_Genes.txt'
    ppi_data ='../GRN_input/TF_PPI.csv'
    panda_output = '/cephyr/NOBACKUP/groups/naiss2023-23-453/GRN_output_2/output_panda.txt'
    lioness_output ='/cephyr/NOBACKUP/groups/naiss2023-23-453/GRN_output_2/output_lioness.txt' 

    panda_obj = Panda(expression_data, motif_data, ppi_data, save_tmp=False,save_memory = False, remove_missing=False, keep_expression_matrix = True)
    with open('/cephyr/NOBACKUP/groups/naiss2023-23-453/GRN_output_2/panda_obj.pkl','wb') as outp:
        pickle.dump(panda_obj,outp, pickle.HIGHEST_PROTOCOL) 

    z = panda_obj.export_panda_results.sort_values(by=["force"], ascending=False)
    z.to_csv('/cephyr/NOBACKUP/groups/naiss2023-23-453/GRN_output_2/output_panda_df.csv', index = False) 
    lioness_obj = Lioness(panda_obj)
    lioness_obj.save_lioness_results(file='/cephyr/NOBACKUP/groups/naiss2023-23-453/GRN_output_2/lioness_total.csv')
    with open('/cephyr/NOBACKUP/groups/naiss2023-23-453/GRN_output_2/lioness_obj.pkl','wb') as outp2:
        pickle.dump(lioness_obj,outp2, pickle.HIGHEST_PROTOCOL)  



if __name__ == "__main__":
    pipeline()
