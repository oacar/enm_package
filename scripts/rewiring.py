import numpy as np
import networkx as nx
import pickle
import argparse
from enm.Enm import *
from enm.utils import *

def run_rewiring(pickle_file, n_sim, random_output_file):
    
    with open(pickle_file,'rb') as f:
        enm = pickle.load(f)
    print('rewiring')
    #Don't run with n>10, takes long
    enm.simulate_rewire(output_name='tmp',save=False, normalized=False,sim_num=n_sim)
    df = pd.DataFrame()
    for i in range(n_sim):
        df_tmp = enm.e_list[i].df
        df_tmp['random_id'] = i
        df = df.append(df_tmp)

    df.to_csv(random_output_file, index=False)
    #with open(snakemake.output['pickle_file'],'wb') as f:
    #    pickle.dump(enm,f, protocol=4)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run PRS analysis')
    parser.add_argument('-p', '--pickle_file', help='pickle file', required=True)
    parser.add_argument('-n', '--n_sim', help='number of simulations', required=True, type=int)
    parser.add_argument('-r', '--random_output_file', help='random output file', required=True)
    args = parser.parse_args()
    run_rewiring(args.pickle_file, args.n_sim, args.random_output_file)