from email import parser
import numpy as np
import networkx as nx
import pickle
import argparse
from enm.Enm import *
from enm.utils import *
import os
import sys
sys.setrecursionlimit(10000)

#figure_path = 'reports/figures/pcc_0603'

def run_prs(network_file, output_path, output_df, output_pickle, cluster_matrix=True, strain_ids_file=None):
    enm = Enm('enm')
    enm.read_network(network_file,sep=',')
    enm.gnm_analysis(normalized=False)
    enm.get_sensor_effector(use_threshold=True)
    print("calculating node positions")
    enm.spring_pos(seed=12)
    print("clustering prs matrix")
    if cluster_matrix:
        enm.cluster_matrix(enm.prs_mat)
    neighbor_degree = []
    for i in enm.graph_gc.nodes:
        neighbor_degree.append(np.average([enm.df.loc[enm.df.orf_name==a,'deg'].values for a in nx.neighbors(enm.graph_gc,i)]))

    #enm.df['neighbor_btw'] = neigbor_btw
    enm.df['neighbor_degree'] = neighbor_degree
    if strain_ids_file is not None and os.path.isfile(strain_ids_file):
        strain_ids = pd.read_csv(strain_ids_file)#snakemake.input['strain_ids_file'])
        enm.df = pd.merge(enm.df , strain_ids, left_on = 'orf_name', right_on='gene1')
    #enm.figure_path=figure_path
    enm.output_path = output_path
    #enm.get_category('data/interim/strain_ids_with_experiment_count_all.csv')

    print("saving files")
    with open(output_pickle,'wb') as f:
        pickle.dump(enm,f, protocol=4)

    enm.df.to_csv(output_df, index=True, index_label='orf_name_id')
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run PRS analysis')
    parser.add_argument('-i', '--network_file', help='network file', required=True)
    parser.add_argument('-o', '--output_path', help='output path', required=True)
    parser.add_argument('-d', '--output_df', help='output df', required=True)
    parser.add_argument('-p', '--output_pickle', help='output pickle', required=True)
    parser.add_argument('-c', '--cluster_matrix', help='cluster matrix', required=False, default=True)
    parser.add_argument('-s', '--strain_ids_file', help='strain ids file', required=False, default=None)
    args = parser.parse_args()
    run_prs(args.network_file, args.output_path, args.output_df, args.output_pickle, args.cluster_matrix, args.strain_ids_file)
    #run_prs(snakemake.input['network_file'], snakemake.params.output_path, snakemake.params.output_df, snakemake.params.output_pickle, snakemake.params.cluster_matrix, snakemake.input['strain_ids_file'])


