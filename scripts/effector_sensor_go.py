import argparse
import numpy as np
import networkx as nx
import pickle
import random
from enm.Enm import *
from enm.utils import *
from utils_python.utils_python import create_goea

#snakemake = {'input': {'gaf': '../data/raw/ontology/sgd.gaf' , 'obo':'../data/raw/ontology/go-basic.obo', "background_file" :'../data/raw/ontology/sgd_costanzogenes', "sgd_info": '../data/raw/ontology/SGD_features.tab','pickle_file_name':'../data/interim_0.2/pcc.pickle'}}
# gaf = snakemake.input['gaf']
# obo = snakemake.input['obo']
# background_file = snakemake.input['background_file']
# sgd_info = snakemake.input['sgd_info']
# pickle_file = snakemake.input['pickle_file_name']
# output_file=snakemake.output['effector_sensor_go_df']

def run_effector_sensor_go(pickle_file, gaf,obo,background_file,sgd_info, sensors_df_fname, effectors_df_fname, effector_sensor_go_df_fname, systematic_gene_name, **kwargs):
    np.random.seed(0)
    random.seed(0)
    with open(pickle_file,'rb') as f:
        enm = pickle.load(f)
    if len(enm.graph_gc.nodes) <3000:
        enm.get_sensor_effector(use_threshold=True, quantile_threshold=0.95)
    else:
        enm.get_sensor_effector(use_threshold=True, quantile_threshold=0.99)
    print('go analysis')
    goea, geneid2name,obodag = create_goea(gaf = gaf, obo_fname=obo, 
                                    background=background_file, sgd_info_tab = sgd_info,
                                    goset=['BP'],
                                    ev_exclude={'ND','IGI','HGI'},
                                    methods = ['fdr'], **kwargs)

    print('get effector sensor go')         
    enm.analyze_components_biology(goea, geneid2name, True, systematic_gene_name)
    enm.analyze_components_biology(goea, geneid2name, False, systematic_gene_name)

    sensors_df = enm.sensors_df.reset_index(drop=True)
    sensors_df['sensor_cluster'] = sensors_df.fillna(value={'sensor_cluster':'Unclustered'}).sensor_cluster.astype("category")

    order = sensors_df.groupby('sensor_cluster').median()['sens'].sort_values().index.values
    sensors_df['sensor_cluster'].cat.reorder_categories(list(order),ordered=True, inplace=True)
    sensors_df['gid'] = sensors_df.groupby('sensor_cluster').ngroup()
    sensors_df['cluster_or_go'] = ['Unclustered' if sensors_df.at[i,'sensor_cluster'] == 'Unclustered'
                                else sensors_df.at[i,'gid'] if pd.isna(sensors_df.at[i,'go_group']) 
                                else sensors_df.at[i,'go_group'] for i in range(sensors_df.shape[0])]
    di = {"cellular response to iron ion starvation":"Iron ion\ntransport" ,
    "mitochondria-nucleus signaling pathway":"Mitochondria nucleus\nsignaling pathway",
    "phenylalanine transport":"Phenylalanine\ntransport",
    "hexose metabolic process":"Hexose metabolic\nprocess",
    "tricarboxylic acid cycle":"Tricarboxylic acid\ncycle"}
    sensors_df.replace({'cluster_or_go': di},inplace=True)

    sensors_df['label'] = ["Unclustered" if sensors_df.at[i,'cluster_or_go']=='Unclustered' 
                                        else f"SC{sensors_df.at[i,'gid']}" if type(sensors_df.at[i,'cluster_or_go'])==np.int64
                                        else f"SC{sensors_df.at[i,'gid']}\n{sensors_df.at[i,'cluster_or_go']}" for i in range(sensors_df.shape[0])]

    sensors_df.to_csv(sensors_df_fname,index=False)
    enm.effectors_df.to_csv(effectors_df_fname,index=False)

    effector_sensor_combined_go_df = pd.DataFrame()
    for i,j in enm.go_groups['sensors_go_groups'].items():
        if j is not None:
            j['cluster_id']=i
            j['cluster_type']='sensor'
            effector_sensor_combined_go_df = effector_sensor_combined_go_df.append(j)
    for i,j in enm.go_groups['effectors_go_groups'].items():
        if j is not None:
            j['cluster_id']=i
            j['cluster_type']='effector'
            effector_sensor_combined_go_df = effector_sensor_combined_go_df.append(j)
    effector_sensor_combined_go_df.to_csv(effector_sensor_go_df_fname,index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run effector sensor go analysis')
    parser.add_argument('--pickle_file', type=str, help='pickle file')
    parser.add_argument('--gaf', type=str, help='gaf file')
    parser.add_argument('--obo', type=str, help='obo file')
    parser.add_argument('--background_file', type=str, help='background file')
    parser.add_argument('--name_id_map', type=str, help='sgd info file', required=False, default=None)
    parser.add_argument('--sensors_df_fname', type=str, help='sensors df fname')
    parser.add_argument('--effectors_df_fname', type=str, help='effectors df fname')
    parser.add_argument('--effector_sensor_go_df_fname', type=str, help='effector sensor go df fname')
    parser.add_argument("--map_column_iloc", type=int, help="map column id", required=False, default=None)
    parser.add_argument("--id_column_iloc", type=int, help="id column id", required=False, default=None)
    args = parser.parse_args()
    kwargs ={}
    if args.map_column_iloc is not None and args.id_column_iloc is not None:
        kwargs = {'map_column_iloc':args.map_column_iloc, 'id_column_iloc':args.id_column_iloc}

    if 'costanzo' in args.pickle_file:
        systematic_gene_name ='Systematic gene name'
    else:
        systematic_gene_name ='orf_name'
    run_effector_sensor_go(args.pickle_file, args.gaf, args.obo, args.background_file, args.name_id_map, args.sensors_df_fname, args.effectors_df_fname, args.effector_sensor_go_df_fname, systematic_gene_name, **kwargs)
