import src
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx
import pickle
from src.enm import *
from src.utils import *

with open('data/interim/pcc_0909/pcc.pickle','rb') as f:
    e_pcc = pickle.load(f)

#e_pcc.df.orf_name.to_csv('data/interim/huri_gc_bg.tsv','\t', index=False, header=False)

#strain_ids = pd.read_csv('data/interim/strain_ids_with_experiment_count_all.csv')

#combined_df = pd.merge(e_pcc.df, strain_ids,left_on='orf_name',right_on='Allele Gene name')

goea, geneid2name = create_goea(gaf = 'data/raw/ontology/sgd.gaf', obo_fname='data/raw/ontology/go-basic.obo', background='data/interim/costanzo_gc_bg.tsv', sgd_info_tab = 'data/raw/ontology/SGD_features.tab')
sensors_pcc = e_pcc.df.loc[e_pcc.df.sens> np.quantile(e_pcc.df.sens,0.99)]
effectors_pcc = e_pcc.df.loc[e_pcc.df.eff> np.quantile(e_pcc.df.eff,0.99)]
degs_pcc = e_pcc.df.loc[e_pcc.df.deg> np.quantile(e_pcc.df.deg,0.99)]
query = e_pcc.df.loc[e_pcc.df.orf_name.isin(degs_pcc.orf_name),'Systematic gene name'].unique()
query_gene_ids = [key for key,value in geneid2name.items() if value in query]
goea_res_all = goea.run_study(query_gene_ids)
goea_res_sig = [r for r in goea_res_all if r.p_fdr_bh <0.05]
go_df_sensor = goea_to_pandas(goea_res_sig, geneid2name)

goea.wr_tsv('data/interim/pcc_0909/pcc_sensor_go.tsv', goea_res_sig)

from goatools.godag_plot import GODagPltVars

pltvars = GODagPltVars()
pltvars.fmthdr = "{GO}"

from goatools.godag_plot import plot_gos, plot_results, plot_goid2goobj

plot_results("reports/figures/pcc_0909/GO_PCC_{NS}.svg", goea_res_sig,
            GODagPltVars= pltvars,
            id2symbol = geneid2name,
        study_items=6,
        items_p_line=3)

gc = e_pcc.graph_gc.copy()
for i,j in gc.nodes(data=True):
    del j['pos']
nx.write_graphml(gc, 'data/interim/graph_files/pcc_full.graphml')


sensor_sub_pcc = get_subnetwork(e_pcc.graph_gc, sensors_pcc.orf_name.values)
nx.set_node_attributes(sensor_sub_pcc, dict(zip(e_pcc.df['orf_name'], e_pcc.df['Systematic gene name'])), 'orf_name') 
is_sensor = dict(zip(sensor_sub_pcc.nodes, [True if i in sensors_pcc.orf_name.values else False  for i in sensor_sub_pcc.nodes]))
is_go = dict(zip(sensor_sub_pcc.nodes, [True if i in go_df_sensor.study_items[7] else False  for val, i  in sensor_sub_pcc.nodes('orf_name')]))
nx.set_node_attributes(sensor_sub_pcc, is_go, 'is_go') 
nx.set_node_attributes(sensor_sub_pcc, is_sensor, 'is_sensor') 
nx.write_graphml(sensor_sub_pcc,'data/interim/graph_files/pcc_sensors.graphml')
