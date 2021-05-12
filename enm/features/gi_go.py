import enm
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx
import pickle
from enm.enm import *
from enm.utils import *

with open('data/interim/gi_1205_0.2/gi_100.pickle','rb') as f:
    e_gi = pickle.load(f)

#e_gi.df.orf_name.to_csv('data/interim/huri_gc_bg.tsv','\t', index=False, header=False)

#strain_ids = pd.read_csv('data/interim/strain_ids_with_experiment_count_all.csv')

#combined_df = pd.merge(e_gi.df, strain_ids,left_on='orf_name',right_on='Allele Gene name')

goea, geneid2name = create_goea(gaf = 'data/raw/ontology/sgd.gaf', obo_fname='data/raw/ontology/go-basic.obo', background='data/interim/costanzo_gc_bg.tsv', sgd_info_tab = 'data/raw/ontology/SGD_features.tab')
sensors_gi = e_gi.df.loc[e_gi.df.sens> np.quantile(e_gi.df.sens,0.99)]
effectors_gi = e_gi.df.loc[e_gi.df.eff> np.quantile(e_gi.df.eff,0.99)]
degs_gi = e_gi.df.loc[e_gi.df.deg> np.quantile(e_gi.df.deg,0.99)]

query_deg = e_gi.df.loc[e_gi.df.orf_name.isin(degs_gi.orf_name),'Systematic gene name'].unique()
query_gene_ids = [key for key,value in geneid2name.items() if value in query_deg]
goea_res_all = goea.run_study(query_gene_ids)
goea_res_sig = [r for r in goea_res_all if r.p_fdr_bh <0.05]
go_df_deg = goea_to_pandas(goea_res_sig, geneid2name)

query_sens = e_gi.df.loc[e_gi.df.orf_name.isin(sensors_gi.orf_name),'Systematic gene name'].unique()
query_gene_ids = [key for key,value in geneid2name.items() if value in query_sens]
goea_res_all = goea.run_study(query_gene_ids)
goea_res_sig = [r for r in goea_res_all if r.p_fdr_bh <0.05]
go_df_sensor = goea_to_pandas(goea_res_sig, geneid2name)

query_eff = e_gi.df.loc[e_gi.df.orf_name.isin(effectors_gi.orf_name),'Systematic gene name'].unique()
query_gene_ids = [key for key,value in geneid2name.items() if value in query_eff]
goea_res_all = goea.run_study(query_gene_ids)
goea_res_sig = [r for r in goea_res_all if r.p_fdr_bh <0.05]
go_df_eff = goea_to_pandas(goea_res_sig, geneid2name)

goea.wr_tsv('data/interim/gi_1205_0.2/gi_sensor_go.tsv', goea_res_sig)

from goatools.godag_plot import GODagPltVars

pltvars = GODagPltVars()
pltvars.fmthdr = "{GO}"

from goatools.godag_plot import plot_gos, plot_results, plot_goid2goobj

plot_results("reports/figures/gi_0909/GO_PCC_{NS}.svg", goea_res_sig,
            GODagPltVars= pltvars,
            id2symbol = geneid2name,
        study_items=6,
        items_p_line=3)

gc = e_gi.graph_gc.copy()
for i,j in gc.nodes(data=True):
    del j['pos']
nx.write_graphml(gc, 'data/interim/graph_files/gi_full.graphml')


sensor_sub_gi = get_subnetwork(e_gi.graph_gc, sensors_gi.orf_name.values)
nx.set_node_attributes(sensor_sub_gi, dict(zip(e_gi.df['orf_name'], e_gi.df['Systematic gene name'])), 'orf_name') 
is_sensor = dict(zip(sensor_sub_gi.nodes, [True if i in sensors_gi.orf_name.values else False  for i in sensor_sub_gi.nodes]))
is_go = dict(zip(sensor_sub_gi.nodes, [True if i in go_df_sensor.study_items[7] else False  for val, i  in sensor_sub_gi.nodes('orf_name')]))
nx.set_node_attributes(sensor_sub_gi, is_go, 'is_go') 
nx.set_node_attributes(sensor_sub_gi, is_sensor, 'is_sensor') 
nx.write_graphml(sensor_sub_gi,'data/interim/graph_files/gi_sensors.graphml')
