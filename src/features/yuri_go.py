import src
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx
import pickle
from src.enm import *
from src.utils import *

with open('data/interim/yuri_0916/yuri.pickle','rb') as f:
    e_yuri = pickle.load(f)

#e_yuri.df.orf_name.to_csv('data/interim/yuri_gc_bg.tsv','\t', index=False, header=False)
sgd  = pd.read_csv('data/raw/ontology/SGD_features.tab','\t',header=None)
pd.merge(sgd,e_yuri.df,  right_on='orf_name', left_on=3)[0].to_csv('data/interim/yuri_gc_bg.tsv',index=False, header=False)

goea, geneid2name = create_goea(gaf = 'data/raw/ontology/sgd.gaf', obo_fname='data/raw/ontology/go-basic.obo', background='data/interim/yuri_gc_bg.tsv', sgd_info_tab = 'data/raw/ontology/SGD_features.tab')

sensors_yuri = e_yuri.df.loc[e_yuri.df.sens> np.quantile(e_yuri.df.sens,0.99)]
query = sensors_yuri.orf_name.unique()
query_gene_ids = [key for key,value in geneid2name.items() if value in query]
goea_res_all = goea.run_study(query_gene_ids)
goea_res_sig = [r for r in goea_res_all if r.p_fdr_bh <0.05]
go_df_sensor = goea_to_pandas(goea_res_sig, geneid2name)

goea.wr_tsv('data/interim/yuri_0916/yuri_sensor_go.tsv', goea_res_sig)

from goatools.godag_plot import GODagPltVars

pltvars = GODagPltVars()
pltvars.fmthdr = "{GO}"

from goatools.godag_plot import plot_gos, plot_results, plot_goid2goobj

plot_results("reports/figures/yuri_0916/GO_YURI_{NS}.svg", goea_res_sig,
            GODagPltVars= pltvars,
            id2symbol = geneid2name,
        study_items=6,
        items_p_line=3, dpi=75)




sensor_sub_yuri = get_subnetwork(e_yuri.graph_gc, sensors_yuri.orf_name.values)
#nx.set_node_attributes(sensor_sub_huri, ensembl_symbol_dict, 'symbol') 
is_sensor = dict(zip(sensor_sub_yuri.nodes, [True if i in sensors_yuri.orf_name.values else False  for i in sensor_sub_yuri.nodes]))
is_go = dict(zip(sensor_sub_yuri.nodes, [True if i in go_df_sensor.study_items[0] else False  for i in sensor_sub_yuri.nodes]))
nx.set_node_attributes(sensor_sub_yuri, is_go, 'is_go') 
nx.set_node_attributes(sensor_sub_yuri, is_sensor, 'is_sensor') 
nx.write_graphml(sensor_sub_yuri,'data/interim/graph_files/yuri_sensors.graphml')

#fig , ax = plt.subplots(figsize=(6,6))
#pos = nx.spring_layout(sensor_sub_yuri,k=0.8,scale=0.1)
#nx.draw(sensor_sub_yuri , ax=ax , pos=pos,
#        node_color = ['orange' if j  is  True else 'white' for i,j in sensor_sub_yuri.nodes('is_sensor')],
#        labels=dict(sensor_sub_yuri.nodes('orf_name')), node_size=200,
#        edge_color='black', edgecolors='black', verticalalignment='bottom'
#        )
#fig.savefig('tmp.png')
