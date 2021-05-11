import src
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx
import pickle
from src.enm import *
from src.utils import *

with open('data/interim/huri_0916/huri.pickle','rb') as f:
    e_huri = pickle.load(f)

e_huri.df.orf_name.to_csv('data/interim/huri_gc_bg.tsv','\t', index=False, header=False)

goea, geneid2name = create_goea_human(gene2go ='data/raw/ontology/gene2go', obo_fname='data/raw/ontology/go-basic.obo', background='data/interim/huri_gc_bg.tsv')
df = e_huri.df
sensors_huri = e_huri.df.loc[e_huri.df.sens>np.quantile(e_huri.df.sens,0.99)]

query = sensors_huri.orf_name.unique()
query_gene_ids = [int(key) for key,value in geneid2name.items() if value in query]

goea_res_all = goea.run_study(query_gene_ids)
goea_res_sig = [r for r in goea_res_all if r.p_fdr_bh < 0.05]

go_df_sensor = goea_to_pandas(goea_res_sig, geneid2name)
goea.wr_tsv('data/interim/huri_0916/huri_go.tsv', goea_res_sig)

import mygene
mg = mygene.MyGeneInfo()
out = mg.querymany(e_huri.df.orf_name,scopes='ensemblgene', fields = 'symbol, unigene,entrezgene,uniprot', species=9606, as_dataframe=True)
ensembl_symbol_dict = out['symbol'].to_dict()
go_df_sensor['symbol'] = go_df_sensor.study_items.apply(lambda x: ', '.join([ensembl_symbol_dict[i] for i in x]))
go_df_sensor.to_csv('data/interim/huri_0916/huri_go_withnames.tsv','\t', index=False)
from goatools.godag_plot import GODagPltVars

pltvars = GODagPltVars()
pltvars.fmthdr = "{GO}"

from goatools.godag_plot import plot_gos, plot_results, plot_goid2goobj

plot_results("reports/figures/huri_0909/GO_HURI_{NS}.svg", goea_res_sig,
            GODagPltVars= pltvars,
            id2symbol = geneid2name,
        study_items=6,
        items_p_line=3)



sensor_sub_huri = get_subnetwork(e_huri.graph_gc, sensors_huri.orf_name.values)
nx.set_node_attributes(sensor_sub_huri, ensembl_symbol_dict, 'symbol') 
is_sensor = dict(zip(sensor_sub_huri.nodes, [True if i in sensors_huri.orf_name.values else False  for i in sensor_sub_huri.nodes]))
is_go = dict(zip(sensor_sub_huri.nodes, [True if i in go_df_sensor.study_items[0] else False  for i in sensor_sub_huri.nodes]))
nx.set_node_attributes(sensor_sub_huri, is_go, 'is_go') 
nx.set_node_attributes(sensor_sub_huri, is_sensor, 'is_sensor') 
nx.write_graphml(sensor_sub_huri,'data/interim/graph_files/huri_sensors.graphml')

