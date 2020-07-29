import src
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx
from src.enm import *
from src.utils import *
from Bio import SeqIO


figure_path = 'reports/figures/huri_0722'
huri = Enm('huri')
huri.read_network('data/interim/huri/huri_cleared.csv')
huri.gnm_analysis()

huri.figure_path=figure_path
huri.output_path = 'data/interim/huri/'
#huri.plot_collectivity()

huri.plot_network_spring()
#huri.plot_vector(sorted=True)
#huri.plot_scatter(x='deg',y='eff',figure_name='deg_eff',figure_extension='pdf')
#huri.plot_scatter(x='deg',y='sens',figure_name='deg_sens',figure_extension='pdf')

huri.simulate_rewire(output_name='rewire_data',save=True)
huri.rewire_df
huri.plot_correlation_density(x='eff',y='deg',figure_extension='pdf')
huri.plot_correlation_density(x='sens',y='deg',correlation='spearman',figure_extension='pdf')
nx.set_node_attributes(huri.graph_gc,dict(zip(list(huri.graph_gc.nodes),list(huri.graph_gc.nodes))),'orf_name')

with open(f"{huri.output_path}/huri.pickle",'wb') as f:
    pickle.dump(huri,f, protocol=4)

#go_df = pd.read_csv('data/interim/go_results/4.tsv','\t')

#huri.plot_network_spring(plot_go=True,go_df_list=[go_df],level_list=[0.2])
