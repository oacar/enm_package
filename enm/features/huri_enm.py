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
huri.gnm_analysis(normalized=False)

huri.figure_path=figure_path
huri.output_path = 'data/interim/huri_1205/'
if os.path.exists(huri.output_path)==False:
    os.mkdir(huri.output_path)
#huri.plot_collectivity()
huri.spring_pos()
#huri.plot_network_spring()
#huri.plot_vector(sorted=True)
#huri.plot_scatter(x='deg',y='eff',figure_name='deg_eff',figure_extension='pdf')
#huri.plot_scatter(x='deg',y='sens',figure_name='deg_sens',figure_extension='pdf')

huri.simulate_rewire(output_name='rewire_data',save=True,normalized=False, simnum=10)
if os.path.exists(f"{huri.output_path}/random_dfs")==False:
    os.mkdir(f"{huri.output_path}/random_dfs")
random_dfs = [huri.e_list[i].df.to_csv(f'{huri.output_path}/random_dfs/rand_{i}.csv') for i in range(10)]
#huri.rewire_df
#huri.plot_correlation_density(x='eff',y='deg',figure_extension='pdf')
#huri.plot_correlation_density(x='sens',y='deg',correlation='spearman',figure_extension='pdf')
nx.set_node_attributes(huri.graph_gc,dict(zip(list(huri.graph_gc.nodes),list(huri.graph_gc.nodes))),'orf_name')

neigbor_btw = []
neighbor_degree = []
for i in huri.graph_gc.nodes:
    neigbor_btw.append(np.average([huri.df.loc[huri.df.orf_name==a,'btw'].values for a in nx.neighbors(huri.graph_gc,i)]))
    neighbor_degree.append(np.average([huri.df.loc[huri.df.orf_name==a,'deg'].values for a in nx.neighbors(huri.graph_gc,i)]))

huri.df['neighbor_btw'] = neigbor_btw
huri.df['neighbor_degree'] = neighbor_degree


with open(f"{huri.output_path}/huri.pickle",'wb') as f:
    pickle.dump(huri,f, protocol=4)

#go_df = pd.read_csv('data/interim/go_results/4.tsv','\t')

#huri.plot_network_spring(plot_go=True,go_df_list=[go_df],level_list=[0.2])
