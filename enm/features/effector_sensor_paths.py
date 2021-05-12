from matplotlib.lines import Line2D
import qgrid
import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx
import numpy as np
import pandas as pd
import pickle
import os
import re
import itertools as itr
from enm.visualize.visualize import plot_correlation_density, plot_vector, plot_lambda_collectivity
from enm.enm import Enm
from enm.utils import *

with open('data/interim/pcc_021521_withclustering//pcc_10_r.pickle','rb') as f:
    e_pcc = pickle.load(f)
e_pcc.get_sensor_effector(use_threshold=True)
sensor_colors = [
"#bfd3e6",
"#9ebcda",
"#8c96c6",
"#8c6bb1",
"#88419d",
"#810f7c",
"#4d004b"]
effector_colors = ['#bae4b3',
'#74c476',
'#238b45']

sensors_pcc = e_pcc.df.loc[e_pcc.df.sens>np.quantile(e_pcc.df.sens,0.99)]
#sensor_sub_pcc = get_subnetwork(e_pcc.graph_gc, sensors_pcc.orf_name.values)
#is_sensor = dict(zip(sensor_sub_pcc.nodes, [True if i in sensors_pcc.orf_name.values else False  for i in sensor_sub_pcc.nodes]))
#nx.set_node_attributes(sensor_sub_pcc, is_sensor, 'is_sensor')
sensor_go_terms_separate = [
    'iron ion\ntransport',
    'phenylalanine\ntransport',
    None, 
    None,
    None,
    None,
    'tricarboxylic\nacid cycle', 
    'hexose metabolic\nprocess',
    'protein folding',
    'mitochondria-nucleus\nsignaling pathway/TCA',
    None,
    None,
    None,
    None,
    None,
    None,
    None,
]
sensor_connected_components = sorted([sorted(list(i)) for i in nx.connected_components(nx.induced_subgraph(e_pcc.graph_gc, sensors_pcc.orf_name))], key=lambda x: x[0])
sensor_go_terms= {}
for i in range(len(sensor_go_terms_separate)):
    for j in sensor_connected_components[i]:
        sensor_go_terms[j] = sensor_go_terms_separate[i]
sensors_pcc = sensors_pcc.merge(pd.DataFrame.from_dict(sensor_go_terms,orient='index',columns=['go_group']),right_index=True,left_on='orf_name')

effector_pcc = e_pcc.df.loc[e_pcc.df.eff>np.quantile(e_pcc.df.eff,0.99)]
#effector_sub_pcc = get_subnetwork(e_pcc.graph_gc, effector_pcc.orf_name.values)
#is_effector = dict(zip(effector_sub_pcc.nodes, [True if i in effector_pcc.orf_name.values else False  for i in effector_sub_pcc.nodes]))
#nx.set_node_attributes(effector_sub_pcc, is_effector, 'is_effector')
effectors_connected_comp = sorted([sorted(list(i)) for i in nx.connected_components(nx.induced_subgraph(e_pcc.graph_gc, effector_pcc.orf_name.tolist()))], key=lambda x:x[0])
effector_pcc.loc[:, 'go_group'] = ['mito' if i in effectors_connected_comp[0] else 'golgi' if i in effectors_connected_comp[1] else 'chromatin' for i in effector_pcc.orf_name]
effector_order = effector_pcc.groupby('go_group').eff.median().sort_values().index.tolist()
sensor_order = sensors_pcc.groupby('go_group').sens.median().sort_values().index.tolist()
source = 'aim10'
target = 'aft1'
e_pcc.get_prs(no_diag=False)
all_paths = []
for source in effector_pcc.orf_name:
    for target in sensors_pcc.dropna().orf_name:
        w, path = e_pcc.get_prs_weighted_path(source,target)
        print(w)
        all_paths.append(path)
flat_list = [item for sublist in all_paths for item in sublist]

all_path_nodes = np.unique(flat_list)
all_path_subnetworks = nx.induced_subgraph(e_pcc.graph_gc, all_path_nodes)

all_path_df = pd.DataFrame({'orf_name':all_path_nodes})
all_path_df['category'] = ['sensor' if i in sensors_pcc.orf_name.values else
                            'effector' if i in effector_pcc.orf_name.values else 
                             'path nodes' for i in all_path_nodes]
eff_sens_combined = effector_pcc.append(sensors_pcc)
all_path_df = pd.merge(all_path_df, eff_sens_combined[['orf_name','go_group']], how='left')
all_path_df['go_group']=[i  if pd.isna(i)==False else 'path nodes' for i in all_path_df.go_group.values]
folder = 'data/interim/eff_sens_paths/'
os.mkdir(folder)

all_path_df.to_csv(f"{folder}/all_paths_df.csv",index=None)
nx.to_pandas_edgelist(all_path_subnetworks).to_csv(f"{folder}/all_paths_network.csv",index=None)


#sensor_pos['ymr057c']=np.array([5.2,-1])
for i,eo in enumerate(effector_order):
    for j, so in enumerate(sensor_order):
        effector_go_group = eo
        sensor_go_group = so
        sensors_sub = nx.induced_subgraph(e_pcc.graph_gc, sensors_pcc.loc[sensors_pcc.go_group==sensor_go_group].orf_name)

        effectors_sub = nx.induced_subgraph(e_pcc.graph_gc, effector_pcc.loc[effector_pcc.go_group==effector_go_group].orf_name)

        sensor_pos, effector_pos, path_init_pos,sub, lmax, wmax = get_path_positions(e_pcc,sensors_sub,effectors_sub)

        fig, ax = plt.subplots()
        sens_plot = nx.draw_networkx_nodes(sensors_sub,pos=sensor_pos,ax=ax, node_size=50,node_color = sensor_colors[j])
        nx.draw_networkx_edges(sensors_sub,pos=sensor_pos,ax=ax)
        eff_plot = nx.draw_networkx_nodes(effectors_sub, pos = effector_pos , ax=ax,node_size=50, node_color=effector_colors[i] )
        nx.draw_networkx_edges(effectors_sub, pos = effector_pos , ax=ax)
        path_plot=nx.draw_networkx_nodes(sub,pos=path_init_pos,alpha=0.8,node_size=50,
                            #  node_size = [prs_mat_df.loc[source,:].to_dict()[i]*10000 for i in sub.nodes],
               node_color = 'black')
        nx.draw_networkx_edges(sub,pos=path_init_pos,alpha=0.8)
        path_plot.set_zorder(2)
        eff_plot.set_zorder(3)
        sens_plot.set_zorder(3)
        ax.set_xlim(0,10)
        ax.set_ylim(-5,15)
        plt.tight_layout()
        filename = "".join(x for x in f'{effector_go_group}_{sensor_go_group}' if x.isalnum())
        plt.savefig(f'reports/figures/eff_sens_path_plots_042721/{filename}.png',bbox_inches='tight')


fig, ax = plt.subplots(3,6,figsize=(24,10))
for i,eo in enumerate(effector_order):
    for j, so in enumerate(sensor_order):
        effector_go_group = eo
        sensor_go_group = so
        sensors_sub = nx.induced_subgraph(e_pcc.graph_gc, sensors_pcc.loc[sensors_pcc.go_group==sensor_go_group].orf_name)
        effectors_sub = nx.induced_subgraph(e_pcc.graph_gc, effector_pcc.loc[effector_pcc.go_group==effector_go_group].orf_name)
        sensor_pos, effector_pos, path_init_pos,sub = get_path_positions(sensors_sub,effectors_sub)
        sens_plot = nx.draw_networkx_nodes(sensors_sub,pos=sensor_pos,ax=ax[i][j], node_size=50,node_color = sensor_colors[j])
        nx.draw_networkx_edges(sensors_sub,pos=sensor_pos,ax=ax[i][j])
        eff_plot = nx.draw_networkx_nodes(effectors_sub, pos = effector_pos , ax=ax[i][j],node_size=50, node_color=effector_colors[i] )
        nx.draw_networkx_edges(effectors_sub, pos = effector_pos , ax=ax[i][j])
        path_plot=nx.draw_networkx_nodes(sub,pos=path_init_pos,alpha=0.8,node_size=50,ax=ax[i][j],
                            #  node_size = [prs_mat_df.loc[source,:].to_dict()[i]*10000 for i in sub.nodes],
               node_color = 'black')
        nx.draw_networkx_edges(sub,pos=path_init_pos,alpha=0.8,ax=ax[i][j])
        path_plot.set_zorder(2)
        eff_plot.set_zorder(3)
        sens_plot.set_zorder(3)
        ax[i][j].set_xlim(0,10)
        ax[i][j].set_ylim(-5,15)
plt.savefig(f'reports/figures/eff_sens_path_plots/combined.png',bbox_inches='tight')

goea , geneid2name = create_goea()
path_go_dict = {effector_order[0]:{},effector_order[1]:{},effector_order[2]:{}}
for i,eo in enumerate(effector_order):
    for j, so in enumerate(sensor_order):
        effector_go_group = eo#effector_order[0]
        sensor_go_group = so#sensor_order[0]
        sensors_sub = nx.induced_subgraph(e_pcc.graph_gc, sensors_pcc.loc[sensors_pcc.go_group==sensor_go_group].orf_name)
        effectors_sub = nx.induced_subgraph(e_pcc.graph_gc, effector_pcc.loc[effector_pcc.go_group==effector_go_group].orf_name)
        sensor_pos, effector_pos, path_init_pos,sub = get_path_positions(sensors_sub,effectors_sub)
        path_df = e_pcc.df.loc[e_pcc.df.orf_name.isin([i for i in sub.nodes][1:-1])]
        print(path_df['Systematic gene name'].tolist())
        #path_go_df = query_goatools(path_df, goea,geneid2name)
        #path_go_dict[eo][so]=path_go_df
