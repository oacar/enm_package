# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.12.0
#   kernelspec:
#     display_name: Python (enm)
#     language: python
#     name: enm
# ---

# %% [markdown] toc=true
# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#Calculate-number-of-clusters-in-rewired" data-toc-modified-id="Calculate-number-of-clusters-in-rewired-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Calculate number of clusters in rewired</a></span></li></ul></div>

# %%
# %config Completer.use_jedi = False

# %%
# %pdb

# %%
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from enm.Enm import Enm, rewire_network
from enm.utils import *
import pickle
import random

#%%
#import argparse
import sys


args = sys.argv
seed = int(args[1])
print(seed)
random.seed(seed)        # or any integer
np.random.seed(seed)
# %%
class smake():
    def __init__(self,input_params):
        self.input = input_params
        self.output = {}


# %%
nw_type='coessentiality'

# %%
snakemake = smake({'pickle_file_name':f'../data/interim/{nw_type}/{nw_type}_enm_object_human.pickle',
                   'pcc_all':f'../data/interim/{nw_type}/{nw_type}_edgelist_human.csv',
                   'gaf':"../data/raw/ontology/goa_human.gaf",
                   'obo':'../data/raw/ontology/go-basic.obo', 
                   'background_file':f'../data/interim/{nw_type}/{nw_type}_go_background_list_human',
                   'sgd_info_tab':"../data/raw/ontology/human_name_id_map"})

# %%
snakemake.params = {'sim_num':50, 'save':True}

# %%
#strain_ids =pd.read_csv(snakemake.input.strain_ids)


# %%
e = Enm('rew')

# %%
e.read_network(snakemake.input['pcc_all'],sep=',')

# %%
gc_rew = rewire_network(e.graph_gc)
e_rew = Enm('rewired')
e_rew.G=gc_rew
e_rew.giant_component()
e_rew.gnm_analysis()
#e_rew.df = pd.merge(e_rew.df,strain_ids, left_on='orf_name',right_on='gene1')
e_rew.get_sensor_effector(True)

# %%
from utils_python.utils_python import create_goea

# %%
goea, geneid2name,a = create_goea(gaf=snakemake.input['gaf'], 
                                  obo_fname=snakemake.input['obo'],
                                  background=snakemake.input['background_file'],
                                 ev_exclude={'ND','IGI','HGI'},
                                 goset=['BP'], sgd_info_tab=snakemake.input['sgd_info_tab'],
                                  project_folder='./',methods=['fdr'], map_column_iloc=1, id_column_iloc=0)

# %%
e_rew.get_sensor_effector(quantile_threshold=0.99)

# %% tags=[]
# e_rew.analyze_components_biology(goea, geneid2name,True ,'orf_name','fdr')
# e_rew.analyze_components_biology(goea, geneid2name,False,'orf_name' , 'fdr')

# %%
sensors_pcc = e_rew.sensors_df

# %%
pos = e_rew.graph_gc.nodes('pos')

# %%
from matplotlib.lines import Line2D

fig, ax = plt.subplots(figsize=(9,9))
legend_elements = [    ]

e_rew.plot_network_spring(ax=ax,
                          node_size=1,
                          node_color='black',

                         edge_color='black',savefig=False)

nx.draw_networkx_nodes(nx.induced_subgraph(e_rew.graph_gc, sensors_pcc.orf_name.tolist()),
                       pos=pos, 
                       node_color='black', alpha=1, node_shape='^')

ax.set_facecolor('white')
# legend_elements.append(
#         Line2D([0], [0], marker='^', color='black', label=f'No GO Enrichment',
#                               markerfacecolor='black', markersize=30, linestyle="None")
#     )
legend_elements.extend(
    [Line2D([0], [0], marker='^', color='black', label='Sensors',
                              markerfacecolor='black', markersize=10, linestyle="None"),
     Line2D([0], [0], marker='o', color='black', label='Non-sensor Genes',
                              markerfacecolor='black', markersize=10, linestyle="None"),
#                    Line2D([0], [0], marker='o', color='black', label='Effectors',
#                               markerfacecolor='black', markersize=10, linestyle="None"),
                   
                       Line2D([0], [0], marker='o', color='black', label= 'High functional similarity',
                              markerfacecolor='black', markersize=0, linestyle="-", alpha=0.5, lw=10),
                   Line2D([0], [0], marker='o', color='red', label= 'Sensor-Sensor edges',
                              markerfacecolor='#018571', markersize=0, linestyle="-",lw=10)
                   #Line2D([0], [0], marker='o', color='blue', label= 'Effector-Effector edges',
    #                          markerfacecolor='#a6611a', markersize=0, linestyle="-")
    ]
)
lgd = ax.legend(handles=legend_elements, fontsize=22,
                loc='center right', bbox_to_anchor=(1.8, 0.5),ncol=1,
               frameon=False)
nx.draw_networkx_edges(nx.induced_subgraph(e_rew.graph_gc, sensors_pcc.orf_name.tolist()),pos=pos, edge_color='red', alpha=0.5)
ax.axis('off')
if snakemake.params['save']:
    plt.savefig(f'../reports/figures/paper_figures_{nw_type}/figs2d.png',bbox_inches='tight',dpi=150)

# %% [markdown]
# # Calculate number of clusters in rewired

# %% tags=[]
with open(snakemake.input['pickle_file_name'],'rb') as f:
    e_pcc = pickle.load(f)

# %% tags=[] jupyter={"outputs_hidden": true}
# %%capture
e_pcc.simulate_rewire(sim_num = snakemake.params['sim_num'])#

# %% tags=[]
from tqdm import tqdm
for i in tqdm(e_pcc.e_list):
    #i.df = pd.merge(i.df , strain_ids , left_on = 'orf_name', right_on='gene1')
    i.get_sensor_effector(quantile_threshold=0.99)
    i.analyze_components_biology(goea, geneid2name, True, 'orf_name', fdr_method='fdr')
    i.analyze_components_biology(goea, geneid2name,False, 'orf_name', fdr_method='fdr')


# %%
def save_rewired_data(e_pcc, idx, path):
    e = e_pcc.e_list[idx]
    df_rew = e.df
    sensors_df_rew = e.sensors_df
    effectors_df_rew = e.effectors_df
    g = e.graph_gc
    MYDIR = (f"{path}/{idx}")
    CHECK_FOLDER = os.path.isdir(MYDIR)

    # If folder doesn't exist, then create it.
    if not CHECK_FOLDER:
        os.makedirs(MYDIR)
        print("created folder : ", MYDIR)

    else:
        print(MYDIR, "folder already exists.")
    df_rew.to_csv(f"{path}/{idx}/df_rew_{idx}.csv",index=False)
    sensors_df_rew.to_csv(f"{path}/{idx}/sensors_df_rew_{idx}.csv",index=False)
    effectors_df_rew.to_csv(f"{path}/{idx}/effectors_df_rew_{idx}.csv",index=False)
    nx.write_edgelist(g, f"{path}/{idx}/g_rew_{idx}.edgelist.gz")


# %%
snakemake.output = {'rewired_data_folder':f"../data/interim/{nw_type}/{nw_type}_rewired_data_{seed}"}

# %%
snakemake.output['rewired_data_folder']

# %% tags=[]
import os
CHECK_FOLDER = os.path.isdir(snakemake.output['rewired_data_folder'])

# If folder doesn't exist, then create it.
if not CHECK_FOLDER:
    os.makedirs(snakemake.output['rewired_data_folder'])


# %% tags=[]
[save_rewired_data(e_pcc,i,snakemake.output['rewired_data_folder']) for i in range(snakemake.params['sim_num'])]

# %%
import glob

sensors_fls = glob.glob(f'{snakemake.output["rewired_data_folder"]}/sensors_df_rew*')

# %%
sensors_df_rew = pd.concat([pd.read_csv(f'{snakemake.output["rewired_data_folder"]}/{idx}/sensors_df_rew_{idx}.csv') for idx in range(snakemake.params['sim_num'])],keys=range(snakemake.params['sim_num'])).reset_index(level=0)

# %%
res = sensors_df_rew.groupby('level_0')['sensor_cluster'].nunique()
res2 = sensors_df_rew.groupby('level_0')['go_group'].nunique()

# %%
#res = [i.sensors_df.dropna(subset=['sensor_cluster']).sensor_cluster.nunique() for i in e_pcc.e_list]
#res2 = [i.sensors_df.dropna(subset=['go_group']).go_group.nunique() for i in e_pcc.e_list]

# %%
idxx = np.argwhere(np.array(res)>0).reshape(1,-1)[0]
idxx2 = np.argwhere(np.array(res2)>0).reshape(1,-1)[0]

# %% [markdown]
# There are 68 rewired networks (out of 500) with a sensor cluster

# %% [markdown]
# This corresponds to 13.6% of rewired networks

# %%
try:
    print(len(idxx)/len(res))
except Exception as e:
    print(len(idxx))

# %% [markdown]
# While 25% of these are GO enriched, that corresponds to only 3.4% of 500 cases

# %%
try:
    print(len(idxx2)/len(res2))
except Exception as e:
    print(len(idxx))

# %%
try:
    print(len(idxx2)/len(idxx))
except Exception as e:
    print(len(idxx))

# %%
rewired_nw_list = [nx.read_edgelist(f'{snakemake.output["rewired_data_folder"]}/{idx}/g_rew_{idx}.edgelist.gz') for idx in range(snakemake.params['sim_num'])]

# %%
# def plot_sensors(rew_list,idx ,sensors_df_list):
#     g =rew_list[idx]
#     sub_orfs =  sensors_df_list.loc[sensors_df_list.level_0==idx].dropna(subset=['sensor_cluster']).orf_name.tolist()
#     #g = e.graph_gc
#     induced_g = nx.induced_subgraph(g,sub_orfs)
#     sub_nw = get_subnetwork(g, sub_orfs, radius= 1)
#     pos_sub = nx.spring_layout(sub_nw)
#     fig, ax_ = plt.subplots() 
#     nx.draw_networkx_nodes(sub_nw,ax=ax_, pos = pos_sub, node_color = ['none' if i in sub_orfs else 'k' for i in sub_nw.nodes])
#     nx.draw_networkx_nodes(nx.induced_subgraph(sub_nw, sub_orfs), pos=pos_sub, node_shape='^', node_color='black')
#     nx.draw_networkx_edges(sub_nw,ax=ax_, pos = pos_sub)
#     nx.draw_networkx_edges(nx.induced_subgraph(sub_nw, sub_orfs), pos=pos_sub, edge_color='red')
#     plt.show()
#     #nx.draw_networkx(sub_nw)

# %%
# def find_antenna_sensors(rew_list,idx ,sensors_df_list):
#     g =rew_list[idx]
#     sub_orfs =  sensors_df_list.loc[sensors_df_list.level_0==idx].dropna(subset=['sensor_cluster']).orf_name.tolist()
#     #g = e.graph_gc
#     induced_g = nx.induced_subgraph(g,sub_orfs)
#     sub_nw = get_subnetwork(g, sub_orfs, radius= 1)
#     all_nodes = sub_nw.nodes
#     sensor_nodes = sub_orfs
#     num_nonsensor = len(np.setdiff1d(all_nodes, sensor_nodes)) / len(list(nx.connected_components(sub_nw)))
#     #print(num_nonsensor)
#     return num_nonsensor==1, len(np.setdiff1d(all_nodes, sensor_nodes))

# %% tags=[]
# for i in idxx:
#     print(i)
#     plot_sensors(rewired_nw_list,i,sensors_df_rew)

# %% [markdown]
# %96 of rewired networks with a sensor cluster have that sensor cluster as antenna

# %%
sensors_df = pd.read_csv(f'../data/interim/{nw_type}/{nw_type}_sensors_df_human.csv')

# %%
fig, ax =plt.subplots(figsize=(5,5))
ax.hist(res,3,color='navajowhite')
ax.set_xlabel('Number of sensor clusters')
ax.set_ylabel('Count')
ax.set_title('Number of sensor clusters\nat 500 rewired networks\ncompared to real network')
ax.axvline(sensors_df['sensor_cluster'].nunique(),c='darkblue',linestyle='-.')
plt.legend(handles = [
    Line2D([0],[0],color='navajowhite',linewidth=10,label='Rewired'),
    Line2D([0],[0],color='darkblue',linestyle='-.', label='Real')
], loc='upper center')
#if snakemake.params['save']:
 #   fig.savefig('reports/figures/paper_figures_supp/rewired_sensor_count.png', bbox_inches='tight',dpi=150)

# %%
fig, ax =plt.subplots(figsize=(5,5))
ax.hist(sensors_df_rew.groupby('level_0')['go_group'].nunique(),3,color='navajowhite')
ax.set_xlabel('Number of GO enriched sensor clusters')
ax.set_ylabel('Count')
ax.set_title(f'Number of GO enriched\nsensor clusters\nat {len(res)} rewired networks\ncompared to real network')
ax.axvline(sensors_df['go_group'].nunique(),c='darkblue',linestyle='-.')
ax.set_xlim(0,)
plt.legend(handles = [
    Line2D([0],[0],color='navajowhite',linewidth=10,label='Rewired'),
    Line2D([0],[0],color='darkblue',linestyle='-.', label='Real')
], loc='upper center')

if snakemake.params['save']:
    fig.savefig(f'../reports/figures/paper_figures_{nw_type}/rewired_go_sensor_count_{seed}.png', bbox_inches='tight',dpi=150)



# %% tags=[]
fig, ax =plt.subplots(figsize=(5,5))
ax.hist(res,3,color='navajowhite')
ax.set_xlabel('Number of sensor clusters')
ax.set_ylabel('Count')
ax.set_title(f'Number of sensor clusters\nat {len(res)} rewired networks\ncompared to real network')
ax.axvline(sensors_df['sensor_cluster'].nunique(),c='darkblue',linestyle='-.')
plt.legend(handles = [
    Line2D([0],[0],color='navajowhite',linewidth=10,label='Rewired'),
    Line2D([0],[0],color='darkblue',linestyle='-.', label='Real')
], loc='upper center')
if snakemake.params['save']:
    fig.savefig(f'../reports/figures/paper_figures_{nw_type}/rewired_sensor_count_{seed}.png', bbox_inches='tight',dpi=150)


# %%
#from IPython import embed; embed()


# %%
#rew_antenna = [find_antenna_sensors(rewired_nw_list,i,sensors_df_rew) for i in idxx]

# %%
# try:
#     print(sum([find_antenna_sensors(rewired_nw_list,i,sensors_df_rew) for i in idxx]) / len(idxx))
#     fig, ax =plt.subplots(figsize=(5,5))
#     ax.hist([j for i,j in rew_antenna if i==True],3,color='navajowhite')
#     ax.axvline(6,c='darkblue',linestyle='-.')
#     ax.set_xlabel('Number of antenna motifs')
#     ax.set_ylabel('Count')
#     ax.set_title('Number of antenna motifs\nat 68 rewired networks with sensor clusters\ncompared to real network')
#     plt.legend(handles = [
#         Line2D([0],[0],color='navajowhite',linewidth=10,label='Rewired'),
#         Line2D([0],[0],color='darkblue',linestyle='-.', label='Real')
#     ], loc='upper center')
# #    if snakemake.params['save']:
#  #       fig.savefig('reports/figures/paper_figures_supp/rewired_antenna_sensor_count.png', bbox_inches='tight',dpi=150)
# except Exception as e:
#     print('no clusters were found')
