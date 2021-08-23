# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.5.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown] toc=true
# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#Import-packages" data-toc-modified-id="Import-packages-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Import packages</a></span></li><li><span><a href="#Read-Enm-pickle-object" data-toc-modified-id="Read-Enm-pickle-object-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>Read Enm pickle object</a></span></li><li><span><a href="#Calculate-node-positions-with-spring-layout" data-toc-modified-id="Calculate-node-positions-with-spring-layout-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>Calculate node positions with spring layout</a></span></li><li><span><a href="#Network-plot" data-toc-modified-id="Network-plot-4"><span class="toc-item-num">4&nbsp;&nbsp;</span>Network plot</a></span><ul class="toc-item"><li><span><a href="#Figure-1B,-left" data-toc-modified-id="Figure-1B,-left-4.1"><span class="toc-item-num">4.1&nbsp;&nbsp;</span>Figure 1B, left</a></span></li><li><span><a href="#Figure-1B,-right" data-toc-modified-id="Figure-1B,-right-4.2"><span class="toc-item-num">4.2&nbsp;&nbsp;</span>Figure 1B, right</a></span></li></ul></li><li><span><a href="#Figure-1C" data-toc-modified-id="Figure-1C-5"><span class="toc-item-num">5&nbsp;&nbsp;</span>Figure 1C</a></span></li><li><span><a href="#Figure-1D" data-toc-modified-id="Figure-1D-6"><span class="toc-item-num">6&nbsp;&nbsp;</span>Figure 1D</a></span></li><li><span><a href="#Figure-3C" data-toc-modified-id="Figure-3C-7"><span class="toc-item-num">7&nbsp;&nbsp;</span>Figure 3C</a></span></li><li><span><a href="#Figure-4B" data-toc-modified-id="Figure-4B-8"><span class="toc-item-num">8&nbsp;&nbsp;</span>Figure 4B</a></span></li><li><span><a href="#Figure-5D/F" data-toc-modified-id="Figure-5D/F-9"><span class="toc-item-num">9&nbsp;&nbsp;</span>Figure 5D/F</a></span></li></ul></div>

# %% [markdown]
# # Import packages

# %%
# %load_ext autoreload
# %autoreload 2
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx
import numpy as np
import pandas as pd
import pickle
import os
import re
import itertools as itr
from enm.Enm import Enm
from enm.utils import *

# %%
#os.chdir('../')
figure_folder = "reports/figures/paper_figures_052521"

# %% [markdown]
# # Read Enm pickle object

# %%
with open(snakemake.input.pickle_file_name,'rb') as f:
    e_pcc = pickle.load(f)



# %% [markdown]
# # Calculate node positions with spring layout

# %% [markdown]
# **This could be different from the ones shown in paper due to random node position calculation**
#
# This does not change any claims in the paper

# %%
#e_pcc.spring_pos(seed=12)

# %%
pos =e_pcc.graph_gc.nodes('pos')

# %% [markdown]
# # Network plot
# ## Figure 1B, left

# %%
fig, ax =plt.subplots(figsize=(5,5))

nx.draw_networkx_nodes(e_pcc.graph_gc,
                           node_size=0.2,
                           alpha=0.5,
                           node_color='k',
                           pos=pos,
                         ax=ax
                           # node_shape=matplotlib.markers.MarkerStyle(marker='o',fillstyle='full')
                           )
nx.draw_networkx_edges(e_pcc.graph_gc,
                           alpha= 0.1,
                           width= 0.1,
                           edge_color='k',
                           pos=pos,
                           label='PCC>0.2',ax=ax)
if snakemake.params['save']:
    plt.savefig(f'{figure_folder}/fig1b_left.png',bbox_inches='tight',dpi=150)

# %% [markdown]
# ## Figure 1B, right

# %%
e_pcc.figure_path=figure_folder
e_pcc.heatmap_annotated( save_figure=snakemake.params['save'] , figure_name = 'fig1b_right')

# %% [markdown]
# # Figure 1C

# %% [markdown]
# `get_clustered_nodes` finds the gene names belonging the outmost, smaller cluster for rows and columns as shown.

# %%
cc = e_pcc.df.iloc[e_pcc.get_clustered_nodes('column'),:]
rc = e_pcc.df.iloc[e_pcc.get_clustered_nodes('row'),:]

# %% [markdown]
# rc: row cluster
#
# cc: column cluster

# %%
df = e_pcc.df
df['eff_norm']=df['eff']/5183
df['sens_norm']=df['sens']/5183
df['cc'] = ['cc' if i in cc.orf_name.tolist() else 'other' for i in df.orf_name.tolist()]
df['rc'] = ['rc' if i in rc.orf_name.tolist() else 'other' for i in df.orf_name.tolist()]


# %%
def run_permutation_test(pooled,sizeZ,sizeY,delta):
    np.random.shuffle(pooled)
    starZ = pooled[:sizeZ]
    starY = pooled[-sizeY:]
    return starZ.mean() - starY.mean()

z = np.array(df.loc[df.rc=='rc','eff_norm'])
y = np.array(df.loc[df.rc!='rc','eff_norm'])
pooled = np.hstack([z,y])
delta = z.mean() - y.mean()
numSamples = 10000
estimates = np.fromiter(map(lambda x: run_permutation_test(pooled,z.size,y.size,delta),range(numSamples)),dtype=float)
diffCount = len(np.where(estimates >= delta)[0])
pval_rc = ((float(diffCount)+1)/(float(numSamples)+1))
pval_rc

# %%
z = np.array(df.loc[df.cc=='cc','sens_norm'])
y = np.array(df.loc[df.cc!='cc','sens_norm'])
pooled = np.hstack([z,y])
delta = z.mean() - y.mean()
numSamples = 10000
estimates = np.fromiter(map(lambda x: run_permutation_test(pooled,z.size,y.size,delta),range(numSamples)),dtype=float)
diffCount = len(np.where(estimates >= delta)[0])
pval_cc = ((float(diffCount)+1)/(float(numSamples)+1))
pval_cc

# %%
import seaborn as sns
color = {'cc': 'orange','other':'k','rc':'orangered'}
fig, ax = plt.subplots(1,2,figsize=(8,4))

sns.boxplot(data=df, x='rc',y='eff_norm',order=['rc','other'],palette=color,ax=ax[0])
ax[0].set_xticklabels(['Distinct\nRow\nCluster', 'Other genes'])
ax[0].set_xlabel('Row clusters')
ax[0].set_ylabel('Effectiveness (a.u.)')
ax[0].plot([0,0, 1,1], [0.006, 0.0065, 0.0065, 0.006], lw=1.5, c='k')
ax[0].text(.5, 0.0065, "***", ha='center', va='bottom', color='k',fontsize='large')
ax[0].set_ylim(-0.0001,0.007)

sns.boxplot(data=df, x='cc',y='sens_norm',order=['cc','other'],palette=color,ax=ax[1])
ax[1].set_xticklabels(['Distinct\nColumn\nCluster', 'Other genes'])
ax[1].set_xlabel('Columns clusters')
ax[1].set_ylabel('Sensitivity (a.u.)')
ax[1].plot([0,0, 1,1], [0.007, 0.0075, 0.0075, 0.007], lw=1.5, c='k')
ax[1].text(.5, 0.0075, "***", ha='center', va='bottom', color='k',fontsize='large')
ax[1].set_ylim(-0.0001,0.0081)

plt.tight_layout()

if snakemake.params['save']:
    plt.savefig(f'{figure_folder}/fig1c.pdf',bbox_inches='tight')

# %% [markdown]
# # Figure 1D

# %%
fig, ax = plt.subplots(figsize=(12,12))
#axs = ax.ravel()
legend_elements = [    ]

#for i in range(len(sensor_order)):
nx.draw_networkx_nodes(e_pcc.graph_gc, ax =ax , pos=pos,node_size=1, node_color='black')
nx.draw_networkx_edges(e_pcc.graph_gc, ax =ax , pos=pos,node_size=1, edge_color='black',alpha=0.2, width=0.1)
nx.draw_networkx_nodes(e_pcc.graph_gc, 
                       nodelist=rc.orf_name.tolist(),
                       ax=ax, 
                       pos=pos,
                       node_color='orangered',
                      edgecolors='black',
                      node_size=100,
                      node_shape='d')

ax.set_facecolor('white')
legend_elements.extend(
    [Line2D([0], [0], marker='d', color='black', label='Row cluster\ngenes',
                              markerfacecolor='orangered', markersize=50, linestyle="None"),
     Line2D([0], [0], marker='o', color='black', label='Other Genes',
                              markerfacecolor='black', markersize=25, linestyle="None"),
#                    Line2D([0], [0], marker='o', color='black', label='Effectors',
#                               markerfacecolor='black', markersize=10, linestyle="None"),
#                    Line2D([0], [0], marker='^', color='black', label='Sensors',
#                               markerfacecolor='black', markersize=10, linestyle="None"),
                       Line2D([0], [0], marker='o', color='black', label= 'High\nfunctional\nsimilarity',
                              markerfacecolor='black', markersize=0, linestyle="-", alpha=0.5,lw=10)
    ]
)
lgd = ax.legend(handles=legend_elements, fontsize=40,loc='center left', bbox_to_anchor=(1.1, 0.5))


#nx.draw_networkx_edges(nx.induced_subgraph(e_pcc.graph_gc, effector_pcc.orf_name.tolist()), ax=ax , pos=pos, edge_color='blue',alpha=0.5)
if snakemake.params['save']:
    plt.savefig(f'{figure_folder}/fig1d_up.png',bbox_inches='tight',dpi=150)

# %%
fig, ax = plt.subplots(figsize=(12,12))
#axs = ax.ravel()
legend_elements = [    ]

#for i in range(len(sensor_order)):
nx.draw_networkx_nodes(e_pcc.graph_gc, ax =ax , pos=pos,node_size=1, node_color='black')
nx.draw_networkx_edges(e_pcc.graph_gc, ax =ax , pos=pos,node_size=1, edge_color='black',alpha=0.2, width=0.1)
nx.draw_networkx_nodes(e_pcc.graph_gc, 
                       nodelist=cc.orf_name.tolist(),
                       ax=ax, 
                       pos=pos,
                       node_color='orange',
                      edgecolors='black',
                      node_size=100,
                      node_shape='D')

ax.set_facecolor('white')
legend_elements.extend([
    Line2D([0], [0], marker='D', color='black', label='Column cluster\ngenes',
                              markerfacecolor='orange', markersize=50, linestyle="None"),
    Line2D([0], [0], marker='o', color='black', label='Other Genes',
                              markerfacecolor='black', markersize=25, linestyle="None"),
                       Line2D([0], [0], marker='o', color='black', label= 'High\nfunctional\nsimilarity',
                              markerfacecolor='black', markersize=0, linestyle="-", alpha=0.5,lw=10),
    ]
)
lgd = ax.legend(handles=legend_elements, fontsize=40,loc='center left', bbox_to_anchor=(1.1, 0.5))


#nx.draw_networkx_edges(nx.induced_subgraph(e_pcc.graph_gc, effector_pcc.orf_name.tolist()), ax=ax , pos=pos, edge_color='blue',alpha=0.5)
if snakemake.params['save']:
    plt.savefig(f'{figure_folder}/fig1d_bottom.png',bbox_inches='tight',dpi=150)

# %% [markdown]
# # Figure 3C

# %%
change_go_group_names = snakemake.params['change_go_group_name']
sensor_go_rename = {
  "cellular response to iron ion starvation":'Iron ion transport' ,
"mitochondria-nucleus signaling pathway":  "Mitochondria-nucleus\nsignaling pathway\nand\nTricarboxylic acid cycle",
"phenylalanine transport":  "Phenylalanine transport",
"hexose metabolic process":  "Hexose metabolic process",
"tricarboxylic acid cycle":  "Tricarboxylic acid cycle"
}

# %%
sensors_pcc = pd.read_csv(snakemake.input.sensors_pcc)
sensor_colors = ["#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33"]
if change_go_group_names:
    sensors_pcc['go_group']=sensors_pcc['go_group'].map(sensor_go_rename,na_action='ignore')
sensor_order = sensors_pcc.groupby('go_group').sens.median().sort_values().index.tolist()

# %%
fig, ax = plt.subplots(figsize=(9,9))
#axs = ax.ravel()
legend_elements = [    ]

#for i in range(len(sensor_order)):
e_pcc.plot_network_spring(ax=ax,
                          node_size=1,
                          node_color='black',
 #                        node_size = [100 if i in sensors_pcc.orf_name.values or i in effector_pcc.orf_name.values else 1 for i in e_pcc.nodes],
                         #node_color = ['red' if i in sensors_pcc.orf_name.values else 'blue' if i in effector_pcc.orf_name.values else 'black' for i in e_pcc.nodes],
                         edge_color='black',savefig=False)
    #                         node_shape=['^' if i in sensors_pcc.orf_name.values else 'v' if i in effector_pcc.orf_name.values else 'o' for i in e_pcc.nodes])
    # nx.draw_networkx_nodes(e_pcc.graph_gc, nodelist=sensors_pcc.orf_name.values, node_size=200, pos=pos,
    #                           node_color='black',
    #                           node_shape='^',edgecolors='black',
    #                           linewidths=1)
nx.draw_networkx_nodes(nx.induced_subgraph(e_pcc.graph_gc, sensors_pcc.orf_name.tolist()),
                       pos=pos, 
                       node_color='black', alpha=1, node_shape='^')

for itr, i in enumerate(sensor_order):
    #print(i, effector_colors[itr])

    orf_names_to_plot = sensors_pcc.loc[sensors_pcc.go_group==i, 'orf_name'].tolist()
    nx.draw_networkx_nodes(e_pcc.graph_gc, nodelist=orf_names_to_plot, node_size=200, pos=pos,
                          node_color=sensor_colors[itr],
                          node_shape='^',edgecolors='black',
                          linewidths=1)
    legend_elements.append(
        Line2D([0], [0], marker='^', color='black', label=f'{i}',
                              markerfacecolor=sensor_colors[itr], markersize=30, linestyle="None")
    )
ax.set_facecolor('white')
legend_elements.append(
        Line2D([0], [0], marker='^', color='black', label=f'No GO Enrichment',
                              markerfacecolor='black', markersize=30, linestyle="None")
    )
legend_elements.extend(
    [Line2D([0], [0], marker='o', color='black', label='Other Genes',
                              markerfacecolor='black', markersize=10, linestyle="None"),
#                    Line2D([0], [0], marker='o', color='black', label='Effectors',
#                               markerfacecolor='black', markersize=10, linestyle="None"),
#                    Line2D([0], [0], marker='^', color='black', label='Sensors',
#                               markerfacecolor='black', markersize=10, linestyle="None"),
                       Line2D([0], [0], marker='o', color='black', label= 'High functional similarity',
                              markerfacecolor='black', markersize=0, linestyle="-", alpha=0.5, lw=10),
                   Line2D([0], [0], marker='o', color='red', label= 'Sensor-Sensor edges',
                              markerfacecolor='#018571', markersize=0, linestyle="-",lw=10)
                   #Line2D([0], [0], marker='o', color='blue', label= 'Effector-Effector edges',
    #                          markerfacecolor='#a6611a', markersize=0, linestyle="-")
    ]
)
#lgd = ax.legend(handles=legend_elements, fontsize=22,loc='center left', bbox_to_anchor=(1.1, 0.5),ncol=5)
nx.draw_networkx_edges(nx.induced_subgraph(e_pcc.graph_gc, sensors_pcc.orf_name.tolist()),pos=pos, edge_color='red', alpha=0.5)
if snakemake.params['save']:
    plt.savefig(f'{figure_folder}/fig3c.png',bbox_inches='tight',dpi=150)

# %%
fig = plt.figure()
figlegend = plt.figure(figsize=(3,1))
ax = fig.add_subplot(111)
#lines = ax.plot(range(10), plt.randn(10), range(10), plt.randn(10))
ax.axis('off')
lgd = ax.legend(handles=legend_elements, fontsize=40, loc='center',ncol=2)
if snakemake.params['save']:
    fig.savefig(f'{figure_folder}/fig3c_legend.png',bbox_inches='tight')

# %%
# antenna_legend = [legend_elements[i] for i in [0,2,3,4,5]]
# antenna_legend.extend(
#     [Line2D([0], [0], marker='o', color='black', label='Sensor to network connection node',
#                               markerfacecolor='black', markersize=30, linestyle="None"),
# #                    Line2D([0], [0], marker='o', color='black', label='Effectors',
# #                               markerfacecolor='black', markersize=10, linestyle="None"),
# #                    Line2D([0], [0], marker='^', color='black', label='Sensors',
# #                               markerfacecolor='black', markersize=10, linestyle="None"),
#                        Line2D([0], [0], marker='o', color='black', label= 'Sensor - nonsensor edges',
#                               markerfacecolor='black', markersize=0, linestyle="-", alpha=0.5, lw=10),
#     Line2D([0], [0], marker='o', color='red', label= 'Sensor - sensor edges',
#                               markerfacecolor='#018571', markersize=0, linestyle="-",lw=10)]
# )
# fig = plt.figure()
# figlegend = plt.figure(figsize=(3,1))
# ax = fig.add_subplot(111)
# #lines = ax.plot(range(10), plt.randn(10), range(10), plt.randn(10))
# ax.axis('off')
# lgd = ax.legend(handles=antenna_legend, fontsize=40, loc='center',ncol=2)
# if snakemake.params['save']:
#     fig.savefig(f'{figure_folder}/fig3f_legend.png',bbox_inches='tight')

# %% [markdown]
# # Figure 4B

# %%
effector_pcc = pd.read_csv(snakemake.input.effector_pcc)
effector_colors = ["#A65628", "#F781BF", "#999999"]
effector_order_orig = effector_pcc.groupby('go_group').eff.median().sort_values().index.tolist()

# %%
if change_go_group_names:
    effector_go_group_map = {
        effector_order_orig[0]:"Chromosome segregation",
        effector_order_orig[1]:"Golgi vesicle transport",
        effector_order_orig[2]:"Respiratory complex assembly"
    }
    effector_pcc['go_group'] = effector_pcc['go_group'].map(effector_go_group_map)
effector_order = effector_pcc.groupby('go_group').eff.median().sort_values().index.tolist()

# %%
fig, ax = plt.subplots(figsize=(10,10))
#axs = ax.ravel()
legend_elements = [    ]

#for i in range(len(sensor_order)):
e_pcc.plot_network_spring(ax=ax,
                          node_size=1,
                          node_color='black',
 #                        node_size = [100 if i in sensors_pcc.orf_name.values or i in effector_pcc.orf_name.values else 1 for i in e_pcc.nodes],
                         #node_color = ['red' if i in sensors_pcc.orf_name.values else 'blue' if i in effector_pcc.orf_name.values else 'black' for i in e_pcc.nodes],
                         edge_color='black',savefig=False)
    #                         node_shape=['^' if i in sensors_pcc.orf_name.values else 'v' if i in effector_pcc.orf_name.values else 'o' for i in e_pcc.nodes])
    # nx.draw_networkx_nodes(e_pcc.graph_gc, nodelist=sensors_pcc.orf_name.values, node_size=200, pos=pos,
    #                           node_color='black',
    #                           node_shape='^',edgecolors='black',
    #                           linewidths=1)

for itr, i in enumerate(effector_order):
   # print(i, effector_colors[itr])
    #print(itr)
    
    orf_names_to_plot = effector_pcc.loc[effector_pcc.go_group==i, 'orf_name'].tolist()
    nx.draw_networkx_nodes(e_pcc.graph_gc, nodelist=orf_names_to_plot, node_size=200, pos=pos,
                          node_color=effector_colors[itr],
                          node_shape='s',edgecolors='black',
                          linewidths=1)
    legend_elements.append(
        Line2D([0], [0], marker='s', color='black', label=f'{i}',
                              markerfacecolor=effector_colors[itr], markersize=20, linestyle="None")
    )
ax.set_facecolor('white')
legend_elements.extend(
    [Line2D([0], [0], marker='o', color='black', label='Other Genes',
                              markerfacecolor='black', markersize=10, linestyle="None"),
#                    Line2D([0], [0], marker='o', color='black', label='Effectors',
#                               markerfacecolor='black', markersize=10, linestyle="None"),
#                    Line2D([0], [0], marker='^', color='black', label='Sensors',
#                               markerfacecolor='black', markersize=10, linestyle="None"),
                       Line2D([0], [0], marker='o', color='black', label= 'High functional similarity',
                              markerfacecolor='black', markersize=0, linestyle="-", alpha=0.5,lw=10),
                   Line2D([0], [0], marker='o', color='blue', label= 'Effector-Effector edges',
                             markerfacecolor='#a6611a', markersize=0, linestyle="-",lw=10)
    ]
)
lgd = ax.legend(handles=legend_elements, fontsize=22,loc='center left', bbox_to_anchor=(1.1, 0.5))


nx.draw_networkx_edges(nx.induced_subgraph(e_pcc.graph_gc, effector_pcc.orf_name.tolist()), ax=ax , pos=pos, edge_color='blue',alpha=0.5)
if snakemake.params['save']:
    plt.savefig(f'{figure_folder}/fig4b.png',bbox_inches='tight',dpi=150)

# %% [markdown]
# # Figure 5D/F

# %%
#define source and target effector/sensor clusters
plot_paths = snakemake.params['plot_paths']
if plot_paths:
    eff_group = 'Chromosome segregation'
    sens_group = "Mitochondria-nucleus\nsignaling pathway\nand\nTricarboxylic acid cycle"
    sub_list = []
    #select source and target gene from respective clusters
    source = 'ctf4'
    target = 'rtg1'
    #calculate source and tartget
    l1 = e_pcc.get_prs_weighted_path(source,target)[1]
    sub_list.extend(l1)
    sub = nx.induced_subgraph(e_pcc.graph_gc, l1)
    node_sub=nx.induced_subgraph(sub,[i for i in l1 if i !=target])

    fig, ax = plt.subplots(figsize=(10,10))
    legend_elements = [    ]
    nx.draw_networkx_nodes(e_pcc.graph_gc, pos=pos, node_size=1, ax=ax, node_color='black')

    nx.draw_networkx_nodes(node_sub,pos=pos,alpha=0.8,
                        #  node_size = [prs_mat_df.loc[source,:].to_dict()[i]*10000 for i in sub.nodes],
        node_color = 'black')

    nx.draw_networkx_edges(sub,pos=pos)
    for itr, i in enumerate(sensor_order):
        if i ==sens_group:
            orf_names_to_plot = sensors_pcc.loc[sensors_pcc.go_group==i, 'orf_name'].tolist()
            sub_list.extend(orf_names_to_plot)

            nx.draw_networkx_nodes(e_pcc.graph_gc, nodelist=orf_names_to_plot, node_size=200, pos=pos,
                                node_color=sensor_colors[itr],
                                node_shape='^',edgecolors='black',
                                linewidths=1)
            nx.draw_networkx_edges(nx.induced_subgraph(e_pcc.graph_gc, orf_names_to_plot), ax=ax , pos=pos, edge_color='red',alpha=0.5)

            legend_elements.append(
                Line2D([0], [0], marker='^', color='black', label=f'Sensors ({i})',
                                    markerfacecolor=sensor_colors[itr], markersize=12, linestyle="None")
            )

    for itr, i in enumerate(effector_order):
        if i == eff_group:
            orf_names_to_plot = effector_pcc.loc[effector_pcc.go_group==i,'orf_name'].tolist()
            sub_list.extend(orf_names_to_plot)

            nx.draw_networkx_nodes(e_pcc.graph_gc, nodelist=orf_names_to_plot, node_size=200, pos=pos,
                                node_color=effector_colors[itr],
                                node_shape='s',edgecolors='black',
                                linewidths=1)
            nx.draw_networkx_edges(nx.induced_subgraph(e_pcc.graph_gc, orf_names_to_plot), ax=ax , pos=pos, edge_color='blue',alpha=0.5)
    ax.set_facecolor('white')
    if snakemake.params['save']:
        plt.savefig(f'{figure_folder}/fig5d.png',bbox_inches='tight',dpi=150)
        
    nx.write_edgelist(nx.induced_subgraph(e_pcc.graph_gc,sub_list),f'{figure_folder}/path1.csv', delimiter=',',data=False)

# %%
if plot_paths:
    eff_group = 'Respiratory complex assembly'
    sens_group = 'Iron ion transport'
    sub_list  = [ ]
    source = 'coa1'
    target = 'fet3'
    l1 = e_pcc.get_prs_weighted_path(source,target)[1]
    sub_list.extend(l1)
    sub = nx.induced_subgraph(e_pcc.graph_gc, l1)
    node_sub=nx.induced_subgraph(sub,[i for i in l1 if i !=target])
    #print(l1)
    fig, ax = plt.subplots(figsize=(10,10))
    legend_elements = [    ]
    nx.draw_networkx_nodes(e_pcc.graph_gc, pos=pos, node_size=1, ax=ax, node_color='black')
    nx.draw_networkx_nodes(node_sub,pos=pos,alpha=1,
                        #  node_size = [prs_mat_df.loc[source,:].to_dict()[i]*10000 for i in sub.nodes],
                        # node_shape = ['^' if i == target else 'o' for i in sub.nodes],
        node_color = [sensor_colors[2]  if i in ['sit1','ftr1'] else 'black' for i in node_sub.nodes])
    nx.draw_networkx_edges(sub,pos=pos)
    for itr, i in enumerate(sensor_order):
        #print(i, effector_colors[itr])
        if i ==sens_group:
            orf_names_to_plot = sensors_pcc.loc[sensors_pcc.go_group==i, 'orf_name'].tolist()
            sub_list.extend(orf_names_to_plot)
            nx.draw_networkx_nodes(e_pcc.graph_gc, nodelist=orf_names_to_plot, node_size=200, pos=pos,
                                node_color=sensor_colors[itr],
                                node_shape='^',edgecolors='black',
                                linewidths=1)
            nx.draw_networkx_edges(nx.induced_subgraph(e_pcc.graph_gc, orf_names_to_plot), ax=ax , pos=pos, edge_color='red',alpha=0.5)

            legend_elements.append(
                Line2D([0], [0], marker='^', color='black', label=f'Sensors ({i})',
                                    markerfacecolor=sensor_colors[itr], markersize=12, linestyle="None")
            )

    for itr, i in enumerate(effector_order):
        if i == eff_group:
            orf_names_to_plot = effector_pcc.loc[effector_pcc.go_group==i,'orf_name'].tolist()
            sub_list.extend(orf_names_to_plot)
            nx.draw_networkx_nodes(e_pcc.graph_gc, nodelist=orf_names_to_plot, node_size=200, pos=pos,
                                node_color=effector_colors[itr],
                                node_shape='s',edgecolors='black',
                                linewidths=1)
            nx.draw_networkx_edges(nx.induced_subgraph(e_pcc.graph_gc, orf_names_to_plot), ax=ax , pos=pos, edge_color='blue',alpha=0.5)

    ax.set_facecolor('white')
    if snakemake.params['save']:
        plt.savefig(f'{figure_folder}/fig5f.png',bbox_inches='tight',dpi=150)
        
    nx.write_edgelist(nx.induced_subgraph(e_pcc.graph_gc,sub_list),f'{figure_folder}/path2.csv', delimiter=',',data=False)

# %%
