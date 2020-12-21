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
#     display_name: Python (enm)
#     language: python
#     name: enm
# ---

# %% [markdown] toc=true
# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#Compare-neighbor-betweenness" data-toc-modified-id="Compare-neighbor-betweenness-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Compare neighbor betweenness</a></span></li><li><span><a href="#Check-high-effectivenes-nodes-for-a-given-degree-to-understand-the-topology" data-toc-modified-id="Check-high-effectivenes-nodes-for-a-given-degree-to-understand-the-topology-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>Check high effectivenes nodes for a given degree to understand the topology</a></span><ul class="toc-item"><li><ul class="toc-item"><li><ul class="toc-item"><li><span><a href="#Sensor:-sps18-and-neighbors" data-toc-modified-id="Sensor:-sps18-and-neighbors-2.0.0.1"><span class="toc-item-num">2.0.0.1&nbsp;&nbsp;</span>Sensor: sps18 and neighbors</a></span></li><li><span><a href="#Sensor:-sps18-and-neighbors" data-toc-modified-id="Sensor:-sps18-and-neighbors-2.0.0.2"><span class="toc-item-num">2.0.0.2&nbsp;&nbsp;</span>Sensor: sps18 and neighbors</a></span></li></ul></li></ul></li></ul></li><li><span><a href="#Work-on-random-generated-graph" data-toc-modified-id="Work-on-random-generated-graph-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>Work on random generated graph</a></span></li></ul></div>

# %%
# %load_ext autoreload
# %autoreload 2
import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx
import numpy as np
import pandas as pd
import pickle
import os
import re
import itertools as itr
from src.visualize.visualize import plot_correlation_density, plot_vector, plot_lambda_collectivity
from src.enm import Enm

# %%
import statsmodels.api as sm


# %%
def plot_prs_subnetwork(gc, enm_obj, attr='eff',node_name='upc2', radius=2,**kwargs):
    if attr == 'eff':
        nx.set_node_attributes(enm_obj.graph_gc,dict(zip(enm_obj.nodes,enm_obj.prs_mat[enm_obj.nodes.index(node_name),:])),attr)
    else:
        nx.set_node_attributes(enm_obj.graph_gc,dict(zip(enm_obj.nodes,enm_obj.prs_mat[:,enm_obj.nodes.index(node_name)])),attr)
    g_sub = nx.ego_graph(enm_obj.graph_gc, node_name,radius)
    return plot_subnetwork(g_sub,attr,**kwargs)


# %%
def plot_subnetwork(g_sub,attr, **kwargs):
    import matplotlib.colors as colors
#['red'  if i==node_name else 'blue' for i in g_sub.nodes ]
    val_range = kwargs.pop('val_range',[val for i, val in g_sub.nodes(attr)])
    vmin = min(val_range)
    vmax = max(val_range)
    cmap = plt.cm.viridis_r
    figsize = kwargs.pop("figsize",(10,10))
    fig,ax = plt.subplots(figsize=figsize)
    #p = nx.spring_layout(g_sub)
    nx.draw_networkx_nodes(g_sub, ax=ax, node_color = [val for i, val in g_sub.nodes(attr)],
           cmap=cmap, vmin=vmin, vmax=vmax,**kwargs)
    if kwargs.pop('with_labels',False):
        nx.draw_networkx_labels(g_sub, ax=ax, **kwargs)
    nx.draw_networkx_edges(g_sub, ax=ax,alpha=0.1,**kwargs)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=colors.LogNorm(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = plt.colorbar(sm)
    title = kwargs.pop('node_name','')
    plt.title(title)
    filename = kwargs.pop('filename',None)
    if filename is not None:
        plt.savefig(filename,bbox_inches='tight')
    
    plt.show()
    #print(g_sub.nodes)
    return val_range



# %%
os.chdir('../')

# %%
data_path = 'data/interim/pcc_prs_0713/'
with open(f'data/interim/pcc.pickle','rb') as f:
    e_pcc = pickle.load(f)



# %%
e_pcc.figure_path = 'reports/figures/pcc_prs_0713/'

# %%
# e_pcc.output_path=data_path
# e_pcc.simulate_rewire(sim_num=2) 
# e_pcc.spring_pos()
# with open(f'data/interim/pcc.pickle','wb') as f:
#     pickle.dump(e_pcc,f)

# %%
e_pcc.plot_network_spring()

# %%
e_pcc.get_prs(no_diag=False)
df=e_pcc.df

# %%
df['norm_eff'] =  df.eff/df.deg

# %%
fig, ax = plt.subplots(1,2, figsize = (10,5))
df.plot.scatter('deg', 'norm_eff', ax=ax[0])
#ax[0].set_xlim(0,25)
ax[0].axhline(0.5,c='red')
df.plot.scatter('deg', 'eff', ax=ax[1])
#ax[1].set_xlim(0,25)


# %%
rand_df = e_pcc.e_list[1].df
rand_df['norm_eff'] =  rand_df.eff/rand_df.deg
fig, ax = plt.subplots(1,2, figsize = (10,5))
rand_df.plot.scatter('deg', 'norm_eff', ax=ax[0])
#ax[0].set_xlim(0,25)
ax[0].axhline(0.5,c='red')
rand_df.plot.scatter('deg', 'eff', ax=ax[1])
#ax[1].set_xlim(0,25)


# %%
fig, ax = plt.subplots()
plt.hist(df.loc[df.loc[:,'deg']==3,'norm_eff'],bins=100)
ax.set_xlabel('Effectiveness')

# %%
fig , ax = plt.subplots(1,2,figsize=(12,5))
df.plot.scatter('deg','btw',ax=ax[0], c = 'eff',colormap='viridis')
ax[0].set_ylim(0,0.04)
df.plot.scatter('eff','btw',ax=ax[1], c = 'deg',colormap='viridis')
ax[1].set_ylim(0,0.04)
plt.tight_layout()
#ax.plot(x, model.coef_*x+model.intercept_)
#ax.axvline(0.25, c='red')

# %%
fig , ax = plt.subplots(1,2,figsize=(12,5))
df.plot.scatter('deg','btw',ax=ax[0], c = 'sens',colormap='viridis')
ax[0].set_ylim(0,0.04)
df.plot.scatter('sens','btw',ax=ax[1], c = 'deg',colormap='viridis')
ax[1].set_ylim(0,0.04)
plt.tight_layout()
#ax.plot(x, model.coef_*x+model.intercept_)
#ax.axvline(0.25, c='red')

# %%
e_pcc.heatmap_annotated()

# %% [markdown]
# # Compare neighbor betweenness

# %%
nx.incidence_matrix(e_pcc.graph_gc,'ybr196c-a')


# %% [markdown]
# # Check high effectivenes nodes for a given degree to understand the topology

# %%
grouped = df.groupby('deg')#['eff'].quantile()
grouped.apply(lambda g: g[g['sens']>np.quantile(g['sens'],0.99)])#.orf_name.values

# %% [markdown]
# #### Sensor: sps18 and neighbors

# %%
vr = plot_subnetwork(e_pcc.graph_gc,e_pcc,attr='eff',node_name='hom3',radius=2,node_size=100,figsize=(8,8),with_labels=True)    
vr

# %%

# %%
df.loc[df.orf_name.isin(['ydr157w', 'rma1', 'agp1', 'bap2', 'stp1', 'hom2', 'stp2', 'hom3', 'ser1', 'ser2', 'dal81'])]

# %%
knn = nx.average_neighbor_degree(e_pcc.graph_gc)
df['knn']=knn.values()

# %% [markdown]
# #### Sensor: sps18 and neighbors

# %%
vr = plot_prs_subnetwork(e_pcc.graph_gc,e_pcc,attr='eff',node_name='sps18',radius=1,node_size=100,figsize=(8,8),with_labels=True)    
vr

# %%
df.loc[df.orf_name.isin(['sps18', 'ynl046w', 'pch2', 'tan1', 'yet2'])]

# %%
node_name = 'aad16'
nx.set_node_attributes(e_pcc.graph_gc,dict(zip(e_pcc.nodes,e_pcc.prs_mat[:,e_pcc.nodes.index(node_name)])),'eff')
g_sub = nx.ego_graph(e_pcc.graph_gc, node_name,2)#['red'  if i==node_name else 'blue' for i in g_sub.nodes ]
nx.draw(g_sub,node_color = [val for i, val in g_sub.nodes('eff')])

g_sub.nodes

# %%
g_sub.nodes('eff')

# %%
df.loc[df.orf_name.isin(['akr2', 'ypr039w','tip41','irc21'])]

# %%
from sklearn.linear_model import LinearRegression
x = df[['eff']]
y = df.btw.values
model = LinearRegression()
model.fit(x, y)

# %%
r_sq = model.score(x,y)

# %%
model.coef_

# %% [markdown]
# # Work on random generated graph
#
#

# %%
gc = nx.Graph()

# %%
gc.add_nodes_from(range(12))

# %%
edges = [(0,1),
        (0,2),
        (1,2),
        (1,3),
        #(1,13),
        (3,4),
        (3,5),
        (2,10),
        (10,9),
        (10,11),
        (0,7),
        (7,8),
        (7,6)]
gc.add_edges_from(edges)

# %%
nx.draw(gc,with_labels=True)

# %%
e_sm = Enm('small')
e_sm.G = gc
e_sm.giant_component()
e_sm.gnm_analysis()

# %%
e_sm.gnm_analysis()

# %%
e_sm.df.loc[:,['orf_name','deg','eff','sens']]

# %%
pos = nx.spring_layout(gc)

# %%
nx.set_node_attributes(e_sm.graph_gc,dict(zip(e_sm.df.orf_name,e_sm.df.eff)),'eff')
plot_subnetwork(e_sm.graph_gc,'eff',with_labels=True, pos=pos, val_range=[0.1,3.5],filename= 'reports/figures/sample_network_0729/eff.png')


# %%
nx.set_node_attributes(e_sm.graph_gc,dict(zip(e_sm.df.orf_name,e_sm.df.sens)),'sens')
plot_subnetwork(e_sm.graph_gc,'sens',with_labels=True,pos=pos,val_range=[0.1,3.5],filename= 'reports/figures/sample_network_0729/sens.png')

# %%
gc.remove_edge(3,5)
gc.add_edge(1,5)

e_sm = Enm('small')
e_sm.G = gc
e_sm.giant_component()
e_sm.gnm_analysis()

# %%
nx.set_node_attributes(e_sm.graph_gc,dict(zip(e_sm.df.orf_name,e_sm.df.eff)),'eff')
plot_subnetwork(e_sm.graph_gc,'eff',with_labels=True, pos=pos, val_range=[0.1,3.5],filename= 'reports/figures/sample_network_0729/eff_ch.png')


# %%
nx.set_node_attributes(e_sm.graph_gc,dict(zip(e_sm.df.orf_name,e_sm.df.sens)),'sens')
plot_subnetwork(e_sm.graph_gc,'sens',with_labels=True,pos=pos,val_range=[0.1,3.5],filename= 'reports/figures/sample_network_0729/sens_ch.png')

# %%
e_sm.df.loc[:,['orf_name','deg','eff','sens']]

# %%

# %%
#gc.remove_edge(3,5)
gc.add_edge(3,13)

e_sm = Enm('small')
e_sm.G = gc
e_sm.giant_component()
e_sm.gnm_analysis()

# %%
pos[13]=np.array([0.77917648,-0.86070557])

# %%
nx.set_node_attributes(e_sm.graph_gc,dict(zip(e_sm.df.orf_name,e_sm.df.eff)),'eff')
plot_subnetwork(e_sm.graph_gc,'eff',with_labels=True, pos=pos, val_range=[0.1,3.5],filename= 'reports/figures/sample_network_0729/eff_ch2.png')


# %%
nx.set_node_attributes(e_sm.graph_gc,dict(zip(e_sm.df.orf_name,e_sm.df.sens)),'sens')
plot_subnetwork(e_sm.graph_gc,'sens',with_labels=True,pos=pos,val_range=[0.1,3.5],filename= 'reports/figures/sample_network_0729/sens_ch2.png')

# %%
e_sm.df.loc[:,['orf_name','deg','eff','sens']]

# %%
e_pcc.df.loc[(e_pcc.df.deg<30 )&( e_pcc.df.eff<10)].plot.scatter('deg','eff')

# %%
nx.set_node_attributes(e_pcc.graph_gc,dict(zip(e_pcc.df.orf_name,e_pcc.df.eff)),'eff')
nx.set_node_attributes(e_pcc.graph_gc,dict(zip(e_pcc.df.orf_name,e_pcc.df.sens)),'sens')

# %%
val  = plot_subnetwork(e_pcc.graph_gc, pos=e_pcc.graph_gc.nodes('pos'),attr='sens', node_size=10)

# %%
e_pcc.plot_network_spring(node_color = e_pcc.df.sens.values,facecolor='white',figsize=(10,10),node_size=50)

# %%
from src.visualize.visualize import plot_network_spring

# %%

# %%
import seaborn as sns

# %%

s = plt.hist([i for i in e_pcc.df.sens if i <10],bins=101)

# %%

s = plt.hist([i for i in e_pcc.df.eff if i <10],bins=101)

# %%
