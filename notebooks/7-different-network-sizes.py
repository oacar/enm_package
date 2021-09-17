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

# %%
import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from enm.utils import *

# %% [markdown]
# # Giant component size change

# %%
cpcc = pd.read_csv('../data/interim/costanzo_pcc_ALL')


# %%
def get_giant_component(G):
    Gc = max([G.subgraph(c).copy() for c in nx.connected_components(G)], key=len)
    return Gc


# %%
ne_ratio_list = []
nv_ratio_list = []

# %%
thr_list =  [0.05,0.1,0.2,0.25,0.3,0.35,0.4]

# %%
for i in np.arange(0.2,0.75,0.05):
    df = cpcc.loc[cpcc.pcc>=i]
    nw = nx.from_pandas_edgelist(df,source='gene1',target='gene2')
    nw_g = get_giant_component(nw)
    n_e = len(nw_g.edges)
    ne_ratio_list.append(n_e/len(nw.edges))
    n_v = len(nw_g.nodes)
    nv_ratio_list.append(n_v/len(nw.nodes))

    print(f"{i} : ne = {n_e}, ne_ratio: {n_e/len(nw.edges)}, nv = {n_v}, nv_ratio = {n_v/len(nw.nodes)}")

# %%
xlab = "Pearson's Correlation Coefficient threshold"

# %%
fig, ax = plt.subplots(1,1,figsize=(8,8))
ax.plot(np.arange(0.2,0.75,0.05),nv_ratio_list,'o-')
ax.set_ylabel('Giant component size to network size')
ax.set_xlabel(xlab)

# %% [markdown]
# # Effector and sensor clusters with GO enrichments

# %%
sensor_df_names = [ ]

# %%
effector_dfs = get_result_dfs('effectors_df', thr_list)

# %%
sensor_dfs =  get_result_dfs('sensors_df', thr_list)

# %%
effector_sensor_go_dfs = get_result_dfs('effector_sensor_combined_go_df',thr_list)

# %%
sensor_dfs[0.2].dropna(subset=['sensor_cluster']).loc[:,'sensor_cluster'].nunique()


# %%
def plot_go_thr_comparison(dfs, col, yaxis,plotname, xlab='PCC Threshold'):
    n_goterms = []
    rat_goterms = []
    n_clusters = []
    n_go_clusters = []
    for i in thr_list:
        df = dfs[i]
        n_goterms.append(df.dropna(subset=['go_group']).shape[0])
        rat_goterms.append(n_goterms[-1]/df.shape[0])
        n_clusters.append(df.dropna(subset=[col]).loc[:,col].nunique())
        n_go_clusters.append(df.dropna(subset=[col]).loc[:,'go_group'].nunique())
        
    fig, axs = plt.subplots(1,2,figsize=(5,2.5))
    axs[0].plot(thr_list, n_clusters, 'o-')
    axs[0].set_ylabel(f'Number of {yaxis} clusters', fontsize=12)
    axs[0].set_xlabel(xlab, fontsize=12)
    
    axs[1].plot(thr_list, n_go_clusters, 'o-')
    axs[1].set_ylabel(f'Number of go enriched\n{yaxis} clusters', fontsize=12)
    axs[1].set_xlabel(xlab, fontsize=12)
    
    plt.tight_layout()
    fig.savefig(f'../reports/figures/paper_figures_supp/{plotname}.png', bbox_inches='tight', dpi=150)
    # axs[2].plot(thr_list, [i/j if i!=0 else 0 for i,j in zip(n_go_clusters,n_clusters) ], 'o-')
    # axs[2].set_ylabel(f'% of go enriched {yaxis} clusters')
    # axs[2].set_xlabel(xlab)


# %%
plot_go_thr_comparison(sensor_dfs,'sensor_cluster', 'sensor',plotname='thr_num_sensor_clusters')

# %%
plot_go_thr_comparison(effector_dfs,'effector_cluster', 'effector',plotname='thr_num_effector_clusters')

# %%
#effector_sensor_go_dfs
go_overlap = pd.DataFrame({thr : [len(np.intersect1d(effector_sensor_go_dfs[thr].loc[effector_sensor_go_dfs[thr].cluster_type=='effector','GO'], effector_sensor_go_dfs[i].loc[effector_sensor_go_dfs[i].cluster_type=='effector','GO'])) for i in thr_list] for thr in sorted(thr_list)})
go_overlap.index = thr_list
def plot_heatmap_overlap(df, filename, title , cmap='YlGnBu' , figsize = (10,6)):
    import seaborn as sns
    mask = np.zeros_like(df)
    mask[np.tril_indices_from(mask, k=-1)] = True
    with sns.axes_style("white"):
        f, ax = plt.subplots(figsize=figsize)
        ax = sns.heatmap(100*df/ np.diag(df),annot=True, mask= mask, square=True, vmax=100, cbar=True, ax=ax, fmt='.1f',cmap=cmap)
        ax.set_title(title)
        
        for t in ax.texts: 
            t.set_text(t.get_text() + "%")
            t.set_size(9)
        ax.set_xlabel(ax.get_xlabel(),fontsize=16)
    plt.savefig(f'../reports/figures/paper_figures_supp/{filename}.png', bbox_inches='tight',dpi=150)    

plot_heatmap_overlap(go_overlap, 'go_overlap_test', 'GO term overlap for effectors\nin different thresholds', figsize=(5,5))


# %%

sensor_overlap = pd.DataFrame({thr : [len(np.intersect1d(sensor_dfs[thr].orf_name, sensor_dfs[i].orf_name)) for i in thr_list] for thr in sorted(thr_list)})
sensor_overlap.index = thr_list

effector_overlap = pd.DataFrame({thr : [len(np.intersect1d(effector_dfs[thr].orf_name, effector_dfs[i].orf_name)) for i in thr_list] for thr in sorted(thr_list)})
effector_overlap.index = thr_list

# %%
plot_heatmap_overlap(sensor_overlap, 'sensor_overlap', 'Sensor overlap for different thresholds', figsize=(5,5))

# %%
plot_heatmap_overlap(effector_overlap, 'effector_overlap', 'Effector overlap for different thresholds', figsize=(5,5))

# %%
