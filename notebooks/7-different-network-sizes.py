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
plot_go_thr_comparison(effector_dfs,'effector_cluster', 'effector')

# %%
effector_sensor_go_dfs[0.2]

# %%
help(t.set_text)

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

sensor_overlap = pd.DataFrame({thr : [len(np.intersect1d(sensor_dfs[thr].orf_name, sensor_dfs[i].orf_name)) for i in thr_list] for thr in sorted(thr_list)})
sensor_overlap.index = thr_list

effector_overlap = pd.DataFrame({thr : [len(np.intersect1d(effector_dfs[thr].orf_name, effector_dfs[i].orf_name)) for i in thr_list] for thr in sorted(thr_list)})
effector_overlap.index = thr_list

mask1 = np.zeros_like(sensor_overlap)
mask1[np.tril_indices_from(mask1, k=-1)] = True
mask2 = np.zeros_like(effector_overlap)
mask2[np.tril_indices_from(mask2,k=-1)] = True
with sns.axes_style("white"):
    f, ax = plt.subplots(1,3,figsize=(16, 6), gridspec_kw={'width_ratios': [10,10, 1]})
    ax[0] = sns.heatmap(100*sensor_overlap / np.diag(sensor_overlap),annot=True, mask=mask1, square=True, vmax=100, cbar=False, ax=ax[0], fmt='.1f',cmap="YlGnBu")
    ax[0].set_title('Sensor overlap for different thresholds')
    for t in ax[0].texts: t.set_text(t.get_text() + "%")
    ax[1] = sns.heatmap(100 * effector_overlap / np.diag(effector_overlap),annot=True, mask=mask2, square=True, ax=ax[1],cbar_ax=ax[2], vmax=100, fmt='.1f',cmap="YlGnBu")
    for t in ax[1].texts: t.set_text(t.get_text() + "%")
    ax[1].set_title('Effector overlap for different thresholds')
    ax[0].text(-0.1, -1, 'A', fontsize=16, fontweight='bold',va='top',ha='right')
    ax[1].text(-0.1, -1, 'B', fontsize=16, fontweight='bold',va='top',ha='right')


plt.tight_layout()
plt.savefig('../reports/figures/paper_figures_supp/figs3.png', bbox_inches='tight',dpi=150)

# %%
import itertools as itr


# %%
def get_pairwise_intersect(dfs, thr_list):
    pairs = itr.combinations(thr_list,2)
    nt = lambda a, b: [len(np.intersect1d(dfs[a].orf_name,effector_dfs[b].orf_name)) , len(dfs[b].orf_name)]
    res = dict([ (t, nt(*t)) for t in pairs ])

    return res


# %% jupyter={"outputs_hidden": true}
get_pairwise_intersect(effector_dfs, thr_list)

# %% jupyter={"outputs_hidden": true}
get_pairwise_intersect(sensor_dfs, thr_list)

# %%
pcc_dfs = get_result_dfs('pcc_df', thr_list,0.2)

# %%
pcc_dfs[i].corr().loc['deg','eff']

# %%
[pcc_dfs[i].corr('pearson').loc['deg','sens'] for i in thr_list]

# %% jupyter={"outputs_hidden": true}
for i in thr_list:
    fig, ax = plt.subplots(figsize=(5,5))
    pcc_dfs[i].plot.scatter('deg','eff',ax=ax)
    ax.set_title(i)

# %% jupyter={"outputs_hidden": true}
for i in thr_list:
    fig, ax = plt.subplots(figsize=(5,5))
    pcc_dfs[i].plot.scatter('deg','sens',ax=ax)
    ax.set_title(i)

# %% jupyter={"outputs_hidden": true}
sensor_dfs[0.3]

# %%

# %% [raw]
# %config Completer.use_jedi = False

# %%
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import seaborn as sns



# %% [markdown]
# # Sensor change in different thresholds

# %%
all_df = pd.read_csv('../data/raw/Data File S3. Genetic interaction profile similarity matrices/cc_ALL.txt','\t',index_col=[0,1], header=[0,1])

# %%
df_long_renamed = all_df.melt(ignore_index=False, col_level=0).reset_index(inplace=False).rename(columns = {'level_0':'gene1','level_1':'Systematic gene name', 'variable':'gene2','value':'pcc'})


# %%
print("PCC threshold","Sensors in the higher threshold network","Sensors in the higher threshold GC","Number of sensors")
for i in range(1,len(thr_list)):
    prev = thr_list[i-1]
    curr = thr_list[i]
    
    sensor_names = sensor_dfs[prev].orf_name
    all_edges = df_long_renamed.loc[(df_long_renamed.pcc>=curr)]
    all_nodes = np.unique([*all_edges.gene1, *all_edges.gene2])
    
    print(prev, sum(sensor_names.isin(all_nodes)), sum(sensor_names.isin(pcc_dfs[curr].orf_name)) , len(sensor_names))


# %%
for i in thr_list:
    sensor_names = sensor_dfs[i].orf_name
    all_edges = df_long_renamed.loc[df_long_renamed.pcc>=i]
    g = nx.from_pandas_edgelist(all_edges,source='gene1',target='gene2',edge_attr='pcc')
    fig, ax = plt.subplots(figsize=(7,7))
    data = [val for u,v,val in g.edges(sensor_names.tolist(), data='pcc')]
    ax.hist(data)
    ax.set_title(i)
    print(i , np.median(data), len(g.nodes))

# %%
exn = pd.read_csv('../../pgsNetwork/analysis/data/derived_data/SGA_data_combined.csv','\t')

# %%
aa = exn.loc[exn['Query allele name'].isin(effector_dfs[0.2].orf_name) & exn['Array allele name'].isin(effector_dfs[0.2].orf_name) ]

# %%
exe_names = np.unique([*aa.loc[aa.data_source=='ExE']['Query allele name'],*aa.loc[aa.data_source=='ExE']['Array allele name']])
nxe_names = np.unique(aa.loc[aa.data_source=='ExN_NxE']['Array allele name'])

# %%
strain_ids = pd.read_csv('../../pgsNetwork/trymake/pgsNetwork//analysis/data/derived_data/strain_ids_with_experiment_count_all.csv')

# %%
strain_ids

# %%
ii = 0.2
effector_dfs[ii].loc[effector_dfs[ii].orf_name.isin(strain_ids.loc[strain_ids.maincat=='essential']['Allele Gene name'])].shape[0]/effector_dfs[ii].shape[0]

# %%
