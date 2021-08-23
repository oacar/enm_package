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
thr_list =  [0.05,0.1,0.2,0.25,0.3,0.35,0.4,0.45,0.5]

# %%
for i in np.arange(0.2,0.75,0.05):
    df = cpcc.loc[cpcc.pcc>i]
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

# %%
sensor_df_names = [ ]


# %%
def get_result_dfs(fname, thr_list, default_thr):
    dfs = {thr: pd.read_csv(f"../data/interim_{thr}/{fname}.csv") for thr in thr_list if thr!=default_thr}
    dfs[default_thr] = pd.read_csv(f"../data/interim/{fname}.csv")
    if default_thr not in thr_list:
        thr_list.insert(0, default_thr)                               
    return dfs


# %%
effector_dfs = get_result_dfs('effectors_df', thr_list, 0.2)

# %%
sensor_dfs =  get_result_dfs('sensors_df', thr_list, 0.2)

# %%
sensor_dfs[0.2].dropna(subset=['sensor_cluster']).loc[:,'sensor_cluster'].nunique()


# %%
def plot_go_thr_comparison(dfs, col, yaxis):
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
        
    fig, axs = plt.subplots(2,2,figsize=(10,10))
    axs[0][0].plot(thr_list, n_goterms,'o-')
    axs[0][0].set_ylabel(f'Number of GO enriched {yaxis}s')

    axs[1][0].plot(thr_list, rat_goterms,'o-')
    axs[1][0].set_ylabel(f'%of {yaxis}s with GO term enrichments in all {yaxis}s')


    axs[0][1].plot(thr_list, n_clusters, 'o-')
    axs[0][1].set_ylabel(f'Number of {yaxis} clusters')
    axs[0][1].set_xlabel(xlab)
    
    axs[1][1].plot(thr_list, n_go_clusters, 'o-')
    axs[1][1].set_ylabel(f'Number of go enriched {yaxis} clusters')
    axs[1][1].set_xlabel(xlab)

# %%
plot_go_thr_comparison(effector_dfs,'effector_cluster', 'effector')

# %%
plot_go_thr_comparison(sensor_dfs,'sensor_cluster', 'sensor')

# %%
[np.intersect1d(sensor_dfs[0.5].orf_name, sensor_dfs[i].orf_name) for i in thr_list]

# %%
[np.intersect1d(effector_dfs[0.25].orf_name, effector_dfs[i].orf_name) for i in thr_list]

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

# %%
[pcc_dfs[i].plot.scatter('deg','eff') for i in thr_list]

# %% jupyter={"outputs_hidden": true}
sensor_dfs[0.3]

# %%
# %config Completer.use_jedi = False

# %%
#for i in thr_list:
i=0.5
for j in thr_list:
    print(i, j , (np.intersect1d(sensor_dfs[i].orf_name, pcc_dfs[j].orf_name)))

# %%
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import seaborn as sns



# %%
def plot_kmeans(data,x,y):
    kmeans = KMeans(init="k-means++", n_clusters=5, n_init=1,
                random_state=0)
    X = np.array(data[y]).reshape(-1,1)
    y_km = kmeans.fit_predict(X)
    data['kmeans'] = y_km
    data.groupby('kmeans').sum().loc[y]
    sns.scatterplot(data=data,x=x,y=y,hue=y_km)


# %%
data =pcc_dfs[0.2].copy()
data['kmeans'] = y_km


# %%
data.groupby('kmeans').sum()['eff'].argmax()

# %%
idx =0.4
plot_kmeans(pcc_dfs[idx], 'deg','eff')

# %%
