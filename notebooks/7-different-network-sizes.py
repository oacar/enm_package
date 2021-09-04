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
thr_list =  [0.05,0.1,0.2,0.25,0.3,0.35,0.4,0.45,0.5]

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
        
    fig, axs = plt.subplots(1,3,figsize=(15,5))
    axs[0].plot(thr_list, n_clusters, 'o-')
    axs[0].set_ylabel(f'Number of {yaxis} clusters')
    axs[0].set_xlabel(xlab)
    
    axs[1].plot(thr_list, n_go_clusters, 'o-')
    axs[1].set_ylabel(f'Number of go enriched {yaxis} clusters')
    axs[1].set_xlabel(xlab)
    
    axs[2].plot(thr_list, [i/j if i!=0 else 0 for i,j in zip(n_go_clusters,n_clusters) ], 'o-')
    axs[2].set_ylabel(f'% of go enriched {yaxis} clusters')
    axs[2].set_xlabel(xlab)


# %%
plot_go_thr_comparison(effector_dfs,'effector_cluster', 'effector')

# %%
plot_go_thr_comparison(sensor_dfs,'sensor_cluster', 'sensor')

# %%

sensor_overlap = pd.DataFrame({thr : [len(np.intersect1d(sensor_dfs[thr].orf_name, sensor_dfs[i].orf_name)) for i in thr_list] for thr in sorted(thr_list)})
sensor_overlap.index = thr_list

# %% [markdown]
# # Sensor overlap between thresholds

# %%
sns.heatmap(sensor_overlap / np.diag(sensor_overlap))

# %% [markdown]
# # Effector overlap between thresholds

# %%

effector_overlap = pd.DataFrame({thr : [len(np.intersect1d(effector_dfs[thr].orf_name, effector_dfs[i].orf_name)) for i in thr_list] for thr in sorted(thr_list)})
effector_overlap.index = thr_list
sns.heatmap(effector_overlap / np.diag(effector_overlap))

# %%
import itertools as itr


# %%
def get_pairwise_intersect(dfs, thr_list):
    pairs = itr.combinations(thr_list,2)
    nt = lambda a, b: [len(np.intersect1d(dfs[a].orf_name,effector_dfs[b].orf_name)) , len(dfs[b].orf_name)]
    res = dict([ (t, nt(*t)) for t in pairs ])

    return res


# %% jupyter={"outputs_hidden": true}
for i in range(1,len(thr_list)):
    prev = thr_list[i-1]
    curr = thr_list[i]
    print(prev,curr)
    curr_sensors = sensor_dfs[curr].orf_name.tolist()
    pcc_dfs[prev]['sensor_on_larger_thr'] = pd.Categorical([1 if i else 0 for i in pcc_dfs[prev].orf_name.isin(curr_sensors)])
    
    fig, ax = plt.subplots(figsize=(6,6))
    pcc_dfs[prev].loc[pcc_dfs[prev].sensor_on_larger_thr==1].plot.scatter('deg','sens',c='sensor_on_larger_thr',cmap='viridis',ax=ax)
    ax.set_title(prev)
    ax.set_ylim(pcc_dfs[prev].sens.min(),pcc_dfs[prev].sens.max())

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
# # Redefine sensors/effectors

# %%
# %load_ext rpy2.ipython

# %%
thr = 0.2
dff = pcc_dfs[thr].drop('sensor_on_larger_thr',axis=1)


# %%
#Import the SignatureTranslatedAnonymousPackage 
from rpy2.robjects.packages import STAP
import rpy2.robjects as ro
#Read the file with the R code snippet
with open('../R/compare_clustering.R', 'r') as f:
    string = f.read()
#Parse using STAP
r_func = STAP(string, "compare_clustering")

# %% language="R"
#
# library(Ckmeans.1d.dp)
# library(tidyverse)
# library(glue)
# compare_clustering <- function(filename, number_of_clusters, column_name){
#     df <- readr::read_csv(filename)
#     res <- Ckmeans.1d.dp(pull(df[column_name]),number_of_clusters)
#   #  plotBIC(res)
#     df$ck = res$cluster
#   #  df$color = df$ck#>quantile(unique(df$ck),quantile_threshold)
# #     ggplot(df, aes_string(x='deg',y=column_name,color="color"))+
# #         geom_point()
#     return(df)
# }

# %%
from enm.utils import *
#goea, geneid2name = create_goea()

# %%
def get_go_enrichment_by_cluster(clustered_df, cluster_column_name, goea, geneid2name):
    go_res = {}
    for i in sorted(clustered_df[cluster_column_name].unique()):
        orf_list = clustered_df.loc[clustered_df[cluster_column_name]==i]
        go_df = query_goatools(orf_list, goea, geneid2name)
        go_res[i] = go_df
    return go_res


# %%
def go_for_a_threshold(thr, column_name,goea, geneid2name):
    #%%R -i dff -i thr -o clustered_df
    clustered_df= convert_rpy2_to_pandas(r_func.compare_clustering(f'../data/interim_{thr}/pcc_df.csv',ro.IntVector((2,50)), column_name))
    go_res = get_go_enrichment_by_cluster(clustered_df, 'ck', goea, geneid2name)
    return clustered_df, go_res


# %% jupyter={"outputs_hidden": true}
clustered_dfs_sens = {}
clustered_dfs_eff = {}

go_res_all_sens = {}
go_res_all_eff = {}

for i in thr_list:
    clustered_df_sens, go_res_sensor = go_for_a_threshold(i, 'sens',goea, geneid2name)
    clustered_df_eff, go_res_eff = go_for_a_threshold(i, 'eff',goea, geneid2name)
    clustered_dfs_sens[i] = clustered_df_sens
    clustered_dfs_eff[i] = clustered_df_eff
    go_res_all_sens[i] = go_res_sensor
    go_res_all_eff[i] = go_res_eff


# %% jupyter={"outputs_hidden": true}
clustered_df_sens, go_res_sensor_02 = go_for_a_threshold(thr, 'sens',goea, geneid2name)
clustered_df_eff, go_res_eff_02 = go_for_a_threshold(thr, 'eff',goea, geneid2name)

# %%
go_res_eff_02[10]

# %% [raw]
# def plot_kmeans(data,x,y):
#     kmeans = KMeans(init="k-means++", n_clusters=5, n_init=1,
#                 random_state=0)
#     X = np.array(data[y]).reshape(-1,1)
#     y_km = kmeans.fit_predict(X)
#     data['kmeans'] = y_km
#     data.groupby('kmeans').sum().loc[y]
#     sns.scatterplot(data=data,x=x,y=y,hue=y_km)

# %%
data =pcc_dfs[0.2].copy()
data['kmeans'] = y_km


# %%
data.groupby('kmeans').sum()['eff'].argmax()

# %%
idx =0.45
plot_kmeans(pcc_dfs[idx], 'deg','sens')

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
