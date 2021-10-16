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

# %%
# %config Completer.use_jedi = False

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
random.seed(4812)        # or any integer
np.random.seed(4813)

# %%
# with open('../data/interim/pcc.pickle', 'rb') as f:
#     enm_ee = pickle.load(f)

# %%
strain_ids =pd.read_csv('../data/interim/strain_ids.csv')


# %%
e = Enm('rew')

# %%
e.read_network('../data/interim/costanzo_pcc_ALL',sep=',')

# %%
gc_rew = rewire_network(e.graph_gc)
e_rew = Enm('rewired')
e_rew.G=gc_rew
e_rew.giant_component()
e_rew.gnm_analysis()
e_rew.df = pd.merge(e_rew.df,strain_ids, left_on='orf_name',right_on='gene1')
e_rew.get_sensor_effector(True)

# %%
goea, geneid2name,a = create_goea(gaf='../data/raw/ontology/sgd.gaf', obo='../data/raw/ontology/go-basic.obo',background='../data/interim/go_background_list',
                               sdg_info_tab='../data/raw/ontology/SGD_features.tab')

# %% tags=[] jupyter={"outputs_hidden": true}
e_rew.analyze_components_biology(goea, geneid2name,True)
e_rew.analyze_components_biology(goea, geneid2name,False)

# %%
e_rew.get_sensor_effector()

# %%
sensors_pcc = e_rew.sensors_df

# %%
pos = e_rew.graph_gc.nodes('pos')

# %%
from matplotlib.lines import Line2D

fig, ax = plt.subplots(figsize=(9,9))
#axs = ax.ravel()
legend_elements = [    ]

#for i in range(len(sensor_order)):
e_rew.plot_network_spring(ax=ax,
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
nx.draw_networkx_nodes(nx.induced_subgraph(e_rew.graph_gc, sensors_pcc.orf_name.tolist()),
                       pos=pos, 
                       node_color='black', alpha=1, node_shape='^')

# for itr, i in enumerate(sensor_order):
#     #print(i, effector_colors[itr])

#     orf_names_to_plot = sensors_pcc.loc[sensors_pcc.go_group==i, 'orf_name'].tolist()
#     nx.draw_networkx_nodes(e_rew.graph_gc, nodelist=orf_names_to_plot, node_size=200, pos=pos,
#                           node_color=sensor_colors[itr],
#                           node_shape='^',edgecolors='black',
#                           linewidths=1)
#     legend_elements.append(
#         Line2D([0], [0], marker='^', color='black', label=f'{i}',
#                               markerfacecolor=sensor_colors[itr], markersize=30, linestyle="None")
#     )
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
lgd = ax.legend(handles=legend_elements, fontsize=22,loc='center left', bbox_to_anchor=(-0.5, -0.1),ncol=5)
nx.draw_networkx_edges(nx.induced_subgraph(e_rew.graph_gc, sensors_pcc.orf_name.tolist()),pos=pos, edge_color='red', alpha=0.5)
plt.savefig(f'../reports/figures/paper_figures_supp/figs1.png',bbox_inches='tight',dpi=150)

# %% [markdown]
# # Calculate number of clusters in rewired

# %% tags=[]
with open('../data/supp/pcc.pickle','rb') as f:
    e_pcc = pickle.load(f)

# %% tags=[] jupyter={"outputs_hidden": true}
e_pcc.simulate_rewire(sim_num = 500)

# %% tags=[]
for i in e_pcc.e_list:
    i.df = pd.merge(i.df , strain_ids , left_on = 'orf_name', right_on='gene1')
    i.get_sensor_effector()
    i.analyze_components_biology(goea, geneid2name)


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


# %% tags=[] jupyter={"outputs_hidden": true}
[save_rewired_data(e_pcc,i,'../data/interim/rewired_data500') for i in range(500)]

# %%
import glob

sensors_fls = glob.glob('../data/interim/rewired_data500/sensors_df_rew*')

# %%
sensors_fls

# %%
sensors_df_rew = pd.concat([pd.read_csv(f'../data/interim/rewired_data500/{idx}/sensors_df_rew_{idx}.csv') for idx in range(500)],keys=range(500)).reset_index(level=0)

# %%
res = sensors_df_rew.groupby('level_0')['sensor_cluster'].nunique()
res2 = sensors_df_rew.groupby('level_0')['go_group'].nunique()

# %%
#res = [i.sensors_df.dropna(subset=['sensor_cluster']).sensor_cluster.nunique() for i in e_pcc.e_list]
#res2 = [i.sensors_df.dropna(subset=['go_group']).go_group.nunique() for i in e_pcc.e_list]

# %%
idxx = np.argwhere(np.array(res)>0).reshape(1,-1)[0]
idxx2 = np.argwhere(np.array(res2)>0).reshape(1,-1)[0]

# %%
len(idxx)

# %%
len(idxx)/len(res)

# %%
len(idxx2)/len(res2)

# %%
len(idxx2)/len(idxx)


# %%
def plot_sensors(e_pcc,idx):
    e =e_pcc.e_list[idx]
    sub_orfs =  e.sensors_df.dropna(subset=['sensor_cluster']).orf_name.tolist()
    g = e.graph_gc
    induced_g = nx.induced_subgraph(g,sub_orfs)
    sub_nw = get_subnetwork(g, sub_orfs, radius= 1)
    pos_sub = nx.spring_layout(sub_nw)
    fig, ax_ = plt.subplots() 
    nx.draw_networkx_nodes(sub_nw,ax=ax_, pos = pos_sub, node_color = ['none' if i in sub_orfs else 'k' for i in sub_nw.nodes])
    nx.draw_networkx_nodes(nx.induced_subgraph(sub_nw, sub_orfs), pos=pos_sub, node_shape='^', node_color='black')
    nx.draw_networkx_edges(sub_nw,ax=ax_, pos = pos_sub)
    nx.draw_networkx_edges(nx.induced_subgraph(sub_nw, sub_orfs), pos=pos_sub, edge_color='red')
    plt.show()
    #nx.draw_networkx(sub_nw)


# %% tags=[] jupyter={"outputs_hidden": true}
for i in idxx:
    print(i)
    plot_sensors(e_pcc,i)

# %%
fig, ax =plt.subplots(figsize=(5,5))
ax.hist(res,3,color='navajowhite')
ax.set_xlabel('Number of sensor clusters')
ax.set_ylabel('Count')
ax.set_title('Number of sensor clusters\nat 500 rewired networks\ncompared to real network')
ax.axvline(9,c='darkblue',linestyle='-.')
plt.legend(handles = [
    Line2D([0],[0],color='navajowhite',linewidth=10,label='Rewired'),
    Line2D([0],[0],color='darkblue',linestyle='-.', label='Real')
], loc='upper center')

fig.savefig('../reports/figures/paper_figures_supp/rewired_sensor_count.png', bbox_inches='tight',dpi=150)


# %%
fig, ax =plt.subplots(figsize=(5,5))
ax.hist(sensors_df_rew.groupby('level_0')['go_group'].nunique(),3,color='navajowhite')
ax.set_xlabel('Number of GO enriched sensor clusters')
ax.set_ylabel('Count')
ax.set_title('Number of GO enriched\nsensor clusters\nat 500 rewired networks\ncompared to real network')
ax.axvline(5,c='darkblue',linestyle='-.')
plt.legend(handles = [
    Line2D([0],[0],color='navajowhite',linewidth=10,label='Rewired'),
    Line2D([0],[0],color='darkblue',linestyle='-.', label='Real')
], loc='upper center')

fig.savefig('../reports/figures/paper_figures_supp/rewired_go_sensor_count.png', bbox_inches='tight',dpi=150)


# %%
res

# %%
fig, ax =plt.subplots(figsize=(5,5))
ax.hist(res,3,color='navajowhite')
ax.set_xlabel('Number of sensor clusters')
ax.set_ylabel('Count')
ax.set_title('Number of sensor clusters\nat 100 rewired networks\ncompared to real network')
ax.axvline(9,c='darkblue',linestyle='-.')
plt.legend(handles = [
    Line2D([0],[0],color='navajowhite',linewidth=10,label='Rewired'),
    Line2D([0],[0],color='darkblue',linestyle='-.', label='Real')
], loc='upper center')

fig.savefig('../reports/figures/paper_figures_supp/rewired_sensor_count.png', bbox_inches='tight',dpi=150)


# %%
