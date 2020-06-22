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
from mlxtend.evaluate import permutation_test 

# Random collectivity
data_path = 'data/interim/pcc_0525/'
with open(f'{data_path}/rewire_data_er.pickle','rb') as f:
    e_er = pickle.load(f)
with open(f'{data_path}/pcc.pickle','rb') as f:
    e_pcc = pickle.load(f)
with open(f'{data_path}/rewire_data_nodegseq.pickle','rb') as f:
    e_nodegseq = pickle.load(f)
figure_path = e_er[0].figure_path = e_pcc.figure_path = e_nodegseq[0].figure_path = 'reports/figures/pcc_0601/'
e_pcc.e_list[0].figure_path=figure_path
e_er[0].plot_collectivity(figure_name = 'coll_er', figure_extension='png')
e_pcc.plot_collectivity(figure_name = 'coll_pcc', figure_extension='png')
e_nodegseq[0].plot_collectivity(figure_name = 'coll_nodegseq', figure_extension='png')
e_pcc.e_list[0].plot_collectivity(figure_name = 'coll_keepingdegseq', figure_extension='png')

plot_lambda_collectivity(e_pcc.gnm.getEigvals(),e_pcc.coll,figure_path,figure_name ='lambda_coll_real')
plot_lambda_collectivity(e_er[0].gnm.getEigvals(),e_er[0].coll,figure_path,figure_name ='lambda_coll_er')
plot_lambda_collectivity(e_pcc.e_list[0].gnm.getEigvals(),e_pcc.e_list[0].coll,figure_path,figure_name ='lambda_coll_keepingdegseq')
#General plots
e_pcc.plot_network_spring()
e_pcc.plot_vector(sorted=True)
e_pcc.plot_scatter(x='deg',y='eff',figure_name='deg_eff',figure_extension='pdf')
e_pcc.plot_scatter(x='deg',y='sens',figure_name='deg_sens',figure_extension='pdf')

# correlation_density
er_df = pd.read_csv(f'{data_path}/rewire_data_er.csv')
nodegseq_df = pd.read_csv(f'{data_path}/rewire_data_nodegseq.csv')
#e_pcc.plot_correlation_density(x='eff',y='deg',figure_extension='png')
plot_correlation_density(e_pcc.df, [e_pcc.rewire_df,nodegseq_df,er_df], x='eff', y='deg', figure_path=figure_path, figure_extension='png')
plot_correlation_density(e_pcc.df, [e_pcc.rewire_df,nodegseq_df,er_df], x='sens', y='deg', figure_path=figure_path, figure_extension='png',correlation='spearman')

#most collective eigenvector plots
e_er[0].plot_vector(sorted=True,figure_name='eig_0_er',figure_extension='png')
e_nodegseq[0].plot_vector(sorted=True,figure_name='eig_0_nodegseq',figure_extension='png')
e_pcc.plot_vector(sorted=True,figure_name='eig_0_pcc',figure_extension='png')
e_pcc.e_list[0].plot_vector(sorted=True,figure_name='eig_0_keepingdegseq',figure_extension='png')


#smalles lambda eigenvector plots
plot_vector(e_er[0].gnm.getEigvecs()[:,0],figure_path=figure_path,sorted=True,figure_name='eig_lambda_er',figure_extension='png',color_id=-1)
plot_vector(e_nodegseq[0].gnm.getEigvecs()[:,0],figure_path=figure_path,sorted=True,figure_name='eig_lambda_nodegseq',figure_extension='png',color_id=-1)
plot_vector(e_pcc.gnm.getEigvecs()[:,0],figure_path=figure_path,sorted=True,figure_name='eig_lambda_pcc',figure_extension='png',color_id=-1)
plot_vector(e_pcc.e_list[0].gnm.getEigvecs()[:,0],figure_path=figure_path,sorted=True,figure_name='eig_lambda_keepingdegseq',figure_extension='png',color_id=-1)

#GO analysis with different collectivity

go_results_folder="data/interim/go_results_excl/"
fls  = os.listdir(go_results_folder)
fls.sort(key=lambda f: int(re.sub('\D', '', f)))
go_res = [pd.read_csv(f"data/interim/go_results_excl/{x}",'\t') for x in fls]
fls_numbers = sorted([int(re.sub('\D','',f)) for f  in fls])
go_res_dd = dict(zip([e_pcc.coll_index_sorted[i] for i in fls_numbers],[len(np.unique(list(itr.chain(*[i.split(', ') for i in gr.study_items.values if pd.isna(i)==False])))) for gr in go_res]))#[e_pcc.coll[i] for i in e_pcc.coll_index_sorted]
for i in range(len(e_pcc.graph_gc.nodes)-1):
    if i not in go_res_dd.keys():
        go_res_dd[i]=0

go_res_df = pd.DataFrame(go_res_dd.items(),columns=['id','number_of_enriched_genes'])
go_res_df = go_res_df.sort_values(by='id')
go_res_df['coll']=e_pcc.coll#[e_pcc.coll[i] for i in e_pcc.coll_index_sorted]
go_res_df['lambda']=e_pcc.gnm.getEigvals()#[i] for i in e_pcc.coll_index_sorted]
#go_res_df#.iloc[8:9,:]

fig, ax = plt.subplots(figsize=(6,6))
ax.scatter((go_res_df['coll']),go_res_df.number_of_enriched_genes, s=10)
ax.set_xlabel('Collectivity')
ax.set_ylabel('Number of GO enriched genes')
#ax.set_xlim(0.1,-0.01)
fig.savefig(f'{figure_path}/coll_goterm.png')

fig, ax = plt.subplots(figsize=(6,6))
ax.scatter((go_res_df['lambda']),go_res_df.number_of_enriched_genes, s=10)#, c=go_res_df.coll,cmap='gray')
ax.set_xlabel('lambda')
ax.set_ylabel('Number of GO enriched genes')
#ax.set_xlim(0.1,-0.01)
fig.savefig(f'{figure_path}/lambda_goterm.png')

fig, ax = plt.subplots(figsize=(6,6))
ax.scatter((go_res_df['id']),go_res_df.number_of_enriched_genes, s=10)
ax.set_xlabel('Modes')
ax.set_ylabel('Number of GO enriched genes')
#ax.set_xlim(0.1,-0.01)
fig.savefig(f'{figure_path}/id_goterm.png')

# Modes with GO enrichment and their collectivity/lambda relationship

go_res_df['go_bool'] = np.where(go_res_df.number_of_enriched_genes>0,True,False)
c1 = go_res_df[go_res_df['go_bool']]
c2 = go_res_df[go_res_df['go_bool']==False]
p_value = permutation_test(c1['lambda'],c2['lambda'],method='approximate',num_rounds=10000,seed=1)
fig, ax = plt.subplots(figsize=(6,6))
#ax.scatter((go_res_df['id']),go_res_df.number_of_enriched_genes, s=10)
ax.set_xlabel('Collectivity')
ax.set_ylabel('Density')
sns.distplot(c1['coll'],ax=ax,color='skyblue',hist=False,label='with GO enrichment')
sns.distplot(c2['coll'],ax=ax,color='red',hist=False,label='no GO enrichment')
#ax.set_xlim(0.1,-0.01)
ax.legend()
fig.savefig(f'{figure_path}/coll_go_hist.png')

p_value = permutation_test(c1['lambda'],c2['lambda'],method='approximate',num_rounds=10000,seed=1)
fig, ax = plt.subplots(figsize=(6,6))
#ax.scatter((go_res_df['id']),go_res_df.number_of_enriched_genes, s=10)
ax.set_xlabel('Lambda')
ax.set_ylabel('Density')
sns.distplot(c1['lambda'],ax=ax,color='skyblue',hist=False,label='with GO enrichment')
sns.distplot(c2['lambda'],ax=ax,color='red',hist=False,label='no GO enrichment')
#ax.set_xlim(0.1,-0.01)
ax.legend()
fig.savefig(f'{figure_path}/lambda_go_hist.png')

p_value = permutation_test(c1['lambda'],c2['lambda'],method='approximate',num_rounds=10000,seed=1)
fig, ax = plt.subplots(figsize=(6,6))
sns.boxplot(x='go_bool',y='lambda',data=go_res_df)

# statistical annotation
x1, x2 = 0, 1   # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
y, h, col = go_res_df['lambda'].max() + 2, 2, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h, f"p = {p_value}", ha='center', va='bottom', color=col)

ax.set_ylabel('Lambda')
ax.set_xlabel('GO enriched')
fig.savefig(f'{figure_path}/lambda_go_boxplot.pdf')

p_value = permutation_test(c1['coll'],c2['coll'],method='approximate',num_rounds=10000,seed=1)
fig, ax = plt.subplots(figsize=(6,6))
sns.boxplot(x='go_bool',y='coll',data=go_res_df)

# statistical annotation
x1, x2 = 0, 1   # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
y, h, col = go_res_df['coll'].max() + 0.005, 0.005, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h, f"p = {p_value}", ha='center', va='bottom', color=col)

ax.set_ylabel('Collectivity')
ax.set_xlabel('GO enriched')
#ax.set_xlim(0.1,-0.01)
#ax.legend()
fig.savefig(f'{figure_path}/collectivity_go_boxplot.pdf')


# Network map with GO
go_df_list =[pd.read_csv(f"{go_results_folder}/{i}.tsv",'\t') for i in range(7)]
strain_ids = pd.read_csv('data/interim/strain_ids_with_experiment_count_all.csv')
orf_name_mapping = dict(zip(strain_ids.iloc[:,1],strain_ids.iloc[:,0]))
#combined_df = pd.merge(e.df, strain_ids, left_on='orf_name',right_on='Allele Gene name')
nx.set_node_attributes(e_pcc.graph_gc,orf_name_mapping, 'orf_name')
e_pcc.plot_network_spring(plot_go=True,go_df_list=go_df_list,plot_legend=True,level_list=[0.1]*7,figure_extension='pdf')

#[go_res_dd[i] for i in range(20)]

