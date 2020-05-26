import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import pickle
import os
import re
import itertools as itr
from src.visualize.visualize import plot_correlation_density, plot_vector
from mlxtend.evaluate import permutation_test 

# Random collectivity
data_path = 'data/interim/pcc_0525/'
with open(f'{data_path}/rewire_data_er.pickle','rb') as f:
    e_er = pickle.load(f)
with open(f'{data_path}/pcc.pickle','rb') as f:
    e_pcc = pickle.load(f)
with open(f'{data_path}/rewire_data_nodegseq.pickle','rb') as f:
    e_nodegseq = pickle.load(f)
figure_path = e_er[0].figure_path = e_pcc.figure_path = e_nodegseq[0].figure_path = 'reports/figures/paper_figures_try/'
e_pcc.e_list[0].figure_path=figure_path
e_er[0].plot_collectivity(figure_name = 'er_coll', figure_extension='pdf')
e_pcc.plot_collectivity(figure_name = 'pcc_coll', figure_extension='pdf')
e_nodegseq[0].plot_collectivity(figure_name = 'nodegseq_coll', figure_extension='pdf')
e_pcc.e_list[0].plot_collectivity(figure_name = 'keepingdegseq_coll', figure_extension='pdf')

# correlation_density
er_df = pd.read_csv(f'{data_path}/rewire_data_er.csv')
nodegseq_df = pd.read_csv(f'{data_path}/rewire_data_nodegseq.csv')
#e_pcc.plot_correlation_density(x='eff',y='deg',figure_extension='pdf')
plot_correlation_density(e_pcc.df, [e_pcc.rewire_df,nodegseq_df,er_df], x='eff', y='deg', figure_path=figure_path, figure_extension='pdf')
plot_correlation_density(e_pcc.df, [e_pcc.rewire_df,nodegseq_df,er_df], x='sens', y='deg', figure_path=figure_path, figure_extension='pdf',correlation='spearman')

#most collective eigenvector plots
e_er[0].plot_vector(sorted=True,figure_name='eig_0_er',figure_extension='pdf')
e_nodegseq[0].plot_vector(sorted=True,figure_name='eig_0_nodegseq',figure_extension='pdf')
e_pcc.plot_vector(sorted=True,figure_name='eig_0_pcc',figure_extension='pdf')
e_pcc.e_list[0].plot_vector(sorted=True,figure_name='eig_0_keepingdegseq',figure_extension='pdf')


#smalles lambda eigenvector plots
plot_vector(e_er[0].gnm.getEigvecs()[:,0],figure_path=figure_path,sorted=True,figure_name='eig_lambda_er',figure_extension='pdf')
plot_vector(e_nodegseq[0].gnm.getEigvecs()[:,0],figure_path=figure_path,sorted=True,figure_name='eig_lambda_nodegseq',figure_extension='pdf')
plot_vector(e_pcc.gnm.getEigvecs()[:,0],figure_path=figure_path,sorted=True,figure_name='eig_lambda_pcc',figure_extension='pdf')
plot_vector(e_pcc.e_list[0].gnm.getEigvecs()[:,0],figure_path=figure_path,sorted=True,figure_name='eig_lambda_keepingdegseq',figure_extension='pdf')

#GO analysis with different collectivity

fls  = os.listdir('data/interim/go_results_excl/')
fls.sort(key=lambda f: int(re.sub('\D', '', f)))
go_res = [pd.read_csv(f"data/interim/go_results_excl/{x}",'\t') for x in fls]
go_res_dd = dict(zip([int(re.sub('\D','',f)) for f  in fls],[len(np.unique(list(itr.chain(*[i.split(', ') for i in gr.study_items.values if pd.isna(i)==False])))) for gr in go_res]))#[e_pcc.coll[i] for i in e_pcc.coll_index_sorted]
for i in range(len(e_pcc.graph_gc.nodes)-1):
    if i not in go_res_dd.keys():
        go_res_dd[i]=0

go_res_df = pd.DataFrame(go_res_dd.items(),columns=['id','number_of_enriched_genes'])
go_res_df = go_res_df.sort_values(by='id')
go_res_df['coll']=[e_pcc.coll[i] for i in e_pcc.coll_index_sorted]
go_res_df['lambda']=[e_pcc.gnm.getEigvals()[i] for i in e_pcc.coll_index_sorted]
#go_res_df#.iloc[8:9,:]

fig, ax = plt.subplots(figsize=(6,6))
ax.scatter((go_res_df['coll']),go_res_df.number_of_enriched_genes, s=10)
ax.set_xlabel('Collectivity')
ax.set_ylabel('Number of GO enriched genes')
#ax.set_xlim(0.1,-0.01)
fig.savefig(f'{figure_path}/coll_goterm.pdf')

fig, ax = plt.subplots(figsize=(6,6))
ax.scatter((go_res_df['lambda']),go_res_df.number_of_enriched_genes, s=10)#, c=go_res_df.coll,cmap='gray')
ax.set_xlabel('lambda')
ax.set_ylabel('Number of GO enriched genes')
#ax.set_xlim(0.1,-0.01)
fig.savefig(f'{figure_path}/lambda_goterm.pdf')

fig, ax = plt.subplots(figsize=(6,6))
ax.scatter((go_res_df['id']),go_res_df.number_of_enriched_genes, s=10)
ax.set_xlabel('Modes')
ax.set_ylabel('Number of GO enriched genes')
#ax.set_xlim(0.1,-0.01)
fig.savefig(f'{figure_path}/id_goterm.pdf')

# Modes with GO enrichment and their collectivity/lambda relationship

go_res_df['go_bool'] = np.where(go_res_df.number_of_enriched_genes>0,True,False)
c1 = go_res_df[go_res_df['go_bool']]
c2 = go_res_df[go_res_df['go_bool']==False]
p_value = permutation_test(c1['coll'],c2['coll'],method='approximate',num_rounds=10000,seed=1)
fig, ax = plt.subplots(figsize=(6,6))
#ax.scatter((go_res_df['id']),go_res_df.number_of_enriched_genes, s=10)
ax.set_xlabel('Collectivity')
ax.set_ylabel('Density')
sns.distplot(c1['coll'],ax=ax,color='skyblue',hist=False,label='with GO enrichment')
sns.distplot(c2['coll'],ax=ax,color='red',hist=False,label='no GO enrichment')
#ax.set_xlim(0.1,-0.01)
ax.legend()
fig.savefig(f'{figure_path}/coll_go_hist.pdf')

p_value = permutation_test(c1['lambda'],c2['lambda'],method='approximate',num_rounds=10000,seed=1)
fig, ax = plt.subplots(figsize=(6,6))
#ax.scatter((go_res_df['id']),go_res_df.number_of_enriched_genes, s=10)
ax.set_xlabel('Lambda')
ax.set_ylabel('Density')
sns.distplot(c1['lambda'],ax=ax,color='skyblue',hist=False,label='with GO enrichment')
sns.distplot(c2['lambda'],ax=ax,color='red',hist=False,label='no GO enrichment')
#ax.set_xlim(0.1,-0.01)
ax.legend()
fig.savefig(f'{figure_path}/lambda_go_hist.pdf')


