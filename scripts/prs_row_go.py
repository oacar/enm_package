from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx
import numpy as np
import pandas as pd
import pickle
from enm.visualize import plot_correlation_density, plot_vector, plot_lambda_collectivity
from enm.Enm import Enm
from enm.utils import *


with open(snakemake.input.pickle_file_name,'rb') as f:
    e_pcc = pickle.load(f)

e_pcc.get_prs(no_diag=False)
e_pcc.get_rwr_mat()


# In[ ]:

rwr_ranked_goa_results = []
ranked_goa_results = []

goea, geneid2name = create_goea(gaf = snakemake.input.gaf, 
                                obo_fname=snakemake.input.obo, 
                                background=snakemake.input.background_file, 
                                sgd_info_tab = snakemake.input.sgd_info)

ranked_goa_results_rows = []
rwr_ranked_goa_results_rows = []
if 'prs_row' in snakemake.output.keys():
    for i in tqdm(e_pcc.nodes):
        goa_df_prs_row = query_goatools(e_pcc.df.loc[e_pcc.df.orf_name.isin(e_pcc.prs_mat_df.loc[i,:].sort_values(ascending=False)[:50].index)],
                goea,geneid2name)
        ranked_goa_results_rows.append(goa_df_prs_row)
    combined_df = combine_data(ranked_goa_results_rows)
    combined_df.to_csv(snakemake.output.prs_row,index=False)
    #with open(snakemake.output.prs_row, 'wb') as f:
    #    pickle.dump(ranked_goa_results_rows,f)
elif 'rwr_row' in snakemake.output.keys():
    for i in tqdm(e_pcc.nodes):
        goa_df_rwr_row = query_goatools(e_pcc.df.loc[e_pcc.df.orf_name.isin(e_pcc.rwr_mat_df.loc[i,:].sort_values(ascending=False)[:50].index)],
                goea,geneid2name)
        rwr_ranked_goa_results_rows.append(goa_df_rwr_row)
    combined_df = combine_data(rwr_ranked_goa_results_rows)
    combined_df.to_csv(snakemake.output.rwr_row,index=False)

# In[ ]:




# In[ ]:


#goea, geneid2name = create_goea(gaf = 'data/raw/ontology/sgd.gaf', obo_fname='data/raw/ontology/go-basic.obo', background='data/interim/costanzo_gc_bg.tsv', sgd_info_tab = 'data/raw/ontology/SGD_features.tab')

#for i in tqdm(e_pcc.nodes):
if 'prs_column' in snakemake.output.keys():
    for i in tqdm(e_pcc.nodes):
        goa_df_prs_col = query_goatools(e_pcc.df.loc[e_pcc.df.orf_name.isin(e_pcc.prs_mat_df.loc[:,i].sort_values(ascending=False)[:50].index)],
                goea,geneid2name)
        ranked_goa_results.append(goa_df_prs_col)

    combined_df = combine_data(ranked_goa_results)
    combined_df.to_csv(snakemake.output.prs_column,index=False)
if 'rwr_column' in snakemake.output.keys():
    for i in tqdm(e_pcc.nodes):
        goa_df_rwr_col = query_goatools(e_pcc.df.loc[e_pcc.df.orf_name.isin(e_pcc.rwr_mat_df.loc[:,i].sort_values(ascending=False)[:50].index)],
                goea,geneid2name)
        rwr_ranked_goa_results.append(goa_df_rwr_col)
    combined_df = combine_data(rwr_ranked_goa_results)
    combined_df.to_csv(snakemake.output.rwr_column,index=False)


# Save to pickle files