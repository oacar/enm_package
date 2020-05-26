import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from src.visualize.visualize import plot_correlation_density
import pickle

#figure_path = ''
# Random collectivity
with open('data/interim/pcc_0525/rewire_data_er.pickle','rb') as f:
    e_er = pickle.load(f)
with open('data/interim/pcc_0525/pcc.pickle','rb') as f:
    e_pcc = pickle.load(f)
with open('data/interim/pcc_0525/rewire_data_nodegseq.pickle','rb') as f:
    e_nodegseq = pickle.load(f)
figure_path = e_er[0].figure_path = e_pcc.figure_path = e_nodegseq[0].figure_path = 'reports/figures/paper_figures_try'
e_er[0].plot_collectivity(figure_name = 'er_coll', figure_extension='pdf')
e_pcc.plot_collectivity(figure_name = 'pcc_coll', figure_extension='pdf')
e_nodegseq[0].plot_collectivity(figure_name = 'nodegseq_coll', figure_extension='pdf')

# correlation_density
er_df = pd.read_csv('data/interim/pcc_0525/rewire_data_er.csv')
nodegseq_df = pd.read_csv('data/interim/pcc_0525/rewire_data_nodegseq.csv')
#e_pcc.plot_correlation_density(x='eff',y='deg',figure_extension='pdf')
plot_correlation_density(e_pcc.df, [e_pcc.rewire_df,nodegseq_df,er_df], x='eff', y='deg', figure_path=figure_path, figure_extension='pdf')
plot_correlation_density(e_pcc.df, [e_pcc.rewire_df,nodegseq_df,er_df], x='sens', y='deg', figure_path=figure_path, figure_extension='pdf',correlation='spearman')
