import src
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx
from src.enm import *
from src.utils import *
from Bio import SeqIO


figure_path = 'reports/figures/pcc_0525'
enm = Enm('enm')
enm.read_network('data/interim/costanzo_pcc_ALL',sep=',')
enm.gnm_analysis()

enm.figure_path=figure_path
enm.output_path = 'data/interim/pcc_0525/'
enm.plot_collectivity()

enm.plot_network_spring()
enm.plot_vector(sorted=True)
enm.plot_scatter(x='deg',y='eff',figure_name='deg_eff',figure_extension='pdf')
enm.plot_scatter(x='deg',y='sens',figure_name='deg_sens',figure_extension='pdf')

enm.simulate_rewire(output_name='rewire_data',save=True)
enm.simulate_rewire(output_name='rewire_data_nodegseq', save=True, nodegseq=True)
enm.simulate_rewire(output_name='rewire_data_er', save=True, nodegseq=True,random_network_type='er')

enm.plot_correlation_density(x='eff',y='deg',figure_extension='pdf')
#enm.plot_correlation_density(x='sens',y='deg',correlation='spearman',figure_extension='pdf')
nx.set_node_attributes(enm.graph_gc,dict(zip(list(enm.graph_gc.nodes),list(enm.graph_gc.nodes))),'orf_name')

go_df = pd.read_csv('data/interim/go_results/4.tsv','\t')

enm.plot_network_spring(plot_go=True,go_df_list=[go_df],level_list=[0.2])
