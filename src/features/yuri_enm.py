import src
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx
from src.enm import *
from src.utils import *
from Bio import SeqIO


figure_path = 'reports/figures/yuri_0519'
yuri = Enm('yuri')
yuri.read_network('data/interim/yuri/yuri_combined.csv')
yuri.gnm_analysis()

yuri.figure_path=figure_path
yuri.output_path = 'data/interim/yuri/'
yuri.plot_collectivity()

yuri.plot_network_spring()
yuri.plot_vector(sorted=True)
yuri.plot_scatter(x='deg',y='eff',figure_name='deg_eff')
yuri.plot_scatter(x='deg',y='sens',figure_name='deg_sens')

yuri.simulate_rewire(output_name='rewire_data',save=True)
yuri.rewire_df
yuri.plot_correlation_density(x='eff',y='deg')
yuri.plot_correlation_density(x='sens',y='deg',correlation='spearman')
