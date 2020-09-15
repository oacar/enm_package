import src
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx
import pickle
from tqdm import tqdm
from src.enm import *
from src.utils import *
from Bio import SeqIO


figure_path = 'reports/figures/gi_0827'
enm = Enm('enm')
enm.read_network('data/interim/SGA_02_withoverlapping.csv',sep=',')
enm.gnm_analysis(normalized=False)

enm.figure_path=figure_path
enm.output_path = 'data/interim/gi_0827/'
enm.spring_pos()

enm.simulate_rewire(output_name='rewire_data',save=True, normalized=False)

neigbor_btw = []
neighbor_degree = []
for i in enm.graph_gc.nodes:
    neigbor_btw.append(np.average([enm.df.loc[enm.df.orf_name==a,'btw'].values for a in nx.neighbors(enm.graph_gc,i)]))
    neighbor_degree.append(np.average([enm.df.loc[enm.df.orf_name==a,'deg'].values for a in nx.neighbors(enm.graph_gc,i)]))

enm.df['neighbor_btw'] = neigbor_btw
enm.df['neighbor_degree'] = neighbor_degree

with open(f"{enm.output_path}/gi.pickle",'wb') as f:
    pickle.dump(enm,f)


