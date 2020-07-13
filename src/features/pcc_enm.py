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


figure_path = 'reports/figures/pcc_0603'
enm = Enm('enm')
enm.read_network('data/interim/costanzo_pcc_ALL',sep=',')
enm.gnm_analysis(normalized=True)

enm.figure_path=figure_path
enm.output_path = 'data/interim/pcc_0603/'

enm.simulate_rewire(output_name='rewire_data',save=True, normalized=True)
enm.simulate_rewire(output_name='rewire_data_nodegseq', save=True, nodegseq=True, normalized=True)
enm.simulate_rewire(output_name='rewire_data_er', save=True, nodegseq=True,random_network_type='er', normalized=True)

with open(f"{enm.output_path}/pcc.pickle",'wb') as f:
    pickle.dump(enm,f)


## ---- 
from src.visualize.visualize import plot_lambda_collectivity

plot_lambda_collectivity(enm.gnm.getEigvals(),enm.coll,figure_path)


#enm
enm.fiedler_evaluation()

adj = nx.adjacency_matrix(gc).todense()
orf_names = self.df.orf_name
g_ig = igraph.Graph.Adjacency((adj > 0).tolist())
g_ig.vs['name'] = orf_names
#g_ig.set_vertec
gc = self.graph_gc
eigvecs = self.gnm.getEigvecs()
figsize = kwargs.pop('figsize',(12,6))

lost_edges = [ ]

for i in tqdm(range(5182)):
    cl1 = [orf_names[i] for i, j in enumerate(eigvecs[:,i]) if j>0]
    cl2 = [orf_names[i] for i, j in enumerate(eigvecs[:,i]) if j<0]
    cg1 = g_ig.induced_subgraph(cl1)
    cg2 = g_ig.induced_subgraph(cl2)
#    cg1 = nx.induced_subgraph(gc, cl1)
#    cg2 = nx.induced_subgraph(gc,cl2)
    lost_edges.append( len(gc.edges)-(len(cg2.es)+len(cg1.es)))

lost_edges[:10]
#enm.lost_edges[:10]
enm.fiedler_evaluation(plot=True)
src.visualize.visualize.plot_fiedler_data(enm,figsize=(8,8))
