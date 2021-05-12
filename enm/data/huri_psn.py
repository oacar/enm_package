import enm
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx
from enm.enm import *
from enm.utils import *


#huri = Enm('huri')
huri_pickle = 'data/interim/huri_1205/huri'
with open(f"{huri_pickle}.pickle",'rb') as f:
    e_huri = pickle.load(f)
#e_huri.get_jaccard_mat()

#with open(f"{huri_pickle}_jaccardmat.pickle",'wb') as f:
#    pickle.dump(e_huri.jaccard_mat, f)

#with open(f"{huri_pickle}_jaccardmat.pickle",'rb') as f:
#    jaccard_mat = pickle.load(f)

psn = jaccard_mat
sigma = 0.2
psn_thr = psn>sigma

psn_g = nx.from_numpy_array(psn_thr)
node_name_map = dict(zip(range(len(e_huri.nodes)),e_huri.nodes))
psn_g = nx.relabel_nodes(psn_g,node_name_map)

e_huri.psn = psn_g
#
with open(f"{huri_pickle}_withpsn.pickle",'wb') as f:
    pickle.dump(e_huri,f, protocol=4)
#huri.read_network('data/interim/huri/huri_cleared.csv')
#huri.gnm_analysis(normalized=False)
