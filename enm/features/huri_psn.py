import src
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx
import pickle
from tqdm import tqdm
from src.enm import *
from src.utils import *



#figure_path = 'reports/figures/pcc_0603'
huri_pickle = 'data/interim/huri_1205/huri'
with open(f"{huri_pickle}_withpsn.pickle",'rb') as f:
    e_huri = pickle.load(f)
enm = Enm('enm')
enm.G = e_huri.psn
enm.giant_component()
enm.gnm_analysis(normalized=False)

#enm.figure_path=figure_path
enm.output_path = 'data/interim/psn_1201/'
if os.path.exists(enm.output_path)==False:
    os.mkdir(enm.output_path)
#enm.get_category('data/interim/strain_ids_with_experiment_count_all.csv')
enm.spring_pos()
if os.path.exists(f"{enm.output_path}/random_dfs")==False:
    os.mkdir(f"{enm.output_path}/random_dfs")
enm.simulate_rewire(output_name='rewire_data',save=True, normalized=False,sim_num=10)
random_dfs = [enm.e_list[i].df.to_csv(f'{enm.output_path}/random_dfs/rand_{i}.csv') for i in range(10)]
#enm.simulate_rewire(output_name='rewire_data_nodegseq', save=True, nodegseq=True, normalized=False, sim_num=10)
#enm.simulate_rewire(output_name='rewire_data_er', save=True, nodegseq=True,random_network_type='er', normalized=False)

neigbor_btw = []
neighbor_degree = []
for i in enm.graph_gc.nodes:
    neigbor_btw.append(np.average([enm.df.loc[enm.df.orf_name==a,'btw'].values for a in nx.neighbors(enm.graph_gc,i)]))
    neighbor_degree.append(np.average([enm.df.loc[enm.df.orf_name==a,'deg'].values for a in nx.neighbors(enm.graph_gc,i)]))

enm.df['neighbor_btw'] = neigbor_btw
enm.df['neighbor_degree'] = neighbor_degree

with open(f"{enm.output_path}/psn.pickle",'wb') as f:
    pickle.dump(enm,f, protocol=4)


