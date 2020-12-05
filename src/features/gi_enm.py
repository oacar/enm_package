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
enm = Enm('enm')
enm.read_network('data/interim/gi/gi_.08.csv')
enm.gnm_analysis(normalized=False)

#enm.figure_path=figure_path
enm.output_path = 'data/interim/gi_1201/'
os.mkdir(enm.output_path)
enm.get_category('data/interim/strain_ids_with_experiment_count_all.csv')
enm.spring_pos()

os.mkdir(f"{enm.output_path}/random_dfs")
enm.simulate_rewire(output_name='rewire_data',save=True, normalized=False,sim_num=10)
random_dfs = [enm.e_list[i].df.to_csv(f'data/interim/gi_1201/random_dfs/rand_{i}.csv') for i in range(10)]
#enm.simulate_rewire(output_name='rewire_data_nodegseq', save=True, nodegseq=True, normalized=False, sim_num=10)
#enm.simulate_rewire(output_name='rewire_data_er', save=True, nodegseq=True,random_network_type='er', normalized=False)

neigbor_btw = []
neighbor_degree = []
for i in enm.graph_gc.nodes:
    neigbor_btw.append(np.average([enm.df.loc[enm.df.orf_name==a,'btw'].values for a in nx.neighbors(enm.graph_gc,i)]))
    neighbor_degree.append(np.average([enm.df.loc[enm.df.orf_name==a,'deg'].values for a in nx.neighbors(enm.graph_gc,i)]))

enm.df['neighbor_btw'] = neigbor_btw
enm.df['neighbor_degree'] = neighbor_degree

with open(f"{enm.output_path}/gi_100.pickle",'wb') as f:
    pickle.dump(enm,f, protocol=4)


##### very strong#####

enm_02 = Enm('enm_02')
enm_02.read_network('data/interim/gi/gi_.2.csv')
enm_02.gnm_analysis(normalized=False)

#enm.figure_path=figure_path
enm_02.output_path = 'data/interim/gi_1205_0.2/'
os.mkdir(enm_02.output_path)
enm_02.get_category('data/interim/strain_ids_with_experiment_count_all.csv')
enm_02.spring_pos()

os.mkdir(f"{enm_02.output_path}/random_dfs")
enm_02.simulate_rewire(output_name='rewire_data',save=True, normalized=False,sim_num=10)
random_dfs = [enm_02.e_list[i].df.to_csv(f'{enm_02.output_path}/random_dfs/rand_{i}.csv') for i in range(10)]
#enm.simulate_rewire(output_name='rewire_data_nodegseq', save=True, nodegseq=True, normalized=False, sim_num=10)
#enm.simulate_rewire(output_name='rewire_data_er', save=True, nodegseq=True,random_network_type='er', normalized=False)

neigbor_btw = []
neighbor_degree = []
for i in enm_02.graph_gc.nodes:
    neigbor_btw.append(np.average([enm_02.df.loc[enm_02.df.orf_name==a,'btw'].values for a in nx.neighbors(enm_02.graph_gc,i)]))
    neighbor_degree.append(np.average([enm_02.df.loc[enm_02.df.orf_name==a,'deg'].values for a in nx.neighbors(enm_02.graph_gc,i)]))

enm_02.df['neighbor_btw'] = neigbor_btw
enm_02.df['neighbor_degree'] = neighbor_degree

with open(f"{enm_02.output_path}/gi_100.pickle",'wb') as f:
    pickle.dump(enm_02,f, protocol=4)

