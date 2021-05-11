#import src
import pandas as pd
import numpy as np
import networkx as nx
import pickle
from tqdm import tqdm
from src.enm import *
from src.utils import *
import matplotlib.pyplot as plt


#figure_path = 'reports/figures/pcc_0603'
enm = Enm('enm')
enm.read_network('data/interim/costanzo_pcc_ALL',sep=',')
enm.gnm_analysis(normalized=False)

#enm.figure_path=figure_path
enm.output_path = 'data/interim/pcc_021521_withclustering/'
enm.get_category('data/interim/strain_ids_with_experiment_count_all.csv')
enm.spring_pos()
print('clustering')
enm.get_sensor_effector(use_threshold=False)

print('go analysis')
goea, geneid2name = create_goea(gaf = 'data/raw/ontology/sgd.gaf', obo_fname='data/raw/ontology/go-basic.obo', 
                                background='data/interim/costanzo_gc_bg.tsv', sgd_info_tab = 'data/raw/ontology/SGD_features.tab')
enm.analyze_components_biology(goea, geneid2name, True)
enm.analyze_components_biology(goea, geneid2name, False)
print('rewiring')
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

with open(f"{enm.output_path}/pcc_10.pickle",'wb') as f:
    pickle.dump(enm,f, protocol=4)


