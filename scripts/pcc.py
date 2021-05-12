import numpy as np
import networkx as nx
import pickle
from enm.enm import *
from enm.utils import *

network_file = snakemake.input['network_file']
output_path = snakemake.params.output_path
#figure_path = 'reports/figures/pcc_0603'
enm = Enm('enm')
enm.read_network(network_file,sep=',')
enm.gnm_analysis(normalized=False)

neigbor_btw = []
neighbor_degree = []
for i in enm.graph_gc.nodes:
    neigbor_btw.append(np.average([enm.df.loc[enm.df.orf_name==a,'btw'].values for a in nx.neighbors(enm.graph_gc,i)]))
    neighbor_degree.append(np.average([enm.df.loc[enm.df.orf_name==a,'deg'].values for a in nx.neighbors(enm.graph_gc,i)]))

enm.df['neighbor_btw'] = neigbor_btw
enm.df['neighbor_degree'] = neighbor_degree

#enm.figure_path=figure_path
enm.output_path = output_path
#enm.get_category('data/interim/strain_ids_with_experiment_count_all.csv')
enm.spring_pos()

with open(f"{enm.output_path}/pcc_10.pickle",'wb') as f:
    pickle.dump(enm,f, protocol=4)


