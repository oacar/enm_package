import numpy as np
import networkx as nx
import pickle
from enm.enm import *
from enm.utils import *

gaf = snakemake.input['gaf']
obo = snakemake.input['obo']
background_file = snakemake.input['background_file']
sgd_info = snakemake.input['sgd_info']
pickle_file = snakemake.input[0]
with open(pickle_file,'rb') as f:
    enm = pickle.load(f)

print('clustering')
enm.get_sensor_effector(use_threshold=False)

print('go analysis')
goea, geneid2name = create_goea(gaf = gaf, obo_fname=obo, 
                                background=background_file, sgd_info_tab = sgd_info)
enm.analyze_components_biology(goea, geneid2name, True)
enm.analyze_components_biology(goea, geneid2name, False)

print('rewiring')
#Don't run
enm.simulate_rewire(output_name=snakemake.output['rewired_dataframe'],save=True, normalized=False,sim_num=snakemake.params.n_sim)

with open(snakemake.output['pickle_file'],'wb') as f:
    pickle.dump(enm,f, protocol=4)