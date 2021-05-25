import numpy as np
import networkx as nx
import pickle
from enm.Enm import *
from enm.utils import *

gaf = snakemake.input['gaf']
obo = snakemake.input['obo']
background_file = snakemake.input['background_file']
sgd_info = snakemake.input['sgd_info']
pickle_file = snakemake.input['pickle_file_name']
with open(pickle_file,'rb') as f:
    enm = pickle.load(f)

enm.get_sensor_effector(use_threshold=True)

print('go analysis')
goea, geneid2name = create_goea(gaf = gaf, obo_fname=obo, 
                                background=background_file, sgd_info_tab = sgd_info)
enm.analyze_components_biology(goea, geneid2name, True)
enm.analyze_components_biology(goea, geneid2name, False)

enm.sensors_df.to_csv(snakemake.output.sensors_df_fname,index=False)
enm.effectors_df.to_csv(snakemake.output.effectors_df_fname,index=False)

effector_sensor_combined_go_df = pd.DataFrame()
for i,j in enm.go_groups['sensors_go_groups'].items():
    if j is not None:
        j['cluster_id']=i
        j['cluster_type']='sensor'
        effector_sensor_combined_go_df = effector_sensor_combined_go_df.append(j)
for i,j in enm.go_groups['effectors_go_groups'].items():
    if j is not None:
        j['cluster_id']=i
        j['cluster_type']='effector'
        effector_sensor_combined_go_df = effector_sensor_combined_go_df.append(j)
effector_sensor_combined_go_df.to_csv(snakemake.output.effector_sensor_combined_go_df,index=False)