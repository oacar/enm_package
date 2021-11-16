import numpy as np
import networkx as nx
import pickle
from enm.Enm import *
from enm.utils import *

#snakemake = {'input': {'gaf': '../data/raw/ontology/sgd.gaf' , 'obo':'../data/raw/ontology/go-basic.obo', "background_file" :'../data/raw/ontology/sgd_costanzogenes', "sgd_info": '../data/raw/ontology/SGD_features.tab','pickle_file_name':'../data/interim_0.2/pcc.pickle'}}
gaf = snakemake.input['gaf']
obo = snakemake.input['obo']
background_file = snakemake.input['background_file']
sgd_info = snakemake.input['sgd_info']
pickle_file = snakemake.input['pickle_file_name']
with open(pickle_file,'rb') as f:
    enm = pickle.load(f)

enm.get_sensor_effector(use_threshold=True)

print('go analysis')
goea, geneid2name,obodag = create_goea(gaf = gaf, obo_fname=obo, 
                                background=background_file, sgd_info_tab = sgd_info)
enm.analyze_components_biology(goea, geneid2name, True)
enm.analyze_components_biology(goea, geneid2name, False)

sensors_df = enm.sensors_df.reset_index(drop=True)
sensors_df['sensor_cluster'] = sensors_df.fillna(value={'sensor_cluster':'Unclustered'}).sensor_cluster.astype("category")

order = sensors_df.groupby('sensor_cluster').median()['sens'].sort_values().index.values
sensors_df['sensor_cluster'].cat.reorder_categories(list(order),ordered=True, inplace=True)
sensors_df['gid'] = sensors_df.groupby('sensor_cluster').ngroup()
sensors_df['cluster_or_go'] = ['Unclustered' if sensors_df.at[i,'sensor_cluster'] == 'Unclustered'
                            else sensors_df.at[i,'gid'] if pd.isna(sensors_df.at[i,'go_group']) 
                            else sensors_df.at[i,'go_group'] for i in range(sensors_df.shape[0])]
di = {"cellular response to iron ion starvation":"Iron ion\ntransport" ,
"mitochondria-nucleus signaling pathway":"Mitochondria nucleus\nsignaling pathway",
"phenylalanine transport":"Phenylalanine\ntransport",
"hexose metabolic process":"Hexose metabolic\nprocess",
"tricarboxylic acid cycle":"Tricarboxylic acid\ncycle"}
sensors_df.replace({'cluster_or_go': di},inplace=True)

sensors_df['label'] = ["Unclustered" if sensors_df.at[i,'cluster_or_go']=='Unclustered' 
                                    else f"SC{sensors_df.at[i,'gid']}" if type(sensors_df.at[i,'cluster_or_go'])==np.int64
                                    else f"SC{sensors_df.at[i,'gid']}\n{sensors_df.at[i,'cluster_or_go']}" for i in range(sensors_df.shape[0])]

sensors_df.to_csv(snakemake.output.sensors_df_fname,index=False)
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
