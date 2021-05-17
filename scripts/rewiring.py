import numpy as np
import networkx as nx
import pickle
from enm.Enm import *
from enm.utils import *

pickle_file = snakemake.input[0]
with open(pickle_file,'rb') as f:
    enm = pickle.load(f)
print('rewiring')
#Don't run with n>10, takes long
enm.simulate_rewire(output_name='tmp',save=False, normalized=False,sim_num=snakemake.params.n_sim)
df = pd.DataFrame()
for i in range(snakemake.params.n_sim):
    df_tmp = enm.e_list[i].df
    df_tmp['random_id'] = i
    df = df.append(df_tmp)

df.to_csv(snakemake.output.pcc_df_random, index=False)
#with open(snakemake.output['pickle_file'],'wb') as f:
#    pickle.dump(enm,f, protocol=4)