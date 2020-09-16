import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx
import numpy as np
import pandas as pd
import plotnine as p9
import pickle 
import os
import re
import itertools as itr
from src.enm import Enm
from src.utils import *

with open(f'data/interim/pcc_0909/pcc.pickle','rb') as f:
    e_pcc = pickle.load(f)
with open(f'data/interim/yuri/yuri.pickle','rb') as f:
    e_yuri = pickle.load(f)
with open(f'data/interim/huri/huri.pickle','rb') as f:
    e_huri = pickle.load(f)

pcc_df = e_pcc.df
yuri_df = e_yuri.df
huri_df = e_huri.df

sensors = pcc_df.loc[pcc_df.sens>np.quantile(pcc_df.sens,0.99)]

pcc_random_ratio_list = get_random_in_to_out_ratio(pcc_df, pcc_df.loc[pcc_df.sens>np.quantile(pcc_df.sens,0.99)], e_pcc.graph_gc)
pcc_ratio = get_in_to_out_edge_ratio(e_pcc.graph_gc, pcc_df.loc[pcc_df.sens>np.quantile(pcc_df.sens,0.99)])

huri_random_ratio_list = get_random_in_to_out_ratio(huri_df, huri_df.loc[huri_df.sens>np.quantile(huri_df.sens,0.99)], e_huri.graph_gc)
huri_ratio = get_in_to_out_edge_ratio(e_huri.graph_gc, huri_df.loc[huri_df.sens>np.quantile(huri_df.sens,0.99)])

yuri_random_ratio_list = get_random_in_to_out_ratio(yuri_df, yuri_df.loc[yuri_df.sens>np.quantile(yuri_df.sens,0.99)], e_yuri.graph_gc)
yuri_ratio = get_in_to_out_edge_ratio(e_yuri.graph_gc, yuri_df.loc[yuri_df.sens>np.quantile(yuri_df.sens,0.99)])

pd.concat(
[
    pd.DataFrame([[('pcc'),pcc_ratio, ('real')]], columns=['data','ratio','type']),
    pd.DataFrame(list(zip(itr.repeat('pcc'),pcc_random_ratio_list, itr.repeat('random'))), columns=['data','ratio','type']),
    pd.DataFrame([[('huri'),huri_ratio, ('real')]], columns=['data','ratio','type']),
    pd.DataFrame(list(zip(itr.repeat('huri'),huri_random_ratio_list, itr.repeat('random'))), columns=['data','ratio','type']),
    pd.DataFrame([[('yuri'),yuri_ratio, ('real')]], columns=['data','ratio','type']),
    pd.DataFrame(list(zip(itr.repeat('yuri'),yuri_random_ratio_list, itr.repeat('random'))), columns=['data','ratio','type'])
        ]
).to_csv('data/interim/sensor_connectivity_0915/sensor_connectivity_df.csv')


