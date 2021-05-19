import pandas as pd
import pickle 
import itertools as itr
from enm.Enm import Enm
from enm.utils import *

with open(snakemake.input[0],'rb') as f:
    e_pcc = pickle.load(f)

pcc_df = e_pcc.df
#sensors = pcc_df.loc[pcc_df.sens>np.quantile(pcc_df.sens,0.99)]

pcc_random_ratio_list = get_random_in_to_out_ratio(pcc_df, pcc_df.loc[pcc_df.sens>np.quantile(pcc_df.sens,0.99)], e_pcc.graph_gc)
pcc_ratio = get_in_to_out_edge_ratio(e_pcc.graph_gc, pcc_df.loc[pcc_df.sens>np.quantile(pcc_df.sens,0.99)])

pd.concat(
[
    pd.DataFrame([[('pcc'),pcc_ratio, ('real')]], columns=['data','ratio','type']),
    pd.DataFrame(list(zip(itr.repeat('pcc'),pcc_random_ratio_list, itr.repeat('random'))), columns=['data','ratio','type']),
        ]
).to_csv(snakemake.output[0])



