import pandas as pd
import numpy as np
import os
#snakemake.input[0]

network_edgelist = pd.read_csv("../data/raw/yuri/Y2H_union.txt",sep='\t',header=None)
network_edgelist['weight'] = 1
network_edgelist.rename(columns={0:"gene1",1:"gene2"}).to_csv('../data/yeast_y2h/y2h_edgelist.csv',index=False)

genelist = [*network_edgelist.iloc[:,0].unique(),*network_edgelist.iloc[:,1].unique()]
sgd = pd.read_csv('../data/raw/ontology/SGD_features.tab',sep='\t',header=None)

sgd.loc[sgd.iloc[:,3].isin(genelist),0].to_csv('../data/yeast_y2h/y2h_background.csv',index=False, header=None)
