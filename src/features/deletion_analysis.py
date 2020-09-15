import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from prody import *
import logging
import seaborn as sns
from scipy.stats import pearsonr
from itertools import chain
import os
from scipy.stats import spearmanr
import pymysql
import itertools as itr
from datetime import datetime
from mlxtend.evaluate import permutation_test
import pickle
import random
from tqdm import tqdm
from src.utils import sequential_deletion
from src.enm import Enm

with open(f'data/interim/pcc_0909/pcc.pickle','rb') as f:
    e_pcc = pickle.load(f)
with open(f'data/interim/yuri/yuri.pickle','rb') as f:
    e_yuri = pickle.load(f)
with open(f'data/interim/huri/huri.pickle','rb') as f:
    e_huri = pickle.load(f)

dd_pcc = sequential_deletion(e_pcc.graph_gc, e_pcc.df,step=1)
dd_huri = sequential_deletion(e_huri.graph_gc, e_huri.df,step=1)
dd_yuri = sequential_deletion(e_yuri.graph_gc, e_yuri.df,step=1)

result_file = 'data/interim/deletion_analysis_0911'
dd_pcc.to_csv(f'{result_file}/pcc_deletion_100.csv')
dd_huri.to_csv(f'{result_file}/huri_deletion_100.csv')
dd_yuri.to_csv(f'{result_file}/yuri_deletion_100.csv')

