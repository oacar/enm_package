import src
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import networkx as nx
import pickle
from src.enm import *
from src.utils import *

from Bio import SeqIO


figure_path = 'reports/figures/pcc_0525'
enm = Enm('enm')
enm.read_network('data/interim/costanzo_pcc_ALL',sep=',')
enm.gnm_analysis()

enm.figure_path=figure_path
enm.output_path = 'data/interim/pcc_0525/'

enm.simulate_rewire(output_name='rewire_data',save=True)
enm.simulate_rewire(output_name='rewire_data_nodegseq', save=True, nodegseq=True)
enm.simulate_rewire(output_name='rewire_data_er', save=True, nodegseq=True,random_network_type='er')

with open(f"{enm.output_path}/pcc.pickle",'wb') as f:
    pickle.dump(enm,f)


