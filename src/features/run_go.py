import pickle
import numpy as np
import pandas as pd
from src.enm import Enm
from src.utils import create_goea
from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.anno.gaf_reader import GafReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from tqdm import tqdm




with open('../data/interim/pcc.pickle','rb') as handle:
    e = pickle.load(handle)

strain_ids = pd.read_csv('../data/interim/strain_ids_with_experiment_count_all.csv')
combined_df = pd.merge(e.df, strain_ids, left_on='orf_name',right_on='Allele Gene name')
combined_df['group']=np.where(combined_df.cat.isna(),'essential','nonessential')

goeaobj, geneid2name = create_goea()
sim_num = len(e.coll)
#go_results = []
for i in tqdm(range(sim_num)):
    eig_id = e.coll_index_sorted[i]
    eig_vec = e.gnm.getEigvecs()[:,eig_id]
    query = combined_df.iloc[np.argsort(np.abs(eig_vec)),:].loc[:,'Systematic gene name'].values[:100]#combined_df.sample(100).loc[:,'Systematic gene name'].values#

    #print(geneid2name)
    query_gene_ids = [key for key,value in geneid2name.items() if value in query]#sgd_info[sgd_info.iloc[:,3].isin(query)].iloc[:,0].values.tolist()

    goea_results_all = goeaobj.run_study(query_gene_ids, prt=None)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]
    goeaobj.wr_tsv(f'../data/interim/go_results_100hinges/{i}.tsv', goea_results_sig , itemid2name = geneid2name)
    #go_results.append(goea_results_sig)
    #s, q, geneid2name = go_findenrichment(None,)
#go = pd.read_csv('tmp.tsv','\t')
#go