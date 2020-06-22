import os 
import networkx as nx
import pandas as pd 
import prody
import numpy as np

from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.anno.gaf_reader import GafReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS


from .enm import *

def create_goea(gaf = '../data/raw/ontology/sgd.gaf', obo_fname = '../data/raw/ontology/go-basic.obo', background='../data/raw/ontology/sgd_costanzogenes', sgd_info_tab='../data/raw/ontology/SGD_features.tab',species='yeast',**kwargs):
    
    #background='../data/raw/ontology/sgd_costanzogenes'
    obodag = GODag(obo_fname)
    objanno = GafReader(gaf)
#    ns2assoc = objanno.get_ns2assc()

    ns2assoc_excl = objanno.get_ns2assc( ev_exclude = {'HGI' , 'IGI'})

    bg = pd.read_csv(background, header=None)
    bg_list = list(bg.iloc[:, 0])  # .to_list()
    geneids = bg_list
    sgd_info = pd.read_csv(sgd_info_tab, sep='\t',header=None)

    geneid2name = pd.Series(sgd_info.iloc[:,3].values,index=sgd_info.iloc[:,0]).to_dict()
    goeaobj = GOEnrichmentStudyNS(
        geneids,  # List of mouse protein-coding genes
        ns2assoc_excl,  # geneid/GO associations
        obodag,  # Ontologies
        propagate_counts=False,
        alpha=0.05,  # default significance cut-off
        methods=['fdr_bh'], prt=None)

    return goeaobj, geneid2name#, objanno, ns2assoc, ns2assoc_excl

def go_findenrichment( query, outname=None, species='yeast',**kwargs):
    goeaobj, geneid2name = create_goea(**kwargs)
    query_gene_ids = [key for key,value in geneid2name.items() if value in query]#sgd_info[sgd_info.iloc[:,3].isin(query)].iloc[:,0].values.tolist()

    goea_results_all = goeaobj.run_study(query_gene_ids, prt=None)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]

#    if outname is not None:
    goeaobj.wr_tsv(outname, goea_results_all, itemid2name=geneid2name)
    go_df = goea_to_pandas(goea_results_sig,geneid2name)
    return goeaobj, goea_results_sig,geneid2name, go_df

def goea_to_pandas(goea_results_sig, geneid2name):
    """ Converts goea object from goatools GO enrichment test to a Pandas dataframe
    :param goea_results_sig: Significant GO term objects
    :type goea_results_sig: list of GOEnrichmentStudy
    :return: Dataframe
    :rtype: Pandas DataFrame
    """
    if len(goea_results_sig) == 0 :
        return None
    go_df_n = pd.DataFrame([[getattr(i,x) for x in i.get_prtflds_default()] for i  in goea_results_sig], columns=goea_results_sig[0].get_prtflds_default())
    orf_names = []
    for i in go_df_n.study_items:
        orf_names.append([geneid2name[_id] for _id in i])
    go_df_n.study_items = orf_names
    return go_df_n 
