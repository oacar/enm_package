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

def create_goea(obo_fname = '../data/raw/ontology/go-basic.obo', background='../data/raw/ontology/sgd_costanzogenes', species='yeast',**kwargs):
    
    #background='../data/raw/ontology/sgd_costanzogenes'
    gaf = '../data/raw/ontology/sgd.gaf'
    obodag = GODag(obo_fname)
    objanno = GafReader(gaf)
    ns2assoc = objanno.get_ns2assc()

    bg = pd.read_csv(background, header=None)
    bg_list = list(bg.iloc[:, 0])  # .to_list()
    geneids = bg_list
    sgd_info = pd.read_csv('../data/raw/ontology/SGD_features.tab', sep='\t',header=None)

    geneid2name = pd.Series(sgd_info.iloc[:,3].values,index=sgd_info.iloc[:,0]).to_dict()
    goeaobj = GOEnrichmentStudyNS(
        geneids,  # List of mouse protein-coding genes
        ns2assoc,  # geneid/GO associations
        obodag,  # Ontologies
        propagate_counts=False,
        alpha=0.05,  # default significance cut-off
        methods=['fdr_bh'], prt=None)
    return goeaobj, geneid2name

def go_findenrichment( query, outname=None,background='../data/raw/ontology/sgd_costanzogenes', species='yeast',**kwargs):
    goeaobj, geneid2name = create_goea()
    query_gene_ids = [key for key,value in geneid2name.items() if value in query]#sgd_info[sgd_info.iloc[:,3].isin(query)].iloc[:,0].values.tolist()

    goea_results_all = goeaobj.run_study(query_gene_ids, prt=None)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]

#    if outname is not None:
    goeaobj.wr_tsv(outname, goea_results_sig, itemid2name=geneid2name)

    return goeaobj, goea_results_sig,geneid2name