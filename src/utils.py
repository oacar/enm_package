import os 
import networkx as nx
import pandas as pd 
import prody
import copy
import numpy as np
import random

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

def create_goea_human(gene2go = '../data/raw/ontology/gene2go', obo_fname = '../data/raw/ontology/go-basic.obo', background='../data/interim/huri_gc_bg.tsv',**kwargs):
    
    #background='../data/raw/ontology/sgd_costanzogenes'
    obodag = GODag(obo_fname)
    from goatools.anno.genetogo_reader import Gene2GoReader

# Read NCBI's gene2go. Store annotations in a list of namedtuples
    objanno = Gene2GoReader(gene2go, taxids=[9606])
#    ns2assoc = objanno.get_ns2assc()

    ns2assoc_excl = objanno.get_ns2assc( ev_exclude = {'HGI' , 'IGI'})

    bg = pd.read_csv(background, header=None)
    bg_list = list(bg.iloc[:, 0])  # .to_list()
    
    import mygene
    mg = mygene.MyGeneInfo()
    out = mg.querymany(bg_list,scopes='ensemblgene', fields = 'unigene,entrezgene,uniprot', species=9606, as_dataframe=True)

    geneids = [int(i) for i in out.loc[bg_list,'entrezgene'].dropna().values.tolist()]
    geneid2name = pd.Series(out['entrezgene'].dropna().index.tolist(),index=[int(i) for i in out.loc[:,'entrezgene'].dropna().values]).to_dict()
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

def network_chance_deletion(list_of_nodes,Gc,out_dict):
    random_hinges = list_of_nodes#df_.sort_values('deg',ascending=False).iloc[0:i,0].tolist()
#    gc_copy = igraph_network(Gc)
    gc_copy=copy.deepcopy(Gc)
    gc_copy.delete_vertices(random_hinges)
    ig = gc_copy
    #Gc_random_hinge_deleted = nx.induced_subgraph(Gc,[n for n in Gc.nodes if n not in random_hinges])
    random_hinge_connected_components = len(ig.components())#nx.number_connected_components(Gc_random_hinge_deleted)
    Gc_hinges_removed = len(ig.components().giant().vs)#max(nx.connected_component_subgraphs(Gc_random_hinge_deleted), key=len)
    out_dict['num_of_comp'].append(random_hinge_connected_components)
    out_dict['gc_size'].append(Gc_hinges_removed)

def sequential_deletion(Gc,df_,step=10):
    deg_sorted_deletion = {'num_of_comp':[], 'gc_size':[]}
    deg_sorted_deletion_rev= {'num_of_comp':[], 'gc_size':[]}
    hinge_sorted_deletion= {'num_of_comp':[], 'gc_size':[]}
    rand_deletion =  {'num_of_comp':[], 'gc_size':[]}
    btw_sorted_deletion =  {'num_of_comp':[], 'gc_size':[]}
    eff_sorted_deletion= {'num_of_comp':[], 'gc_size':[]}
    hinge_sorted_deletion_rev= {'num_of_comp':[], 'gc_size':[]}
    sens_sorted_deletion =  {'num_of_comp':[], 'gc_size':[]}
    eigenvec_centr_sorted_deletion =  {'num_of_comp':[], 'gc_size':[]}
    closeness_centr_sorted_deletion =  {'num_of_comp':[], 'gc_size':[]}
    sens_sorted_deletion_rev =  {'num_of_comp':[], 'gc_size':[]}

    rand_nodes = random.sample([n for n in Gc.nodes],len(Gc.nodes)-1)
    range_ = range(1,len(Gc.nodes)-1,step)
    Gc = igraph_network(Gc, orf_names=df_.orf_name.values)
    #m#ax_node_num_rand_hinge_comp_list = [ ]
    for i in tqdm(range_):
        #degree
#        print(i)
        network_chance_deletion(df_.sort_values('deg',ascending=False).iloc[0:i,0].tolist(), Gc, deg_sorted_deletion)
#        print('Degree deletion successfull')

        #degree rev
        network_chance_deletion(df_.sort_values('deg',ascending=True).iloc[0:i,0].tolist(), Gc, deg_sorted_deletion_rev)

        #eigenvec_hinge
#        network_chance_deletion(df_.reindex(df_[0].abs().sort_values().index).iloc[0:i,0].tolist(), Gc, hinge_sorted_deletion)

        #random
        network_chance_deletion(rand_nodes[0:i], Gc, rand_deletion)

        #effectiveness
        network_chance_deletion(df_.reindex(df_['eff'].abs().sort_values(ascending=False).index).iloc[0:i,0].tolist(), Gc, eff_sorted_deletion)

        #eigenvec_hinge descending
#        network_chance_deletion(df_.reindex(df_[0].abs().sort_values(ascending=False).index).iloc[0:i,0].tolist(), Gc, hinge_sorted_deletion_rev)

        #sensitivity
        network_chance_deletion(df_.reindex(df_['sens'].abs().sort_values(ascending=False).index).iloc[0:i,0].tolist(), Gc, sens_sorted_deletion)

        #sensitivity rev
#        network_chance_deletion(df_.reindex(df_['sens'].abs().sort_values(ascending=True).index).iloc[0:i,0].tolist(), Gc, sens_sorted_deletion_rev)

        #betweenness
        network_chance_deletion(df_.reindex(df_['btw'].abs().sort_values(ascending=False).index).iloc[0:i,0].tolist(), Gc,btw_sorted_deletion)

        #eigenvector centrality
        network_chance_deletion(df_.reindex(df_['eigenvec_centr'].abs().sort_values(ascending=False).index).iloc[0:i,0].tolist(), Gc,eigenvec_centr_sorted_deletion)

        #closeness
        network_chance_deletion(df_.reindex(df_['closeness_centr'].abs().sort_values(ascending=False).index).iloc[0:i,0].tolist(), Gc,closeness_centr_sorted_deletion)

    dd = {'deg_sorted_deletion':deg_sorted_deletion,
#     'hinge_sorted_deletion':hinge_sorted_deletion,
     'rand_deletion':rand_deletion,
     'eff_sorted_deletion':eff_sorted_deletion,
     'sens_sorted_deletion':sens_sorted_deletion,
     'btw_sorted_deletion':btw_sorted_deletion,
#     'hinge_sorted_deletion_rev':hinge_sorted_deletion_rev,
          'closeness_centr_sorted_deletion':closeness_centr_sorted_deletion,
#          'sens_sorted_deletion_rev':sens_sorted_deletion_rev,
          'eigenvec_centr_sorted_deletion':eigenvec_centr_sorted_deletion,
        'deg_sorted_deletion_rec' : deg_sorted_deletion_rev,
        'range':range_}
    dd_df  = pd.DataFrame.from_dict({(outerKey, innerKey): values for outerKey, innerDict in dd.items() if outerKey!='range' for innerKey, values in innerDict.items()})
    dd_df['range'] = range_
    dd_df = dd_df.melt('range')
    return dd_df

def sample_nodes(sample_space, size=1):
    return   np.random.choice(sample_space, size=size, replace=False)

def get_degree_distribution(df):
    return df.groupby('deg').count().to_dict()['orf_name']   

def sample_nodes_with_degree(gnm_df , nodes_df ):
    deg_dist = get_degree_distribution(nodes_df)
    sampled_nodes = []
    for deg, count in deg_dist.items():
        sample_space = gnm_df.loc[(gnm_df.deg == deg)&(gnm_df.orf_name.isin(nodes_df.orf_name.values)==False),'orf_name'].values
        nds = sample_nodes(sample_space,size=count)
        sampled_nodes.extend(nds)
#        print(nds, deg, count, len(nds))
    return gnm_df.loc[gnm_df.orf_name.isin(sampled_nodes)]

def get_in_to_out_edge_ratio(G, nodes_df):
    flat_list_ego = np.unique([item for sublist in [list(nx.ego_graph(G,i,radius=1).nodes) for i in nodes_df.orf_name.values] for item in sublist])
    #rat = len([i for i in flat_list_ego if i in sensors['orf_name'].values ])/len([i for i in flat_list_ego if i not in sensors['orf_name'].values ])
    rat = len([i for i in flat_list_ego if i in nodes_df['orf_name'].values ])/len(flat_list_ego)
    return rat

def get_random_in_to_out_ratio(gnm_df, df , G):
    random_ratio_list = []
    for i in range(100):
        rand_nodes = sample_nodes_with_degree(gnm_df, df)
    #    print(len(rand_nodes))
        random_ratio_list.append(get_in_to_out_edge_ratio(G, rand_nodes))
    return random_ratio_list

def get_subnetwork(gc, subset):
    neighbors =[[n for n in nx.neighbors(gc, i)] for i in subset]
    flat_list = [item for sublist in neighbors for item in sublist]
    flat_list.extend(subset)
    sub_gc = nx.induced_subgraph(gc, flat_list).copy()
    for (n,d) in sub_gc.nodes(data=True):
        del d["pos"]
    return sub_gc


def get_maximal_subsets(sets):
    sets = sorted(map(set, sets), key=len,reverse=True)
    maximal_subsets = []
    for s in sets:
        if not any(maximal_subset.issuperset(s) for maximal_subset in maximal_subsets):
            maximal_subsets.append(s)

    return maximal_subsets
 
