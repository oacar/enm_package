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
    """Get a query dataframe with set of orfs and return significant go results by first creating goa analysis object and running enrichment analysis
    use this function instead of `query_goatools` to create goatools object and write results to file directly

    :param query: dataframe for query genes
    :type query: pandas dataframe
    :param outname: tsv filename to write results, if None the results will be written to stdio, defaults to None
    :type outname: str, optional
    :param species: which species is under question. depreceated, defaults to 'yeast'
    :type species: str, optional
    """
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


def query_goatools(query, goea,geneid2name):
    """get query dataframe and goa files and return enrichments

    :param query: query gene dataframe
    :type query: pandas dataframe
    :param goea: goa object that will be used to run gene ontology analysis using GOAtools
    :type goea: goatools.goea.go_enrichment_ns.GOEnrichmentStudyNS
    :param geneid2name: dictionary to map geneids to gene names. needed to convert systematic names to intended names
    :type geneid2name: dict
    :return: enrichment dataframe 
    :rtype: pandas dataframe
    """
    query_gene_ids = [key for key,value in geneid2name.items() if value in query.loc[:,'Systematic gene name'].unique()]
    goea_res_all = goea.run_study(query_gene_ids)
    goea_res_sig = [r for r in goea_res_all if r.p_fdr_bh <0.1]
    go_df_sensor = goea_to_pandas(goea_res_sig, geneid2name)
    return go_df_sensor

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
    btw_edges = len(nx.induced_subgraph(G,nodes_df.orf_name).edges)
    flat_list_ego = np.unique([item for sublist in [list(nx.ego_graph(G,i,radius=1).nodes) for i in nodes_df.orf_name.values] for item in sublist])
    total_edges = len(nx.induced_subgraph(G,flat_list_ego).edges)- len(nx.induced_subgraph(G,np.setdiff1d(flat_list_ego,nodes_df.orf_name)).edges)
    #rat = len([i for i in flat_list_ego if i in sensors['orf_name'].values ])/len([i for i in flat_list_ego if i not in sensors['orf_name'].values ])
    rat = btw_edges/total_edges #len([i for i in flat_list_ego if i in nodes_df['orf_name'].values ])/len(flat_list_ego)
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

def jaccard_index(G, a,b):
    a_neigh = [i for i in nx.neighbors(G, a)]
    b_neigh = [i for i in nx.neighbors(G, b)]
    union = np.union1d(a_neigh,b_neigh)
    intersect = np.intersect1d(a_neigh,b_neigh)
    jaccard = len(intersect)/len(union)
    return jaccard

def get_max_path(enm, target_cluster, source_cluster):
    """Calculate max information carrying path between 2 clusters by finding all possible nodes and returns based on maximum distance

    :param enm: this will be used to calculate prs weighted paths
    :type enm: Enm object
    :param target_cluster: cluster which contains possible target nodes
    :type target_cluster: networkx graph
    :param source_cluster: cluster which contains possible source nodes
    :type source_cluster: networkx graph
    """
    wmin = np.inf
    lmin = 0
    for source in source_cluster.nodes:
        for target in target_cluster.nodes:
            w, l1 = enm.get_prs_weighted_path(source,target)
            #print(w)
            num_of_sources = len([i for i in l1 if i in source_cluster.nodes])
            if num_of_sources > 1 :
                continue
            if w < wmin :
                lmin=l1
    #if lmin == 0 :
     #   lmin = l1
    return wmin, lmin
def get_path_positions(enm, sensors_sub, effectors_sub):
    """
        calculate positions for effectors and sensors
        choose a random source and target form these, respectively
        get PRS path between source and target
        put the nodes on the path linearly in between clusters
        change the source and target positions based on the initial cluster positions

        :param enm: this will be used to calculate prs weighted paths
        :type enm: Enm object
    """
    sensor_pos = nx.spring_layout(sensors_sub,center=(5,0))
    effector_pos = nx.spring_layout(effectors_sub, center=(5,10),scale=2)
    wmin, l1 = get_max_path(enm, sensors_sub, effectors_sub)
#    source = random.choice([i for i in effector_pos.keys()])
#    target = random.choice([i for i in sensor_pos.keys()])
    #e_pcc.prs_mat_df=pd.DataFrame(e_pcc.prs_mat,columns=e_pcc.nodes, index=e_pcc.nodes)
#    l1 = e_pcc.get_prs_weighted_path(source,target)[1]
    #l1.extend(['fum1', 'fbp1'])
    sub = nx.induced_subgraph(enm.graph_gc, l1)
    path_init_pos=dict(zip(l1[1:-1],np.array(tuple(zip(np.repeat(5,(len(l1)-2)), np.linspace(effector_pos[l1[0]][1]-1,sensor_pos[l1[-1]][1],len(l1))[1:-1])))))
    
    for i in sub.nodes:
        if i in effector_pos.keys():
            path_init_pos[i]=effector_pos[i]
        if i in sensor_pos.keys():
            path_init_pos[i]=sensor_pos[i]
    return sensor_pos, effector_pos, path_init_pos,sub, l1, wmin