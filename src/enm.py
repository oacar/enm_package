import os 
import networkx as nx
import pandas as pd 
import prody
import numpy as np
import igraph

from .visualize import *
#from .utilities import *

class Enm():
    """This object is wrapper around prody GNM object and networkx analysis
    """
    def __init__(self, name):
        """Constructor

        :param name: name of the object for future reference, can be arbitrary
        :type name: string
        """
        self.name = name
        self.figure_path = f"../reports/figures/"
        self.rewired = False
    def print_name(self):
        """This prints name variable
        """
        print(self.name)

    def spring_pos(self):
        """Create spring layout of the giant component network and assign positions to nodes
        """
        try:
            pos = nx.spring_layout(self.graph_gc,k=0.6,scale=4,iterations=200)
            nx.set_node_attributes(self.graph_gc,pos,'pos')
        except AttributeError:
            raise('Giant component is not set yet. First call read network or giant component')

    def read_network(self,path,**kwargs):
        """Read network file and assign it to object

        :param path: Network file path
        :type path: string
        
        """
        sep = kwargs.pop('sep',None)
        _, ext = os.path.splitext(path)
        #return fname, ext
        if ext =='.csv' or sep ==',':
            nw_read = pd.read_csv(path)
            # Create network ##########
            G = nx.from_pandas_edgelist(nw_read, source='gene1', target='gene2')
        elif ext == '.gpickle':
            G = nx.read_gpickle(path)
        elif ext =='.tsv' or sep =='\t':
            nw_read = pd.read_csv(path,sep='\t')
            # Create network ##########
            G = nx.from_pandas_edgelist(nw_read, source='gene1', target='gene2')
        elif ext=='.gml':
            G = nx.read_gml(path)
        self.G= G
        self.giant_component()

    def giant_component(self):
        """From the graph variable create giant component and assing nodes and degree to Enm object
        """

        Gc = max([self.G.subgraph(c).copy() for c in nx.connected_components(self.G)], key=len)
        self.graph_gc = Gc
        nodes = [_ for _ in self.graph_gc.nodes()]
        self.nodes= nodes
        degree = [deg for id, deg in list(self.graph_gc.degree)]
        self.degree = degree


    def laplacian_matrix(self):
        """get Laplacian matrix of the giant component. Wrapper around networkx.laplacian_matrix
        """
        self.L=nx.laplacian_matrix(self.graph_gc, weight=None).todense()

    def get_gnm(self):
        """Calculate GNM modes and collectivity values
        """
        if hasattr(self, 'L') is False:
            self.laplacian_matrix()
        gnm = prody.GNM()
        gnm.setKirchhoff(self.L)

        gnm.calcModes(n_modes=None, zeros=False)
        #gnm.nodes = nodes
        #gnm.degrees = degree
        self.gnm=gnm
        #showCrossCorr(gnm)
        #sqf_orig = prody.calcSqFlucts(self.gnm)
        coll = prody.calcCollectivity(self.gnm)
        self.coll=coll
        coll_index_sorted = sorted(range(len(self.coll)), key=lambda k: self.coll[k], reverse=True)

        self.coll_index_sorted=coll_index_sorted

    def get_prs(self):
        """Calculate Perturbation response matrix
        """
        try:
            prs_mat, _, _ = prody.calcPerturbResponse(self.gnm, suppress_diag=True)
        except Exception as e:
            raise(AttributeError('GNM is not calculated yet. Call get_gnm() first'))
        self.prs_mat = prs_mat

    def create_df(self):
        """Create an overall dataframe to store network related and GNM related values
        """
        df = pd.DataFrame()
        df['orf_name']=self.nodes
        df['deg']=self.degree
        eff_orig = np.sum(self.prs_mat, axis=1)
        sens_orig = np.sum(self.prs_mat, axis=0)
        eigvecs_df = pd.DataFrame(self.gnm.getEigvecs()[:,self.coll_index_sorted[:8]],columns=[f'eig_{i}' for i in range(8)])

        df_ = pd.merge(df,eigvecs_df,left_index=True,right_index=True)
        df_['eff'] = eff_orig
        df_['sens'] = sens_orig
        eigenvector_centr = nx.eigenvector_centrality_numpy(self.graph_gc)
        closeness_centr = nx.closeness_centrality(self.graph_gc)
        df_['btw']= betweenness_nx(self.graph_gc,normalized=True)
        df_['trans']= list((nx.clustering(self.graph_gc)).values())
        df_['eigenvec_centr'] = [eigenvector_centr[i] for i in eigenvector_centr]
        df_['closeness_centr'] = [closeness_centr[i] for i in closeness_centr]
        
        self.df=df_


    def gnm_analysis(self):
        """Wrapper to run gnm, prs and create_df
        """
        self.get_gnm()
        self.get_prs()
        self.create_df()
        

    def plot_network_spring(self,**kwargs):
        """Plot network with spring layout
        """
        Gc = self.graph_gc
        if nx.get_node_attributes(Gc,'pos')=={}:
            self.spring_pos()
        
        plot_network_spring(Gc, self.figure_path, **kwargs)

    def plot_collectivity(self,**kwargs):
        """Plot collectivity values and color most collective modes
        """
        plot_collectivity(self.coll,self.coll_index_sorted,self.figure_path, **kwargs)

    def plot_scatter(self,x,y,**kwargs):
        """Plot scatter plot of 2 columns in df

        :param x: first variable name
        :type x: string
        :param y: second variable name
        :type y: string
        """
        plot_scatter(self.df,x,y, self.figure_path, **kwargs)

    def simulate_rewire(self, rewired=False, rewire_df_name=None,arr_name=None, **kwargs):
        """Wrapper around enm.simulate_rewire function out of the class
        """
        arr, rewire_df = simulate_rewire(self.graph_gc,rewired, rewire_df_name,arr_name, **kwargs)
        self.arr = arr
        self.rewire_df = rewire_df

    def plot_correlation_density(self,x,y, **kwargs):
        """Plot correlation density between x and y. Requires rewire_df

        :param x: first variable name
        :type x: string
        :param y: second variable name
        :type y: string
        """
        try:
            plot_correlation_density(self.df,self.rewire_df,x,y,self.figure_path, **kwargs)
        except Exception as e:
            raise(AttributeError('rewire_df is not calculated yet. Call simulate_rewire() first'))
        

    def heatmap_annotated(self,**kwargs):
        """Plot PRS heatmap with clustering dendrograms
        """
        heatmap_annotated(self.prs_mat,self.figure_path)



def betweenness_ig(g, normalized=False):
    n = g.vcount()
    btw = g.betweenness()

    if normalized:
        norm = [2 * a / (n * n - 3 * n + 2) for a in btw]
        return norm
    else:
        return btw


def betweenness_nx(g, normalized=False):
    adj = nx.adjacency_matrix(g).todense()
    g_ig = igraph.Graph.Adjacency((adj > 0).tolist())
    return betweenness_ig(g_ig, normalized)


def rewire_network(Gc, **kwargs):
    """This is wrapper around networkx's rewire functions for my purposes. It can rewire by keeping degree, or if not it will generate a network with same number of nodes and arbitrary number of edges using barabasi_albert_graph or erdos_renyi_graph

    :param Gc: Network to be rewired    
    :type Gc: networkx object
    :return: Rewired network
    :rtype: networkx object
    """
    nodegseq = kwargs.pop('nodegseq',False)
    if nodegseq:
        random_network_type = kwargs.pop('random_network_type','ba')
        if random_network_type == 'ba':
            Gc_rewired= nx.barabasi_albert_graph(n=len(Gc.nodes),m=7)
        elif random_network_type == 'er':
            Gc_rewired = nx.erdos_renyi_graph(100,0.1)
    else:
        Gc_rewired = Gc.copy()
        swp_count = nx.connected_double_edge_swap(Gc_rewired, 2 * len(Gc_rewired.nodes))
    return Gc_rewired


#def rewire_network():

        

def simulate_rewire(Gc,rewired=False, rewire_df_name=None,arr_name=None, **kwargs):
    """This function reads rewired network GNM data or calls rewire function

    :param Gc: network to be rewired
    :type Gc: networkx object
    :param rewired: Is there a rewired GNM data available, defaults to False
    :type rewired: bool, optional
    :param rewire_df_name: If rewired is True, related dataframe path should be given, defaults to None
    :type rewire_df_name: string, optional
    :param arr_name: If rewired is True, related numpy array txt file path should be given, defaults to None
    :type arr_name: string, optional
    :raises ValueError: If rewired is True and rewire_df_name or arr_name is not given, raises error
    :return: arr
    :rtype: numpy array
    :return: rewire_df
    :rtype: pandas dataframe    
    """ 
    from scipy.stats import pearsonr, spearmanr
    from tqdm import tqdm

    if rewired:
            # rewire_df_name = kwargs.pop('rewire_df_name',None)
            # arr_name = kwargs.pop('arr_name',None)
            if rewire_df_name is None or arr_name is None:
                raise ValueError('Rewired dataframe path or arr_name is not given ')
            rewire_df = pd.read_csv(rewire_df_name)
            arr = np.loadtxt(arr_name) 
    else:
        

        sim_num=kwargs.pop('sim_num',10)
        arr=np.zeros((len(Gc.nodes),len(Gc.nodes),int(sim_num)))
        # pearsonr(eff,degree)[0]
        rewire_df = pd.DataFrame(columns=['eff_deg_pearson', 'sens_deg_pearson', 'eff_deg_spearman', 'sens_deg_spearman',
                                        'eff_btw_pearson', 'sens_btw_pearson', 'eff_btw_spearman', 'sens_btw_spearman',
                                        'eff_trans_pearson', 'sens_trans_pearson', 'eff_trans_spearman',
                                        'sens_trans_spearman'])
        eff_corr_list = []
        sens_corr_list = []
        eff_hist_list = []
        sens_hist_list = []
        sqf_corr_list_p = []
        sqf_corr_list_s = []


        for i in tqdm(range(int(sim_num))):
            try:
                Gc_rew = rewire_network(Gc,**kwargs)
                enm_rew = Enm('rewired')

                enm_rew.graph = Gc_rew
                enm_rew.giant_component()
                enm_rew.gnm_analysis()
                res = enm_rew.prs_mat
                #Gc_rew = enm_rew.graph_gc
                degree = enm_rew.degree
            except Exception as e:
                pass

            arr[:,:,i]=res
            eff_rew = enm_rew.df.eff.values
            sens_rew = enm_rew.df.sens.values
            # eff_hist_list.append(eff_rew)
            # sens_hist_list.append(eff_rew)
            betweenness = enm_rew.df.btw.values#betweenness_nx(Gc_rew, normalized=True)
            clustering_coeff = enm_rew.df.eff_trans_pearson.values
            rewire_df_itr = pd.DataFrame([[pearsonr(eff_rew, degree)[0], pearsonr(sens_rew, degree)[0],
                                        spearmanr(eff_rew, degree)[0], spearmanr(sens_rew, degree)[0],
                                        pearsonr(eff_rew, betweenness)[0], pearsonr(sens_rew, betweenness)[0],
                                        spearmanr(eff_rew, betweenness)[0], spearmanr(sens_rew, betweenness)[0],
                                        pearsonr(eff_rew, clustering_coeff)[0], pearsonr(sens_rew, clustering_coeff)[0],
                                        spearmanr(eff_rew, clustering_coeff)[0],
                                        spearmanr(sens_rew, clustering_coeff)[0]]],
                                        columns=['eff_deg_pearson', 'sens_deg_pearson', 'eff_deg_spearman',
                                                'sens_deg_spearman',
                                                'eff_btw_pearson', 'sens_btw_pearson', 'eff_btw_spearman',
                                                'sens_btw_spearman',
                                                'eff_trans_pearson', 'sens_trans_pearson', 'eff_trans_spearman',
                                                'sens_trans_spearman'])
            rewire_df = rewire_df.append(rewire_df_itr)
            # if i % 1 == 0:
            #     print(i)
        #rewire_df.to_csv(outpath+'/rewire_df.csv')

    return arr, rewire_df
