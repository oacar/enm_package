import os 
import networkx as nx
import pandas as pd 
import prody
import igraph
import numpy as np

class Enm():
    def __init__(self, name):
        self.name = name
        self.figure_path = f"../figures/"
    def print_name(self):
        print(self.name)

    def spring_pos(self):
        try:
            pos = nx.spring_layout(self.graph_gc,k=0.6,scale=4,iterations=200)
            nx.set_node_attributes(self.graph_gc,pos,'pos')
        except AttributeError:
            raise('Giant component is not set yet. First call read network or giant component')
    def read_network(self,path,**kwargs):
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
        Gc = max([self.G.subgraph(c).copy() for c in nx.connected_components(self.G)], key=len)
        self.graph_gc = Gc
        nodes = [_ for _ in self.graph_gc.nodes()]
        self.nodes= nodes
        degree = [deg for id, deg in list(self.graph_gc.degree)]
        self.degree = degree


    def laplacian_matrix(self):
        self.L=nx.laplacian_matrix(self.graph_gc, weight=None).todense()

    def get_gnm(self):
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
        prs_mat, _, _ = prody.calcPerturbResponse(self.gnm, suppress_diag=True)
        self.prs_mat = prs_mat

    def create_df(self):
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

        self.get_gnm()
        self.get_prs()
        self.create_df()
        
    def foo(self):
        return 'foo'#

    def plot_network_spring(self, **kwargs):
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        Gc = self.graph_gc
        spring_pos=Gc.nodes(data='pos')
        legend_elements = [Line2D([0], [0], marker='o', color='red', label='Genes',
                                markerfacecolor='red', markersize=10,linestyle="None"),
            Line2D([0], [0], marker='o', color='black', label='PCC>0.2',
                                markerfacecolor='black', markersize=0,linestyle="-")
        ]
        fig, ax = plt.subplots(figsize=(5,5))
        nx.draw_networkx_nodes(Gc, #node_color=['#00b4d9' if x in hng.iloc[:,0].tolist() else (1,0,0,0.01) for x in Gc.nodes],
                                #node_size =[100 if x in hng.iloc[:,0].tolist() else 10 for x in Gc.nodes],
                                #alpha=0.1,
                                node_size=kwargs.pop('node_size',0.2),
                                #alpha=0.5,
                                node_color=kwargs.pop('node_color','white'),
                                pos=spring_pos,
                                label='Genes',
                            #node_shape=matplotlib.markers.MarkerStyle(marker='o',fillstyle='full')
                                )
        nx.draw_networkx_edges(Gc,
                alpha=kwargs.pop('edge_alpha',0.2),
                width=kwargs.pop('edge_width',0.1),
                edge_color=kwargs.pop('edge_color','white'),
                pos=spring_pos,
                label='PCC>0.2')
        ax.set_facecolor(kwargs.pop('facecolor',"#000000"))
        #plot_go_contours(Gc,ax,go_df_list[i],1,clabels=True,level=0.01)
        plot_legend = kwargs.pop('plot_legend',False)
        if plot_legend:
            plt.legend(handles=legend_elements,fontsize=14,loc='lower right')
        #
        #plt.title(f'Costanzo 2016 profile similarity network',fontsize=20)
        #plt.legend()
        fig.savefig(f"{self.figure_path}/{kwargs.pop('figure_name','plot')}.{kwargs.pop('figure_extension','png')}") 
        return ax 
        # plt.close()
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