import networkx as nx
import numpy as np


def plot_network_spring(Gc, figure_path, **kwargs):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    
    spring_pos=Gc.nodes(data='pos')
    node_color=kwargs.pop('node_color','white')
    edge_color=kwargs.pop('edge_color','white')
    legend_elements = [Line2D([0], [0], marker='o', color=node_color, label=kwargs.pop('node_label','Genes'),
                            markerfacecolor=node_color, markersize=10,linestyle="None"),
        Line2D([0], [0], marker='o', color=edge_color, label=kwargs.pop('edge_label','PCC>0.2'),
                            markerfacecolor=edge_color, markersize=0,linestyle="-")
    ]
    fig, ax = plt.subplots(figsize = kwargs.pop('figsize',(5,5)))#figsize=(5,5))
    nx.draw_networkx_nodes(Gc,
                            node_size=kwargs.pop('node_size',0.2),
                            #alpha=0.5,
                            node_color=node_color,
                            pos=spring_pos,
                            label='Genes',
                        #node_shape=matplotlib.markers.MarkerStyle(marker='o',fillstyle='full')
                            )
    nx.draw_networkx_edges(Gc,
            alpha=kwargs.pop('edge_alpha',0.2),
            width=kwargs.pop('edge_width',0.1),
            edge_color=edge_color,
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
    fig.savefig(f"{figure_path}/{kwargs.pop('figure_name','plot')}.{kwargs.pop('figure_extension','png')}") 
    return ax 
    # plt.close()