import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
import scipy.cluster.hierarchy as sch
from scipy.stats import pearsonr, spearmanr, gaussian_kde
figsize = (5, 5)
#matplotlib.rcParams['pdf.fonttype'] = 42
#matplotlib.rcParams['ps.fonttype'] = 42
params = {'legend.fontsize': 'x-large',
          'figure.figsize': figsize,
          'axes.labelsize': 'x-large',
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 'x-large',
          'ytick.labelsize': 'x-large',
          'pdf.fonttype': 42,
          'ps.fonttype': 42}
plt.rcParams.update(params)


# TODO Add eigenvector plotting functions, with shading capability
# TODO Add GO kernel plotting function

def plot_vector(data, figure_path, sorted=False, x_label='Nodes', y_label='Eigenvector',**kwargs):
    mode_colors = kwargs.pop('mode_colors', [
                             'orange', 'blue', 'lightblue', 'tab:brown', 'darkgreen', 'm', 'crimson','darkred'])

    fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (5, 5)))
    color_id = kwargs.pop('color_id', 0)

    if sorted:
        ax.plot(np.sort(data), '-o',
                c=mode_colors[color_id], zorder=1)
    else:
        ax.plot(data, '-o', c=mode_colors[color_id])
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    #plt.title('Most collective mode')
    ax.axhline(0, linestyle='--', c='r')

    ax.autoscale(False)
    if sorted:
        # (range(len(data)), key=lambda k: np.abs(data)[k], reverse=False)
        index_sorted = np.argsort(np.abs(np.sort(data)))
    #    print(index_sorted[:round(len(data)*0.05)])
        if color_id == -1:
            ax.fill_between(index_sorted[:round(len(data)*0.05)], ax.get_ylim()[0], ax.get_ylim()[1],  # where=y > threshold,
                        color='grey', alpha=0.5)  # , transform=ax.get_xaxis_transform())
        else:
            ax.fill_between(index_sorted[:round(len(data)*0.05)], ax.get_ylim()[0], ax.get_ylim()[1],  # where=y > threshold,
                        color='grey', alpha=0.5)  # , transform=ax.get_xaxis_transform())


    plt.tight_layout()
    fig.savefig(
        f"{figure_path}/{kwargs.pop('figure_name','eigenplot')}.{kwargs.pop('figure_extension','png')}")

    return ax


def plot_network_spring(Gc, figure_path, plot_go=False, go_df_list=None, level_list=[0.1],ax=None, savefig=False,**kwargs):
    """Plots networkx object using spring_layout and a legend for nodes and edges

    :param Gc:  The network to plot
    :type Gc: networkx object
    :param figure_path: Folder to save plotted figure
    :type figure_path: string
    :return: returns Axes for downstream pipeline
    :rtype: Axes
    """
    mode_colors = kwargs.pop('mode_colors', [
                             'orange', 'blue', 'lightblue', 'tab:brown', 'darkgreen', 'm', 'crimson'])


    if plot_go and go_df_list is None:
        plot_go=False
        print('GO dataframe list is not given with kw=go_df_list. GO contours are not plotted')

    spring_pos = Gc.nodes(data='pos')
    node_color = kwargs.pop('node_color', 'white')
    edge_color = kwargs.pop('edge_color', 'white')
    legend_elements = kwargs.pop('legend_elements',None)
    plot_legend = kwargs.pop('plot_legend',False)
    if legend_elements is None and plot_legend:
        print('ss')
        legend_elements = [Line2D([0], [0], marker='o', color=node_color, label=kwargs.pop('node_label', 'Genes'),
                              markerfacecolor=node_color, markersize=10, linestyle="None"),
                       Line2D([0], [0], marker='o', color=edge_color, label=kwargs.pop('edge_label', 'PCC>0.2'),
                              markerfacecolor=edge_color, markersize=0, linestyle="-")
                       ]

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop(
         'figsize', (5, 5)))  # figsize=(5,5))
    nx.draw_networkx_nodes(Gc,
                           node_size=kwargs.pop('node_size', 0.2),
                           # alpha=0.5,
                           node_color=node_color,
                           pos=spring_pos,
                           label='Genes', ax=ax, **kwargs
                           # node_shape=matplotlib.markers.MarkerStyle(marker='o',fillstyle='full')
                           )
    nx.draw_networkx_edges(Gc,
                           alpha=kwargs.pop('edge_alpha', 0.2),
                           width=kwargs.pop('edge_width', 0.1),
                           edge_color=edge_color,
                           pos=spring_pos,
                           label='PCC>0.2',ax=ax, **kwargs)
    ax.set_facecolor(kwargs.pop('facecolor', "#000000"))

    if plot_go:
        for i,go_df in enumerate(go_df_list):
            plot_go_contours(Gc,ax,go_df,1,color=mode_colors[i],clabels=False,level=level_list[i])
            legend_elements.append(
            Line2D([0], [0], marker='o', color=mode_colors[i], label=f'{go_df.name[0]}',
                           markersize=0,linestyle="-")
                    )
    # plot_go_contours(Gc,ax,go_df_list[i],1,clabels=True,level=0.01)
    plot_legend = kwargs.pop('plot_legend', False)
    if savefig:
        if plot_legend:
            lgd = ax.legend(handles=legend_elements, fontsize=14,loc='center left', bbox_to_anchor=(1.0, 0.5))
        
     #       frame = lgd.get_frame()
    #        frame.set_color('black')
      #      for text in lgd.get_texts():
       #         text.set_color("white") 
            #plt.title(f'Costanzo 2016 profile similarity network',fontsize=20)
        # plt.legend()
            plt.savefig(
                f"{figure_path}/{kwargs.pop('figure_name','network_plot')}.{kwargs.pop('figure_extension','png')}",bbox_extra_artists=(lgd,),bbox_inches='tight',**kwargs)
        else:
            plt.savefig(
                f"{figure_path}/{kwargs.pop('figure_name','network_plot')}.{kwargs.pop('figure_extension','png')}",bbox_inches='tight',**kwargs)

    return ax
    # plt.close()


def plot_collectivity(coll, coll_index_sorted, figure_path, x_label='Modes', y_label='Collectivity', **kwargs):
    """[summary]

    :param coll: list of collectivity values calculated by prody.calcCollectivity
    :type coll: list
    :param coll_index_sorted: Indexes of collectivity values sorted in descending order
    :type coll_index_sorted: list
    :param figure_path: Folder to save figure
    :type figure_path: string
    :return: returns Axes object for downstream pipeline
    :rtype: Axes
    """

    mode_colors = kwargs.pop('mode_colors', [
                             'orange', 'blue', 'lightblue', 'tab:brown', 'darkgreen', 'm', 'crimson'])
    fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (5, 5)))
    plt.scatter(coll_index_sorted[0], [
        coll[i] for i in [coll_index_sorted[0]]], s=20, color=mode_colors[0], zorder=2)
    plt.scatter(0, coll[0], s=20, c='darkred', zorder=2)
    plt.plot(range(len(coll)), coll, c='tab:cyan', linewidth=0.5, zorder=1)

    # plt.axhline(np.quantile(coll,0.99),c='r')
    #plt.title('GNM modes collectivity')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.tight_layout()
    plt.savefig(
        f"{figure_path}/{kwargs.pop('figure_name','coll')}.{kwargs.pop('figure_extension','png')}", transparent=False)
    return ax


def plot_scatter(df, x, y, figure_path, **kwargs):
    """This is a wrapper around data frame to make easier scatter plots. It uses seaborn and matplotlib

    :param df: a data frame to plot scatter for
    :type df: pandas data frame
    :param x: pandas data frame column for x axis
    :type x: string
    :param y: pandas data frame column for y axis
    :type y: string
    :param figure_path: folder to save figure. automatically saves figure
    :type figure_path: string
    :return: returns Axes object for downstream pipeline
    :rtype: Axes
    """

    hue = kwargs.pop('hue', None)
# TODO ============== Add categorical hue
    fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (5, 5)))
    sns.set_style("ticks")
    if hue is None:
        # ,hue='cat_4',hue_order=sorted(combined_df.cat_4.unique()))
        sns.scatterplot(x=x, y=y, data=df, edgecolor=None)
    else:
        sns.scatterplot(x=x, y=y, data=df, edgecolor=None,
                        hue=hue, hue_order=sorted(df[hue].unique()))
    # plt.ylim(0,0.005)
    # plt.xlim(0,)
    plt.ylabel(kwargs.pop('y_label', y))
    plt.xlabel(kwargs.pop('x_label', x))

    # sns.axes_style
    # plt.grid(False)
    # plt.axis('off')
    #plt.axhline(np.quantile(df_['eff'],0.99),c='r',label='99% Effectiveness')
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.tight_layout()
    fig.savefig(
        f"{figure_path}/{kwargs.pop('figure_name','scatter')}.{kwargs.pop('figure_extension','png')}", transparent=True)
    return ax


def plot_correlation_density(df, rewire_df, x, y, figure_path, **kwargs):
    """Plot density of correlation values between x and y. Needs rewire_df to be defined

    :param df: data frame to calculate correlations
    :type df: pandas dataframe
    :param rewire_df: dataframe for simulated values. needs to have correlation between x and y
    :type rewire_df: pandas dataframe
    :param x: column name for first variable, sens or eff
    :type x: string
    :param y: column name for second variable, deg, btw or trans
    :type y: string
    :param figure_path: folder to save figure. automatically saves figure
    :type figure_path: string
    :return: returns Axes object for downstream pipeline
    :rtype: Axes
    """
    figsize = kwargs.pop('figsize',(5,5))
    fig, ax = plt.subplots(figsize=figsize)
    #plt.hist([x for x,y in deg_sens_rew], label='Randomized')
    #sns.distplot([x for x,y in deg_sens_rew],label='Randomized', kde=False)

    # Creating another Y axis
    second_ax = ax.twinx()

    # Plotting kde without hist on the second Y axis
    #sns.distplot([x for x,y in deg_sens_rew],label='Randomized', ax=second_ax, kde=True, hist=False)
    # sns.distplot(rewire_df_nodegseq.sens_deg_spearman,label='Randomized2',hist=False)
    correlation = kwargs.pop('correlation', 'pearson')
    #print(correlation)
    labels = kwargs.pop("labels",['Keeping degree','Barabasi-Albert', 'Erdos-renyi'])

    for itr,t in enumerate(rewire_df):
        sns.distplot(t[f'{x}_{y}_{correlation}'], label=labels[itr],
                     ax=second_ax, kde=True, hist=False)
    # Removing Y ticks from the second axis
    second_ax.set_yticks([])
    # sns.distplot()
    color_real = kwargs.pop('color_real', 'r')
    label_real = kwargs.pop('label_real', 'Real')
    if correlation is 'pearson':
        plt.axvline(pearsonr(df[x], df[y])[0], c=color_real, label=label_real)
    else:
        plt.axvline(spearmanr(df[x], df[y])[0], c=color_real, label=label_real)
    # plt.ylabel('Counts')
    ax.yaxis.set_label_text('Density')
    ax.xaxis.set_label_text(f"{np.where(x=='eff','Effectiveness','Sensitivity')} {np.where(y=='deg','Degree',np.where(y=='btw', 'Betweenness','Transitivity'))} correlation")

    ax.yaxis.set_label_position("left")
    lgd = second_ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
#    plt.legend()
#    plt.tight_layout()
    figure_name = kwargs.pop('figure_name',f'correlation_density_{x}_{y}')
    fig.savefig(f"{figure_path}/{figure_name}.{kwargs.pop('figure_extension','png')}", transparent=False,bbox_extra_artists=(lgd,),bbox_inches='tight')
    return ax

def heatmap_annotated(prs_mat, prs_mat_cl_orig, figure_path, row_linkage, col_linkage, row_colors = None, col_colors=None, save_figure=False, **kwargs):
    from seaborn.matrix import ClusterGrid
    from seaborn import heatmap
    fig = plt.figure(figsize=kwargs.pop('figsize', (8, 8)))
    # Add an axes at position rect [left, bottom, width, height]
    ncol = 4 if row_colors is not None else 3
    width_ratios = [1,.2,4,.5] if row_colors is not None else [1,4,.5]
    nrow = 4 if col_colors is not None else 3
    height_ratios = [1,.2,4,.5] if col_colors is not None else [1,4,.5]
    from matplotlib import gridspec
    gs = gridspec.GridSpec(nrow, ncol, width_ratios = width_ratios, height_ratios = height_ratios)
    if row_colors is not None and col_colors is not None:
        ax_row_colors = fig.add_subplot(gs[2,1])
        ax_col_colors = fig.add_subplot(gs[1,2])

        ax_row_dend = fig.add_subplot(gs[2,0])
        ax_col_dend = fig.add_subplot(gs[0,2])

        ax_heatmap = fig.add_subplot(gs[2,2])
        ax_row_data = fig.add_subplot(gs[2,3])
        ax_col_data = fig.add_subplot(gs[3,2])
    else:
        ax_row_dend = fig.add_subplot(gs[1,0])
        ax_col_dend = fig.add_subplot(gs[0,1])

        ax_heatmap = fig.add_subplot(gs[1,1])
        ax_row_data = fig.add_subplot(gs[1,2])
        ax_col_data = fig.add_subplot(gs[2,1])


    ax_row_data.set_axis_off()
    ax_col_data.set_axis_off()
    ax_row_dend.set_axis_off()
    ax_col_dend.set_axis_off()
    # orientation='left' is reponsible for making the
    # dendrogram appear to the left
    Z1_cl = sch.dendrogram(row_linkage, orientation='left',
            link_color_func=lambda k: 'black', ax = ax_row_dend)
    # top side dendogram
    #Y = sch.linkage(D, method='single')
    Z2_cl = sch.dendrogram(col_linkage, color_threshold=0, 
                        link_color_func=lambda k: 'black',
                    ax = ax_col_dend)


    #axmatrix = fig.add_axes([0.3, 0.1, 0.6, 0.6])
    idx1_cl = Z1_cl['leaves']
    idx2_cl = Z2_cl['leaves']
    if row_colors is not None and col_colors is not None:
        matrix, cmap = ClusterGrid.color_list_to_matrix_and_cmap(row_colors, idx1_cl,axis=0)
        heatmap(np.flip(matrix), cmap=cmap, cbar=False, ax = ax_row_colors, xticklabels=False, yticklabels=False)

        matrix, cmap = ClusterGrid.color_list_to_matrix_and_cmap(col_colors, idx2_cl,axis=1)
        heatmap(matrix, cmap=cmap, cbar=False, ax = ax_col_colors, xticklabels=False, yticklabels=False)
    #prs_mat = e_pcc.prs_mat
    #prs_mat = prs_mat[idx1, :]
    #prs_mat = prs_mat[:, idx2]
    prs_mat_cl = prs_mat_cl_orig[idx1_cl, :]
    prs_mat_cl = prs_mat_cl[:, idx2_cl]
    # the actual heat-map
    im = ax_heatmap.matshow(prs_mat_cl, aspect='auto', #norm=mpl.colors.PowerNorm(gamma=2), 
                        origin='lower', cmap="YlGnBu")


    # xticks to the right (x-axis)
    #ax_heatmap.set_xticks(range(40))
    ax_heatmap.set_xticklabels(idx1_cl, minor=False)
    ax_heatmap.xaxis.set_label_position('bottom')
    ax_heatmap.xaxis.tick_bottom()

    plt.xticks(rotation=-90, fontsize=8)  # ,colors='black')

    # xticks to the right (y-axis)
    #ax_heatmap.set_yticks(range(40))
    ax_heatmap.set_yticklabels(idx2_cl, minor=False)
    ax_heatmap.yaxis.set_label_position('right')
    ax_heatmap.yaxis.tick_right()
    #ax_heatmap.set_axis_off()
    # to add the color bar
    # axcolor = fig.add_axes([0.94, 0.1, 0.02, 0.6])
    ax_colorbar = fig.add_subplot(gs[0, 0])
    ax_colorbar.set_axis_off()
    row_data = np.mean(prs_mat[idx1_cl,:], axis=1)#np.mean(prs_mat_cl, axis=1)
        # ,orientation=u'vertical')
    ax_row_data.plot((row_data), range(len(row_data)), '-')
    ax_row_data.set_ylim(0,len(row_data))

    col_data = np.mean(prs_mat[:,idx2_cl], axis=0)#np.mean(prs_mat_cl, axis=0)
        # ,orientation=u'vertical')
    ax_col_data.plot(range(len(col_data)), col_data, '-')
    ax_col_data.set_xlim(0,len(col_data))
    #plt.axis('off')
    #plt.axis('off')
    ax_heatmap.set_xticks([])
    ax_heatmap.set_yticks([])
    plt.colorbar(im, ax=ax_colorbar)
    
    # plt.show()
    outname = f"{figure_path}/{kwargs.pop('figure_name','prs_heatmap')}.{kwargs.pop('figure_extension','png')}"
    if save_figure:
        plt.savefig(outname)

        # plt.figure()
        # plt.plot(range(len(col_data)),col_data,'-')
        # plt.xticks([])
        # plt.savefig(outname+'_coldata.png',dpi=100)

        # plt.figure()
        # plt.plot(range(len(row_data)),np.flip(row_data),'-')
        # plt.xticks([])
        # plt.savefig(outname+'_rowdata.png',dpi=100)

    else:
        plt.show()

    return row_data, col_data
    # return ax

def plot_go_contours(Gc, ax,go_df, k =1,clabels=False,level=1e-6,pos=None,**kwargs):
    color_ = kwargs.pop('color','#00000F')
    if pos is None:
        pos = dict(Gc.nodes.data('pos'))
    
    #x, y = [x, y for x,y in Gc.nodes.data('pos')]
    min_pos = np.min([pos[key] for key in pos],axis=0)
    max_pos = np.max([pos[key] for key in pos],axis=0)
    labels = nx.get_node_attributes(Gc, 'orf_name')
    labels_dict = {k: v for v, k in labels.items()}
    for i in  range(1):#range(np.min([go_df.shape[0],5])):
        nodes = go_df.iloc[i,:].study_items.split(', ')
        if len(nodes)<5:
            nodes = go_df.iloc[i+1,:].study_items.split(', ')
           
        nodes_indices = [labels_dict[node] for node in nodes if node in labels_dict.keys()]
        X,Y,Z = create_contour(pos, nodes_indices, k, max_pos, min_pos)
        C = ax.contour(X, Y, Z, [level],colors=color_)#, colors=[tuple(process_colors[n_process, :])], alpha=1)
        if clabels:
            fmt = {}
            strs = [go_df.iloc[i,:]['name']]
            for l, s in zip(C.levels, strs):
                fmt[l] = s
    #        print(i)
            plt.clabel(C, C.levels, inline=False, fmt=fmt, fontsize=18,use_clabeltext=True)

def create_contour(pos, nodes_indices, k, max_pos, min_pos):
    pos3 = {idx: pos[node_index] for idx, node_index in enumerate(nodes_indices)}
#        print(pos3)
    pos3 = np.vstack(list(pos3.values()))
    #pos3 = remove_outliers(pos3,k)
    kernel = gaussian_kde(pos3.T)
    [X, Y] = np.mgrid[min_pos[0]:max_pos[0]:100j,min_pos[1]:max_pos[1]:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    Z = np.reshape(kernel(positions).T, X.shape)
    
    return X,Y,Z
def remove_outliers(arr, k):
    mu, sigma = np.mean(arr, axis=0), np.std(arr, axis=0, ddof=1)
    return arr[np.all(np.abs((arr - mu) / sigma) < k, axis=1)]

def plot_fiedler_data(enm,figure_name = "lost_edges_node_counts", figure_extension = 'png', **kwargs):
    figsize = kwargs.pop('figsize',(6,12))

    lost_edges=enm.lost_edges
    node_counts_c1 = enm.node_counts_c1
    fig,(ax1,ax2,ax3) = plt.subplots(3,2,figsize=figsize)

    ax1[0].plot(lost_edges,'o')
    ax1[0].set_xlabel('Modes')
    ax1[0].set_ylabel('Lost Edges')
    ax1[1].plot(node_counts_c1/len(enm.graph_gc.nodes),'o')
    ax1[1].set_xlabel('Modes')
    ax1[1].set_ylabel('Smaller Cluster size')
    #ax1[2].plot(node_counts_c1,lost_edges,'o')
    #ax1[2].set_xlabel('Smaller Cluster size')
    #ax1[2].set_ylabel('Lost Edges')

    ax2[0].plot(enm.coll,lost_edges,'o')
    ax2[0].set_xlabel('Collectivity')
    ax2[0].set_ylabel('Lost Edges')
    ax2[1].plot(enm.coll,node_counts_c1,'o')
    ax2[1].set_xlabel('Collectivity')
    ax2[1].set_ylabel('Smaller cluster size')


    ax3[0].plot(enm.gnm.getEigvals(),lost_edges,'o')
    ax3[0].set_xlabel('Lambda')
    ax3[0].set_ylabel('Lost Edges')

    ax3[1].plot(enm.gnm.getEigvals(),node_counts_c1,'o')
    ax3[1].set_xlabel('Lambda')
    ax3[1].set_ylabel('Smaller cluster size')


    fig.tight_layout()
    fig.savefig(f'{enm.figure_path}/{figure_name}.{figure_extension}',bbox_to_anchor='tight')
    plt.close()

def plot_lambda_collectivity(lambda_, coll_,figure_path,**kwargs):
    figure_name = kwargs.pop('figure_name','lambda_coll')
    figure_extension = kwargs.pop('figure_extension','png')
    fig,ax = plt.subplots(figsize=kwargs.pop('figsize',(7,7)))
    ax.scatter(lambda_, coll_, c='tab:cyan', linewidth=0.5, zorder=1)
    ax.set_xlabel('Lambda')
    ax.set_ylabel('Collectivity')
    plt.tight_layout()
    fig.savefig(f"{figure_path}/{figure_name}.{figure_extension}")
    plt.close()

