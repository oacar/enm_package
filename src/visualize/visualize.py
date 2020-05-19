import copy
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
import scipy.cluster.hierarchy as sch
from scipy.spatial import distance
from scipy.stats import pearsonr, spearmanr
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

def plot_vector(data, figure_path, sorted=False, **kwargs):
    mode_colors = kwargs.pop('mode_colors', [
                             'orange', 'blue', 'lightblue', 'tab:brown', 'darkgreen', 'm', 'crimson'])

    fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (5, 5)))
    if sorted:
        ax.plot(np.sort(data), '-o',
                c=mode_colors[kwargs.pop('color_id', 0)], zorder=1)
    else:
        ax.plot(data, '-o', c=mode_colors[kwargs.pop('color_id', 0)])
    ax.set_xlabel(kwargs.pop('x_label', 'Nodes'))
    ax.set_ylabel(kwargs.pop('y_label', 'Eigenvector'))
    #plt.title('Most collective mode')
    ax.axhline(0, linestyle='--', c='r')

    ax.autoscale(False)
    if sorted:
        # (range(len(data)), key=lambda k: np.abs(data)[k], reverse=False)
        index_sorted = np.argsort(np.abs(np.sort(data)))
    #    print(index_sorted[:round(len(data)*0.05)])
        ax.fill_between(index_sorted[:round(len(data)*0.05)], ax.get_ylim()[0], ax.get_ylim()[1],  # where=y > threshold,
                        color='grey', alpha=0.5)  # , transform=ax.get_xaxis_transform())

    fig.savefig(
        f"{figure_path}/{kwargs.pop('figure_name','eigenplot')}.{kwargs.pop('figure_extension','png')}")

    return ax


def plot_network_spring(Gc, figure_path, **kwargs):
    """Plots networkx object using spring_layout and a legend for nodes and edges

    :param Gc:  The network to plot
    :type Gc: networkx object
    :param figure_path: Folder to save plotted figure
    :type figure_path: string
    :return: returns Axes for downstream pipeline
    :rtype: Axes
    """

    spring_pos = Gc.nodes(data='pos')
    node_color = kwargs.pop('node_color', 'white')
    edge_color = kwargs.pop('edge_color', 'white')
    legend_elements = [Line2D([0], [0], marker='o', color=node_color, label=kwargs.pop('node_label', 'Genes'),
                              markerfacecolor=node_color, markersize=10, linestyle="None"),
                       Line2D([0], [0], marker='o', color=edge_color, label=kwargs.pop('edge_label', 'PCC>0.2'),
                              markerfacecolor=edge_color, markersize=0, linestyle="-")
                       ]
    fig, ax = plt.subplots(figsize=kwargs.pop(
        'figsize', (5, 5)))  # figsize=(5,5))
    nx.draw_networkx_nodes(Gc,
                           node_size=kwargs.pop('node_size', 0.2),
                           # alpha=0.5,
                           node_color=node_color,
                           pos=spring_pos,
                           label='Genes',
                           # node_shape=matplotlib.markers.MarkerStyle(marker='o',fillstyle='full')
                           )
    nx.draw_networkx_edges(Gc,
                           alpha=kwargs.pop('edge_alpha', 0.2),
                           width=kwargs.pop('edge_width', 0.1),
                           edge_color=edge_color,
                           pos=spring_pos,
                           label='PCC>0.2')
    ax.set_facecolor(kwargs.pop('facecolor', "#000000"))
    # plot_go_contours(Gc,ax,go_df_list[i],1,clabels=True,level=0.01)
    plot_legend = kwargs.pop('plot_legend', False)
    if plot_legend:
        plt.legend(handles=legend_elements, fontsize=14, loc='lower right')
    #
    #plt.title(f'Costanzo 2016 profile similarity network',fontsize=20)
    # plt.legend()
    fig.savefig(
        f"{figure_path}/{kwargs.pop('figure_name','network_plot')}.{kwargs.pop('figure_extension','png')}")
    return ax
    # plt.close()


def plot_collectivity(coll, coll_index_sorted, figure_path, **kwargs):
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
    fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (4, 4)))
    plt.scatter(coll_index_sorted[:7], [
                coll[i] for i in coll_index_sorted[:7]], s=20, color=mode_colors, zorder=2)
    plt.scatter(0, coll[0], s=20, c='darkred', zorder=2)
    plt.plot(range(len(coll)), coll, c='tab:cyan', linewidth=0.5, zorder=1)

    # plt.axhline(np.quantile(coll,0.99),c='r')
    #plt.title('GNM modes collectivity')
    plt.xlabel(kwargs.pop('xlabel', 'Modes'))
    plt.ylabel(kwargs.pop('ylabel', 'Collectivity'))
    plt.savefig(
        f"{figure_path}/{kwargs.pop('figure_name','coll')}.{kwargs.pop('figure_extension','png')}", transparent=True)
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
    fig, ax = plt.subplots()
    #plt.hist([x for x,y in deg_sens_rew], label='Randomized')
    #sns.distplot([x for x,y in deg_sens_rew],label='Randomized', kde=False)

    # Creating another Y axis
    second_ax = ax.twinx()

    # Plotting kde without hist on the second Y axis
    #sns.distplot([x for x,y in deg_sens_rew],label='Randomized', ax=second_ax, kde=True, hist=False)
    # sns.distplot(rewire_df_nodegseq.sens_deg_spearman,label='Randomized2',hist=False)
    correlation = kwargs.pop('correlation', 'pearson')
    sns.distplot(rewire_df[f'{x}_{y}_{correlation}'],
                 label='Randomized', ax=second_ax, kde=True, hist=False)

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
    ax.yaxis.set_label_text('Counts')
    ax.xaxis.set_label_text(f"{np.where(x=='eff','Effectiveness','Sensitivity')} {np.where(y=='deg','Degree',np.where(y=='btw', 'Betweenness','Transitivity'))} correlation")

    ax.yaxis.set_label_position("left")
    plt.legend()
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    figure_name = kwargs.pop('figure_name',f'correlation_density_{x}_{y}')
    fig.savefig(f"{figure_path}/{figure_name}.{kwargs.pop('figure_extension','png')}", transparent=True)
    return ax


def heatmap_annotated(prs_mat, figure_path, **kwargs):

    # sch.set_link_color_palette(['b'])
    
    # Dendrogram that comes to the left
    fig = plt.figure(figsize=kwargs.pop('figsize', (8, 8)))
    # Add an axes at position rect [left, bottom, width, height]
    quantile_threshold = kwargs.pop('quantile_threshold', 0.95)
    method = kwargs.pop('method', 'ward')
    q99 = np.quantile(prs_mat, quantile_threshold)
    prs_mat_cl = copy.deepcopy(prs_mat)
    prs_mat_cl[prs_mat_cl > q99] = q99
    # , optimal_ordering=True)
    row_linkage = sch.linkage(distance.pdist(prs_mat_cl), method=method)
    # ,optimal_ordering=True)
    col_linkage = sch.linkage(distance.pdist(prs_mat_cl.T), method=method)
    ax1 = fig.add_axes([0.09, 0.1, 0.2, 0.6])
   # Y = sch.linkage(D, method='centroid')
    # orientation='left' is reponsible for making the
    # dendrogram appear to the left
    Z1 = sch.dendrogram(row_linkage, orientation='left',
                        link_color_func=lambda k: 'black')
    ax1.set_xticks([])
    ax1.set_yticks([])
    plt.axis('off')
    # top side dendogram
    ax2 = fig.add_axes([0.3, 0.71, 0.6, 0.2])
    #Y = sch.linkage(D, method='single')
    Z2 = sch.dendrogram(col_linkage, color_threshold=0,
                        link_color_func=lambda k: 'black')
    ax2.set_xticks([])
    ax2.set_yticks([])
    plt.axis('off')
    # main heat-map
    axmatrix = fig.add_axes([0.3, 0.1, 0.6, 0.6])
    idx1 = Z1['leaves']
    idx2 = Z2['leaves']
    prs_mat = prs_mat[idx1, :]
    prs_mat = prs_mat[:, idx2]
    prs_mat_cl = prs_mat_cl[idx1, :]
    prs_mat_cl = prs_mat_cl[:, idx2]
    # the actual heat-map
    im = axmatrix.matshow(prs_mat_cl, aspect='auto',
                          origin='lower', cmap="YlGnBu")
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # xticks to the right (x-axis)
    axmatrix.set_xticks(range(40))
    axmatrix.set_xticklabels(idx1, minor=False)
    axmatrix.xaxis.set_label_position('bottom')
    axmatrix.xaxis.tick_bottom()

    plt.xticks(rotation=-90, fontsize=8)  # ,colors='black')

    # xticks to the right (y-axis)
    axmatrix.set_yticks(range(40))
    axmatrix.set_yticklabels(idx2, minor=False)
    axmatrix.yaxis.set_label_position('right')
    axmatrix.yaxis.tick_right()

    # to add the color bar
    # axcolor = fig.add_axes([0.94, 0.1, 0.02, 0.6])
    plt.axis('off')
    ax_rowdata = fig.add_axes([0.9, 0.09, 0.1, 0.63])
    row_data = np.mean(prs_mat, axis=1)
    # ,orientation=u'vertical')
    ax_rowdata.plot((row_data), range(len(row_data)), '-')
    plt.axis('off')

    ax_coldata = fig.add_axes([0.28, -0.0, 0.64, 0.1])
    col_data = np.mean(prs_mat, axis=0)
    # ,orientation=u'vertical')
    ax_coldata.plot(range(len(col_data)), col_data, '-')
    plt.axis('off')
    ax_colorbar = fig.add_axes([0.0, 0.71, 0.1, .2])
    plt.colorbar(im, ax=ax_colorbar)
    plt.axis('off')
    # plt.show()
    outname = f"{figure_path}/{kwargs.pop('figure_name','correlation_density_{x}_{y}')}.{kwargs.pop('figure_extension','png')}"
    if outname is not None:
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

    # return ax
