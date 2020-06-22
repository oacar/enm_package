#enm.plot_correlation_density(x='eff',y='deg',figure_extension='pdf')
#enm.plot_correlation_density(x='sens',y='deg',correlation='spearman',figure_extension='pdf')
#nx.set_node_attributes(enm.graph_gc,dict(zip(list(enm.graph_gc.nodes),list(enm.graph_gc.nodes))),'orf_name')

## GO Analysis for sensors
#go_df = pd.read_csv('data/interim/go_results/4.tsv','\t')
sensors = enm.df.loc[enm.df.loc[:,'sens'] > np.quantile(enm.df.loc[:,'sens'],0.99)].orf_name.values
strain_ids = pd.read_csv('data/interim/strain_ids_with_experiment_count_all.csv')
sensors_orf_name = strain_ids.loc[strain_ids['Allele Gene name'].isin(sensors)]['Systematic gene name'].unique()
go_path = 'data/raw/ontology/'
#goeaobj, geneid2name = create_goea()
goeaobj, goea_results_sig,geneid2name, go_df= go_findenrichment(query = sensors_orf_name, obo_fname=f'{go_path}/goslim_generic.obo', gaf = f'{go_path}/sgd.gaf',background=f'{go_path}/sgd_costanzogenes',sgd_info_tab=f'{go_path}/SGD_features.tab')
#enm.plot_network_spring(plot_go=True,go_df_list=[go_df],level_list=[0.2])
