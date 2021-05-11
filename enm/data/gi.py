import pandas as pd

gi = pd.read_csv('~/pgsNetwork/analysis/data/derived_data/SGA_data_combined.csv',sep='\t')
gi_strong = gi.loc[(gi['Genetic interaction score (ε)']< -.08) & (gi['P-value']<0.05)]
gi_strong['gene1'] = gi_strong['Query allele name']
gi_strong['gene2'] = gi_strong['Array allele name']

gi_verystrong =gi_strong.loc[(gi['Genetic interaction score (ε)']< -.2) & (gi['P-value']<0.05)]
 
gi_strong.to_csv('data/interim/gi/gi_.08.csv')
gi_verystrong.to_csv('data/interim/gi/gi_.2.csv')
