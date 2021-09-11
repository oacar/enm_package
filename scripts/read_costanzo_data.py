import numpy as np
import pandas as pd
input_file = snakemake.input[0]#"data/raw/Data File S3. Genetic interaction profile similarity matrices/cc_ALL.txt"
output_file=snakemake.output[0]
thr = float(snakemake.params['threshold'])

df = pd.read_csv(input_file,'\t',index_col=[0,1], header=[0,1])

df_long = df.melt(ignore_index=False, col_level=0)

df_long_renamed = df_long.reset_index(inplace=False).rename(columns = {'level_0':'gene1','level_1':'Systematic gene name', 'variable':'gene2','value':'pcc'})

df_long_renamed.loc[df_long_renamed.pcc>=thr,['gene1','gene2','pcc']].to_csv(output_file, index=False,float_format="%.6f")
#df = df.rename(columns = {'index':'new column name'})

df_long_renamed.iloc[:,[0,1]].drop_duplicates().to_csv(snakemake.output[1],index=False)

sgd = pd.read_csv(snakemake.input[1],'\t',header=None)

#sgd.iloc[:,[0,3]]
#Save background list for later GO enrichment analysis
sgd.loc[sgd.iloc[:,3].isin(df_long_renamed['Systematic gene name']),0].to_csv(snakemake.output[2],index=False, header=None)
