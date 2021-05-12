import numpy as np
import pandas as pd
input_file = snakemake.input[0]#"data/raw/Data File S3. Genetic interaction profile similarity matrices/cc_ALL.txt"
output_file=snakemake.output[0]

df = pd.read_csv(input_file,'\t',index_col=[0,1], header=[0,1])

df_long = df.melt(ignore_index=False, col_level=0)

df_long_renamed = df_long.reset_index(inplace=False).rename(columns = {'level_0':'gene1','level_1':'gene1_sys', 'variable':'gene2','value':'pcc'})

df_long_renamed.loc[df_long_renamed.pcc>=0.2,['gene1','gene2','pcc']].to_csv(output_file, index=False,float_format="%.6f")
#df = df.rename(columns = {'index':'new column name'})
