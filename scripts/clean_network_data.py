import argparse
import pandas as pd
#input_file = snakemake.input[0]#"data/raw/Data File S3. Genetic interaction profile similarity matrices/cc_ALL.txt"
#output_file=snakemake.output[0]

def read_costanzo_data(input_file, sgd_file, background_output, nw_output, strain_ids_output, thr):
#    thr = float(snakemake.params['threshold'])
    df = pd.read_csv(input_file,'\t',index_col=[0,1], header=[0,1])

    df_long = df.melt(ignore_index=False, col_level=0)

    df_long_renamed = df_long.reset_index(inplace=False).rename(columns = {'level_0':'gene1','level_1':'Systematic gene name', 'variable':'gene2','value':'pcc'})

    df_long_renamed.loc[df_long_renamed.pcc>=thr,['gene1','gene2','pcc']].to_csv(nw_output, index=False,float_format="%.6f")
    #df = df.rename(columns = {'index':'new column name'})

    df_long_renamed.iloc[:,[0,1]].drop_duplicates().to_csv(strain_ids_output,index=False)

    sgd = pd.read_csv(sgd_file,'\t',header=None)

    #sgd.iloc[:,[0,3]]
    #Save background list for later GO enrichment analysis
    sgd.loc[sgd.iloc[:,3].isin(df_long_renamed['Systematic gene name']),0].to_csv(background_output,index=False, header=None)

def read_y2h_data(input_file, sgd_file, background_output, nw_output):
    network_edgelist = pd.read_csv(input_file,sep='\t',header=None)
    network_edgelist['weight'] = 1
    network_edgelist.rename(columns={0:"gene1",1:"gene2"}).to_csv(nw_output,index=False)

    genelist = [*network_edgelist.iloc[:,0].unique(),*network_edgelist.iloc[:,1].unique()]
    sgd = pd.read_csv(sgd_file,sep='\t',header=None)

    sgd.loc[sgd.iloc[:,3].isin(genelist),0].to_csv(background_output,index=False, header=None)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Reads in a costanzo data file and outputs a network and background file')

    parser.add_argument('--input_type', type=str, default='costanzo', help='Type of input file. Currently only costanzo or y2h is supported')
    parser.add_argument('-i', '--input', help='Input file', required=True)
    parser.add_argument('-s', '--sgd', help='SGD file', required=True)
    parser.add_argument('-b', '--background', help='Background file', required=True)
    parser.add_argument('-n', '--network', help='Network file', required=True)
    parser.add_argument('-t', '--threshold', help='Threshold', required=False, default=0.2, type=float)
    parser.add_argument('-st', '--strain_ids', help='Strain ids file', required=False, default=None)
    
    args = parser.parse_args()
    if args.input_type == 'costanzo':
        read_costanzo_data(args.input, args.sgd, args.background, args.network, args.strain_ids, args.threshold)
    elif args.input_type == 'y2h':
        read_y2h_data(args.input, args.sgd, args.background, args.network)