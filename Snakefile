
OUTPUT_PATH= "data/interim"
RAW_INPUT_PATH = "data/raw"

rule read_costanzo_data:
    input: f"{RAW_INPUT_PATH}/Data File S3. Genetic interaction profile similarity matrices/cc_ALL.txt"
    params: 
        output_path = OUTPUT_PATH
    output: f"{OUTPUT_PATH}/costanzo_pcc_ALL" 
    script: "scripts/read_costanzo_data.py"

rule create_enm_object:
    input:
        network_file = "data/interim/costanzo_pcc_ALL" ,
        gaf= f"{RAW_INPUT_PATH}/ontology/sgd.gaf",
        obo= f"{RAW_INPUT_PATH}/ontology/go-basic.obo",
        background_file = "data/interim_bak/costanzo_gc_bg.tsv",
        sgd_info = f"{RAW_INPUT_PATH}/ontology/SGD_features.tab"
    params:
        output_path = OUTPUT_PATH
    output: 
        pickle_file= f"{OUTPUT_PATH}/pcc.pickle",
        df_filename= f"{OUTPUT_PATH}/pcc_df.csv"
    script: "scripts/pcc.py"

rule rewire_network:
    input: "{params.output_path}/pcc.pickle"
    params:
        n_sim = 10,
        output_path = OUTPUT_PATH
    output: 
        pickle_file = f"{OUTPUT_PATH}/pcc_{{params.n_sim}}.pickle",
        rewired_dataframe = "rewired_data"
    script: "scripts/rewiring.py"

rule sensor_in_to_out_ratio:
    input: "{params.output_path}/pcc.pickle"
    params:
        output_path = OUTPUT_PATH
    output: f"{OUTPUT_PATH}/sensor_connectivity_df.csv"
    script: "scripts/connectivity.py"