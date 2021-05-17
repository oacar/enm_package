
OUTPUT_PATH= "data/interim"
RAW_INPUT_PATH = "data/raw"
N_SIM=100
PICKLE_FILE_NAME = "data/interim/pcc.pickle"

#rule clean:
#    shell: "rm -rf data/interim/"

rule read_costanzo_data:
    input: f"{RAW_INPUT_PATH}/Data File S3. Genetic interaction profile similarity matrices/cc_ALL.txt"
    params: 
        output_path = OUTPUT_PATH
    output: 
        f"{OUTPUT_PATH}/costanzo_pcc_ALL" ,
        f"{OUTPUT_PATH}/strain_ids.csv" 
    script: "scripts/read_costanzo_data.py"

rule create_enm_object:
    input:
        network_file = f"{OUTPUT_PATH}/costanzo_pcc_ALL" ,
        strain_ids_file = f"{OUTPUT_PATH}/strain_ids.csv"
    params:
        output_path = OUTPUT_PATH
    output: 
        pickle_file= PICKLE_FILE_NAME,
        df_filename= f"{OUTPUT_PATH}/pcc_df.csv"
    script: "scripts/pcc.py"

rule rewire_network:
    input: PICKLE_FILE_NAME
    params:
        n_sim = N_SIM,
    output: 
        pcc_df_random = f"{OUTPUT_PATH}/pcc_df_random_{N_SIM}.csv"
    script: "scripts/rewiring.py"

rule sensor_in_to_out_ratio:
    input: PICKLE_FILE_NAME
    params:
        output_path = OUTPUT_PATH
    output: f"{OUTPUT_PATH}/sensor_connectivity_df.csv"
    script: "scripts/connectivity.py"

rule go_analyze:
    input: f"{OUTPUT_PATH}/figure6_data_051621/rows/prs/{{file}}.csv"
    output: f"{OUTPUT_PATH}/figure6_data_051621/rows/prs/{{file}}.csv.terms"
    run: "analyze.pl ~/enm_package/data/raw/sgd.gaf 5183 ~/enm_package/data/raw/ontology/go-basic.obo {input}" 

rule effector_sensor_go:
    input:
        pickle_file_name= PICKLE_FILE_NAME,
        gaf= f"{RAW_INPUT_PATH}/ontology/sgd.gaf",
        obo= f"{RAW_INPUT_PATH}/ontology/go-basic.obo",
        background_file = "data/interim_bak/costanzo_gc_bg.tsv",
        sgd_info = f"{RAW_INPUT_PATH}/ontology/SGD_features.tab"
    output:
        sensors_df_fname = f"{OUTPUT_PATH}/sensors_df.csv",
        effectors_df_fname = f"{OUTPUT_PATH}/effectors_df.csv",
    script: "scripts/effector_sensor_go.py"