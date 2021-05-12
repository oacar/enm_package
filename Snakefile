
OUTPUT_PATH= "data/interim/",
rule read_costanzo_data:
    input: "data/raw/Data File S3. Genetic interaction profile similarity matrices/cc_ALL.txt"
    output: "data/interim/costanzo_pcc_ALL" 
    script: "scripts/read_costanzo_data.py"

rule create_enm_data:
    input:
        network_file = "data/interim/costanzo_pcc_ALL" ,
        gaf= "data/raw/ontology/sgd.gaf",
        obo= "data/raw/ontology/go-basic.obo",
        background_file = "data/interim_bak/costanzo_gc_bg.tsv",
        sgd_info = "data/raw/ontology/SGD_features.tab"
    params:
        output_path = OUTPUT_PATH
    output: 
        pickle_file= "{params.output_path}/pcc.pickle"
    script: "scripts/pcc.py"

rule rewire_network:
    input: "{params.output_path}/pcc.pickle"
    params:
        n_sim = 10,
        output_path = OUTPUT_PATH
    output: 
        pickle_file = "{params.output_path}/pcc_{params.n_sim}.pickle",
        rewired_dataframe = "rewired_data"
    script: "scripts/rewiring.py"