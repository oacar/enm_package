
from signal import NSIG


OUTPUT_PATH= "data/interim"
RAW_INPUT_PATH = "data/raw"
N_SIM=100
PICKLE_FILE_NAME = "data/interim/pcc.pickle"

#rule clean:
#    shell: "rm -rf data/interim/"

# rule all:
#     input: 
#     "reports/01-Fig1bcd_3c_4b_5df-052421.html",
#     "reports/02-Figure2-051321.html",
#     "reports/03-Fig3abde_4acd_5b-051821.html"


rule read_costanzo_data:
    input: 
        f"{RAW_INPUT_PATH}/Data File S3. Genetic interaction profile similarity matrices/cc_ALL.txt",
        f"{RAW_INPUT_PATH}/ontology/SGD_features.tab"
    params: 
        output_path = OUTPUT_PATH
    output: 
        f"{OUTPUT_PATH}/costanzo_pcc_ALL" ,
        f"{OUTPUT_PATH}/strain_ids.csv" ,
        f"{OUTPUT_PATH}/go_background_list"
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
        background_file = f"{OUTPUT_PATH}/go_background_list",
        sgd_info = f"{RAW_INPUT_PATH}/ontology/SGD_features.tab"
    output:
        sensors_df_fname = f"{OUTPUT_PATH}/sensors_df.csv",
        effectors_df_fname = f"{OUTPUT_PATH}/effectors_df.csv",
        effector_sensor_combined_go_df = f"{OUTPUT_PATH}/effector_sensor_combined_go_df.csv"
    script: "scripts/effector_sensor_go.py"

rule prs_row_go:
    input: 
        pickle_file_name= PICKLE_FILE_NAME,
        gaf= f"{RAW_INPUT_PATH}/ontology/sgd.gaf",
        obo= f"{RAW_INPUT_PATH}/ontology/go-basic.obo",
        background_file = f"{OUTPUT_PATH}/go_background_list",
        sgd_info = f"{RAW_INPUT_PATH}/ontology/SGD_features.tab"
    output: 
        prs_row=f"{OUTPUT_PATH}/prs_ranked_goa_rows.csv"
    script: "scripts/prs_row_go.py"
rule prs_col_go:
    input: 
        pickle_file_name= PICKLE_FILE_NAME,
        gaf= f"{RAW_INPUT_PATH}/ontology/sgd.gaf",
        obo= f"{RAW_INPUT_PATH}/ontology/go-basic.obo",
        background_file = f"{OUTPUT_PATH}/go_background_list",
        sgd_info = f"{RAW_INPUT_PATH}/ontology/SGD_features.tab"
    output: 
        prs_column=f"{OUTPUT_PATH}/prs_ranked_goa_columns.csv"
    script: "scripts/prs_row_go.py"
rule rwr_row_go:
    input: 
        pickle_file_name= PICKLE_FILE_NAME,
        gaf= f"{RAW_INPUT_PATH}/ontology/sgd.gaf",
        obo= f"{RAW_INPUT_PATH}/ontology/go-basic.obo",
        background_file = f"{OUTPUT_PATH}/go_background_list",
        sgd_info = f"{RAW_INPUT_PATH}/ontology/SGD_features.tab"
    output: 
        rwr_row=f"{OUTPUT_PATH}/rwr_ranked_goa_rows.csv"
    script: "scripts/prs_row_go.py"
rule rwr_col_go:
    input: 
        pickle_file_name= PICKLE_FILE_NAME,
        gaf= f"{RAW_INPUT_PATH}/ontology/sgd.gaf",
        obo= f"{RAW_INPUT_PATH}/ontology/go-basic.obo",
        background_file = f"{OUTPUT_PATH}/go_background_list",
        sgd_info = f"{RAW_INPUT_PATH}/ontology/SGD_features.tab"
    output: 
        rwr_column=f"{OUTPUT_PATH}/rwr_ranked_goa_columns.csv",
    script: "scripts/prs_row_go.py"
rule prs_rwr_compare:
    input: 
        rwr_column=f"{OUTPUT_PATH}/rwr_ranked_goa_columns.csv",
        prs_column=f"{OUTPUT_PATH}/prs_ranked_goa_columns.csv",
        rwr_row=f"{OUTPUT_PATH}/rwr_ranked_goa_rows.csv",
        prs_row=f"{OUTPUT_PATH}/prs_ranked_goa_rows.csv"
    #output: 
    #    plot = f"{OUTPUT_PATH}/figure5B.pdf"
    #script: "scripts/figure5B.py"


rule figure2:
    input:
        pcc_df="data/interim/pcc_df.csv",
        pcc_df_random=f"data/interim/pcc_df_random_{N_SIM}.csv"
    conda:
        "r_env.yml"
    output: "reports/02-Figure2-051321.html"
    script:
        "notebooks/02-Figure2-051321.Rmd"

rule figure3_4_5:
    input:
        pcc_df="data/interim/pcc_df.csv",
        sensor_connectivity_df = "data/interim/sensor_connectivity_df.csv",
        sensors_pcc = "data/interim/sensors_df.csv",
        effector_pcc = "data/interim/effectors_df.csv",
        rwr_column=f"{OUTPUT_PATH}/rwr_ranked_goa_columns.csv",
        prs_column=f"{OUTPUT_PATH}/prs_ranked_goa_columns.csv",
        rwr_row=f"{OUTPUT_PATH}/rwr_ranked_goa_rows.csv",
        prs_row=f"{OUTPUT_PATH}/prs_ranked_goa_rows.csv"
        #pcc_df_random=f"data/interim/pcc_df_random_{N_SIM}.csv"
    conda:
        "r_env.yml"
    output: "reports/03-Fig3abde_4acd_5b-051821.html"
    script:
        "notebooks/03-Fig3abde_4acd_5b-051821.Rmd"

rule figure_networks:
    input: 
        pickle_file_name= PICKLE_FILE_NAME,
        sensors_pcc = "data/interim/sensors_df.csv",
        effector_pcc = "data/interim/effectors_df.csv"
    log:
        # optional path to the processed notebook
        notebook="reports/01-Fig1bcd_3c_4b_5df-052421.ipynb"
    output:
        notebook="reports/01-Fig1bcd_3c_4b_5df-052421.ipynb"
    notebook: "notebooks/01-Fig1bcd_3c_4b_5df-052421.ipynb"

rule figure_networks_html:
    input: 
        "reports/01-Fig1bcd_3c_4b_5df-052421.ipynb"
    output: "reports/01-Fig1bcd_3c_4b_5df-052421.html"
    shell: "jupyter nbconvert {input} --to html"
