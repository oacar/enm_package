
from signal import NSIG


OUTPUT_PATH= "data/interim"
RAW_INPUT_PATH = "data/raw"
N_SIM=100
PICKLE_FILE_NAME = f"{OUTPUT_PATH}/pcc.pickle"
SAVE_FIGURES = True
#rule clean:
#    shell: "rm -rf data/interim/"

rule all:
    input: 
        "reports/01-Fig1bcd_3c_4b_5df-052421.html",
        "reports/02-Figure2-051321.html",
        "reports/03-Fig3abde_4acd-051821.html",
        "reports/04-Signaling-related-effector-sensors.ipynb",
        #"reports/05-Figs2.ipynb",
        "reports/08-Fig5abc_figs5.ipynb"


rule read_costanzo_data:
    input: 
        f"{RAW_INPUT_PATH}/Data File S3. Genetic interaction profile similarity matrices/cc_ALL.txt",
        f"{RAW_INPUT_PATH}/ontology/SGD_features.tab"
    params: 
        output_path = OUTPUT_PATH

    conda:
        "enm_snakemake.yml"
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
        output_path = OUTPUT_PATH,
        cluster_matrix = True
    conda:
        "enm_snakemake.yml"
    output: 
        pickle_file= PICKLE_FILE_NAME,
        df_filename= f"{OUTPUT_PATH}/pcc_df.csv"
    script: "scripts/pcc.py"

rule rewire_network:
    input: PICKLE_FILE_NAME
    params:
        n_sim = N_SIM,
    conda:
        "enm_snakemake.yml"
    output: 
        pcc_df_random = f"{OUTPUT_PATH}/pcc_df_random_{N_SIM}.csv"
    script: "scripts/rewiring.py"

rule sensor_in_to_out_ratio:
    input: PICKLE_FILE_NAME
    params:
        output_path = OUTPUT_PATH
    conda:
        "enm_snakemake.yml"
    output: f"{OUTPUT_PATH}/sensor_connectivity_df.csv"
    script: "scripts/connectivity.py"

rule effector_sensor_go:
    input:
        pickle_file_name= PICKLE_FILE_NAME,
        gaf= f"{RAW_INPUT_PATH}/ontology/sgd.gaf",
        obo= f"{RAW_INPUT_PATH}/ontology/go-basic.obo",
        background_file = f"{OUTPUT_PATH}/go_background_list",
        sgd_info = f"{RAW_INPUT_PATH}/ontology/SGD_features.tab"
    conda:
        "enm_snakemake.yml"
    output:
        sensors_df_fname = f"{OUTPUT_PATH}/sensors_df.csv",
        effectors_df_fname = f"{OUTPUT_PATH}/effectors_df.csv",
        effector_sensor_combined_go_df = f"{OUTPUT_PATH}/effector_sensor_combined_go_df.csv"
    script: "scripts/effector_sensor_go.py"



rule figure2:
    input:
        pcc_df=f"{OUTPUT_PATH}/pcc_df.csv",
        pcc_df_random = f"{OUTPUT_PATH}/pcc_df_random_{N_SIM}.csv"
    params:
        save=SAVE_FIGURES
    output: "reports/02-Figure2-051321.html"
    script:
        "notebooks/02-Figure2-051321.Rmd"

rule figure3_4:
    input:
        pcc_df=f"{OUTPUT_PATH}/pcc_df.csv",
        sensor_connectivity_df = f"{OUTPUT_PATH}/sensor_connectivity_df.csv",
        sensors_pcc = f"{OUTPUT_PATH}/sensors_df.csv",
        effector_pcc = f"{OUTPUT_PATH}/effectors_df.csv",
    params:
        save=SAVE_FIGURES
    output: "reports/03-Fig3abde_4acd-051821.html"
    script:
        "notebooks/03-Fig3abde_4acd-051821.Rmd"


rule figure_networks:
    input: 
        pickle_file_name= PICKLE_FILE_NAME,
        sensors_pcc = f"{OUTPUT_PATH}/sensors_df.csv",
        effector_pcc = f"{OUTPUT_PATH}/effectors_df.csv"
    params:
        save=SAVE_FIGURES
    log:
        # optional path to the processed notebook
        notebook="reports/01-Fig1bcd_3c_4b_5df-052421.ipynb"
    conda:
        "enm_snakemake.yml"
    output:
        notebook="reports/01-Fig1bcd_3c_4b_5df-052421.ipynb"
    notebook: "notebooks/01-Fig1bcd_3c_4b_5df-052421.ipynb"

rule figure_networks_html:
    input: 
        "reports/01-Fig1bcd_3c_4b_5df-052421.ipynb"
    conda:
        "enm_snakemake.yml"
    output: "reports/01-Fig1bcd_3c_4b_5df-052421.html"
    shell: "jupyter nbconvert {input} --to html"


rule figs2:
    input: 
        strain_ids = 'data/interim/strain_ids.csv',
        pcc_all = 'data/interim/costanzo_pcc_ALL',
        pickle_file_name= PICKLE_FILE_NAME,
        gaf= f"{RAW_INPUT_PATH}/ontology/sgd.gaf",
        obo= f"{RAW_INPUT_PATH}/ontology/go-basic.obo",
        background_file = f"{OUTPUT_PATH}/go_background_list",
        sgd_info = f"{RAW_INPUT_PATH}/ontology/SGD_features.tab"
    output: 
        rewired_data_folder = 'data/interim/rewired_data',
        notebook="reports/05-Figs2.ipynb"
    logs:
        notebook="reports/05-Figs2.ipynb"
    conda:
        "enm_snakemake.yml"
    params:
        save=SAVE_FIGURES,
        sim_num =10,
        figure_folder = 'reports/figures/paper_figures_supp/'
    notebook: "notebooks/05-Figs2.ipynb"


rule figure_5_s5:
    input: 
        gaf= f"{RAW_INPUT_PATH}/ontology/sgd.gaf",
        obo= f"{RAW_INPUT_PATH}/ontology/go-basic.obo",
        background_file = f"{OUTPUT_PATH}/go_background_list",
        sgd_info = f"{RAW_INPUT_PATH}/ontology/SGD_features.tab",
        pickle_file_name= PICKLE_FILE_NAME,
        sensors_pcc = f"{OUTPUT_PATH}/sensors_df.csv",
        effector_pcc = f"{OUTPUT_PATH}/effectors_df.csv"
    params:
        save=SAVE_FIGURES
    log:
        # optional path to the processed notebook
        notebook="reports/08-Fig5abc_figs5.ipynb"
    conda:
        "enm_snakemake.yml"
    output:
        notebook="reports/08-Fig5abc_figs5.ipynb",
        ec1 = 'data/interim/eff_sens_path1.csv',
        ec2 = 'data/interim/eff_sens_path2.csv',
        ec3 = 'data/interim/eff_sens_path3.csv',
        combined_data_for_colors = 'data/interim/eff_sens_combined_for_coloring.csv'

    notebook: "notebooks/08-Fig5abc_figs5.ipynb"

rule eff_sens_signaling:
    input: 
        gaf= f"{RAW_INPUT_PATH}/ontology/sgd.gaf",
        obo= f"{RAW_INPUT_PATH}/ontology/go-basic.obo",
        background_file = f"{OUTPUT_PATH}/go_background_list",
        sgd_info = f"{RAW_INPUT_PATH}/ontology/SGD_features.tab",
        sensors_pcc = f"{OUTPUT_PATH}/sensors_df.csv",
        effector_pcc = f"{OUTPUT_PATH}/effectors_df.csv"
    log:
        # optional path to the processed notebook
        notebook="reports/04-Signaling-related-effector-sensors.ipynb"
    conda:
        "enm_snakemake.yml"
    output:
        notebook="reports/04-Signaling-related-effector-sensors.ipynb",
        sensors_signaling_df = 'data/interim/signaling_related_sensors.csv'
    notebook: "notebooks/04-Signaling-related-effector-sensors.ipynb"