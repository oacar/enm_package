
from signal import NSIG


OUTPUT_PATH= "data/interim"
NW_TYPE = ['costanzo','y2h']
RAW_INPUT_PATH = "data/raw"
N_SIM=10
#PICKLE_FILE_NAME = f"{OUTPUT_PATH}/pcc.pickle"
SAVE_FIGURES = False
#rule clean:
#    shell: "rm -rf data/interim/"

rule all:
    input:
        expand(f"{OUTPUT_PATH}/{{NW_TYPE}}/{{NW_TYPE}}_effector_sensor_combined_go_df.csv", NW_TYPE=NW_TYPE),
        f"{OUTPUT_PATH}/coessentiality/coessentiality_effector_sensor_combined_go_df_human.csv",
        f"{OUTPUT_PATH}/roguev/roguev_effector_sensor_combined_go_df_pombe.csv",
        expand(f"{OUTPUT_PATH}/{{NW_TYPE}}/{{NW_TYPE}}_df_random_10.csv", NW_TYPE=NW_TYPE),
        f"{OUTPUT_PATH}/coessentiality/coessentiality_df_random_10_human.csv",
        f"{OUTPUT_PATH}/roguev/roguev_df_random_10.csv",
        

# rule all:
#     input: 
#         "reports/01-Fig1bcd_3c_4b_5df-052421.html",
#         "reports/02-Figure2-051321.html",
#         "reports/03-Fig3abde_4acd-051821.html",
#         "reports/04-Signaling-related-effector-sensors.html",
#         #"reports/05-Figs2.ipynb",
#         "reports/08-Fig5abc_figs5.ipynb"
rule create_raw_costanzo:
    input:
        f"{RAW_INPUT_PATH}/Data File S3. Genetic interaction profile similarity matrices/cc_ALL.txt",
    output:
        f"{RAW_INPUT_PATH}/costanzo/costanzo_raw"
    shell:
        "mkdir -v -p data/raw/costanzo && cp '{input}' '{output}'"
rule create_raw_y2h:
    input:
        f"{RAW_INPUT_PATH}/yuri/Y2H_union.txt"
    output:
        f"{RAW_INPUT_PATH}/y2h/y2h_raw"
    shell:"mkdir -p data/raw/y2h && cp '{input}' '{output}'"

rule create_raw_yeast_coex:
    input:
        f"{RAW_INPUT_PATH}/yeast_coex/yeast_AggNet.hdf5"
    output:
        f"{RAW_INPUT_PATH}/yeast_coex/yeast_coex_raw"
    shell:
        "mkdir -p data/raw/yeast_coex && cp '{input}' '{output}'"       

rule clean_pombe_gi_data:
    input: 
        f"{RAW_INPUT_PATH}/{{nw_type}}/raw_data/Dataset_S1.txt",
    params: 
        output_path = OUTPUT_PATH,
        threshold = 0.2,
        
    conda:
        "enm_snakemake.yml"
    output: 
        f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_edgelist.csv" ,
        f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_go_background_list"
    shell:
        "python3 scripts/clean_network_data.py --input_type {wildcards.nw_type} -i {input[0]} -n {output[0]} -b {output[1]} -t {params.threshold}"

rule clean_coessentiality_network_data:
    input: 
        f"{RAW_INPUT_PATH}/{{nw_type}}/{{nw_type}}_raw.txt",
        f"{RAW_INPUT_PATH}/ontology/coessentiality_bg_uniprot_mapping"
    params: 
        output_path = OUTPUT_PATH,
        threshold = 0.2,
        strain_ids_file = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_strain_ids.csv" ,
    conda:
        "enm_snakemake.yml"
    output: 
        f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_edgelist_human.csv" ,
        f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_go_background_list_human"
    shell:
        "python3 scripts/clean_network_data.py --input_type {wildcards.nw_type} -i {input[0]} -s {input[1]} -n {output[0]} -st {params.strain_ids_file} -b {output[1]} -t {params.threshold}"

rule clean_network_data:
    input: 
        f"{RAW_INPUT_PATH}/{{nw_type}}/{{nw_type}}_raw",
        f"{RAW_INPUT_PATH}/ontology/SGD_features.tab"
    params: 
        output_path = OUTPUT_PATH,
        threshold = 0.2,
        strain_ids_file = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_strain_ids.csv" ,
    conda:
        "enm_snakemake.yml"
    output: 
        f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_edgelist.csv" ,
        f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_go_background_list"
    shell:
        "python3 scripts/clean_network_data.py --input_type {wildcards.nw_type} -i {input[0]} -s {input[1]} -n {output[0]} -st {params.strain_ids_file} -b {output[1]} -t {params.threshold}"


# rule read_y2h_data:
#     input: 
#         f"{RAW_INPUT_PATH}/yuri/Y2H_union.txt",
#         f"{RAW_INPUT_PATH}/ontology/SGD_features.tab"
#     params: 
#         output_path = OUTPUT_PATH,
#     conda:
#         "enm_snakemake.yml"
#     output: 
#         f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_edgelist.csv" ,
#         f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_go_background_list"
#     shell:
#         "python3 scripts/clean_network_data.py --input_type {nw_type} -i {input[0]} -s {input[1]} -n {output[0]}  -b {output[2]} "
# rule read_costanzo_data:
#     input: 
#         f"{RAW_INPUT_PATH}/Data File S3. Genetic interaction profile similarity matrices/cc_ALL.txt",
#         f"{RAW_INPUT_PATH}/ontology/SGD_features.tab"
#     params: 
#         output_path = OUTPUT_PATH,
#         threshold = 0.2
#     conda:
#         "enm_snakemake.yml"
#     output: 
#         f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_edgelist.csv" ,
#         f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_strain_ids.csv" ,
#         f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_go_background_list"
#     shell:
#         "python3 scripts/clean_network_data.py --input_type {nw_type} -i {input[0]} -s {input[1]} -n {output[0]} -st {output[1]} -b {output[2]} -t {params.threshold}"
rule create_enm_object_human:
    input:
        network_file = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_edgelist_{{species}}.csv",
    params:
        strain_ids_file = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_strain_ids.csv",
        output_path = f"{OUTPUT_PATH}/{{nw_type}}",
        cluster_matrix = False
    conda:
        "enm_snakemake.yml"
    output: 
        pickle_file= f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_enm_object_{{species}}.pickle",
        df_filename= f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_df_{{species}}.csv"
    shell: "python3 scripts/run_prs.py --network_file {input.network_file} --strain_ids_file {params.strain_ids_file} --output_path {params.output_path} --cluster_matrix {params.cluster_matrix} --output_pickle {output.pickle_file} --output_df {output.df_filename}"
rule create_enm_object:
    input:
        network_file = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_edgelist.csv",
    params:
        strain_ids_file = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_strain_ids.csv",
        output_path = f"{OUTPUT_PATH}/{{nw_type}}",
        cluster_matrix = False
    conda:
        "enm_snakemake.yml"
    output: 
        pickle_file= f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_enm_object.pickle",
        df_filename= f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_df.csv"
    shell: "python3 scripts/run_prs.py --network_file {input.network_file} --strain_ids_file {params.strain_ids_file} --output_path {params.output_path} --cluster_matrix {params.cluster_matrix} --output_pickle {output.pickle_file} --output_df {output.df_filename}"

rule rewire_network_human:
    input: f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_enm_object_human.pickle"
    params:
        n_sim = N_SIM,
    conda:
        "enm_snakemake.yml"
    output: 
        pcc_df_random = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_df_random_{N_SIM}_human.csv"
    shell: "python3 scripts/rewiring.py --pickle_file {input[0]} --random_output_file {output.pcc_df_random} --n_sim {params.n_sim}"
rule rewire_network:
    input: f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_enm_object.pickle"
    params:
        n_sim = N_SIM,
    conda:
        "enm_snakemake.yml"
    output: 
        pcc_df_random = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_df_random_{N_SIM}.csv"
    shell: "python3 scripts/rewiring.py --pickle_file {input[0]} --random_output_file {output.pcc_df_random} --n_sim {params.n_sim}"

# rule sensor_in_to_out_ratio:
#     input: PICKLE_FILE_NAME
#     params:
#         output_path = OUTPUT_PATH
#     conda:
#         "enm_snakemake.yml"
#     output: f"{OUTPUT_PATH}/sensor_connectivity_df.csv"
#     script: "scripts/connectivity.py"

rule effector_sensor_go_human:
    input:
        pickle_file_name= f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_enm_object_{{species}}.pickle",
        gaf= f"{RAW_INPUT_PATH}/ontology/goa_{{species}}.gaf",
        obo= f"{RAW_INPUT_PATH}/ontology/go-basic.obo",
        background_file = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_go_background_list_{{species}}"
    params:
        mapping = f"{RAW_INPUT_PATH}/ontology/{{species}}_name_id_map"
    conda:
        "enm_snakemake.yml"
    output:
        sensors_df_fname = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_sensors_df_{{species}}.csv",
        effectors_df_fname = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_effectors_df_{{species}}.csv",
        effector_sensor_combined_go_df = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_effector_sensor_combined_go_df_{{species}}.csv"
    shell:
        "python3 scripts/effector_sensor_go.py --pickle_file {input.pickle_file_name} --gaf {input.gaf} --obo {input.obo} --background_file {input.background_file} --name_id_map {params.mapping} --sensors_df_fname {output.sensors_df_fname} --effectors_df_fname {output.effectors_df_fname} --effector_sensor_go_df_fname {output.effector_sensor_combined_go_df} --map_column_iloc 1 --id_column_iloc 0"

rule effector_sensor_go_pombe:
    input:
        pickle_file_name= f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_enm_object.pickle",
        gaf= f"{RAW_INPUT_PATH}/ontology/pombase.gaf",
        obo= f"{RAW_INPUT_PATH}/ontology/go-basic.obo",
        background_file = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_go_background_list"
        #sgd_info = None#f"{RAW_INPUT_PATH}/ontology/SGD_features.tab"
    conda:
        "enm_snakemake.yml"
    output:
        sensors_df_fname = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_sensors_df_pombe.csv",
        effectors_df_fname = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_effectors_df_pombe.csv",
        effector_sensor_combined_go_df = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_effector_sensor_combined_go_df_pombe.csv"
    shell:
        "python3 scripts/effector_sensor_go.py --pickle_file {input.pickle_file_name} --gaf {input.gaf} --obo {input.obo} --background_file {input.background_file}  --sensors_df_fname {output.sensors_df_fname} --effectors_df_fname {output.effectors_df_fname} --effector_sensor_go_df_fname {output.effector_sensor_combined_go_df}"

rule effector_sensor_go:
    input:
        pickle_file_name= f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_enm_object.pickle",
        gaf= f"{RAW_INPUT_PATH}/ontology/sgd.gaf",
        obo= f"{RAW_INPUT_PATH}/ontology/go-basic.obo",
        background_file = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_go_background_list",
        sgd_info = f"{RAW_INPUT_PATH}/ontology/SGD_features.tab"
    conda:
        "enm_snakemake.yml"
    output:
        sensors_df_fname = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_sensors_df.csv",
        effectors_df_fname = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_effectors_df.csv",
        effector_sensor_combined_go_df = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_effector_sensor_combined_go_df.csv"
    shell:
        "python3 scripts/effector_sensor_go.py --pickle_file {input.pickle_file_name} --gaf {input.gaf} --obo {input.obo} --background_file {input.background_file} --name_id_map {input.sgd_info}  --sensors_df_fname {output.sensors_df_fname} --effectors_df_fname {output.effectors_df_fname} --effector_sensor_go_df_fname {output.effector_sensor_combined_go_df}"



# rule figure2:
#     input:
#         pcc_df=f"{OUTPUT_PATH}/pcc_df.csv",
#         pcc_df_random = f"{OUTPUT_PATH}/pcc_df_random_{N_SIM}.csv"
#     params:
#         save=SAVE_FIGURES
#     output: "reports/02-Figure2-051321.html"
#     script:
#         "notebooks/02-Figure2-051321.Rmd"

rule figure3_4_human:
    input:
        pcc_df=f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_df_human.csv",
#        sensor_connectivity_df = f"{OUTPUT_PATH}/sensor_connectivity_df.csv",
        sensors_pcc = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_sensors_df_human.csv",
        effector_pcc = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_effectors_df_human.csv",
    params:
        save=False
    output: "reports/03-Fig3abde_4acd-{nw_type}_human.html"
    script:
        "notebooks/03-Fig3abde_4acd-051821.Rmd"
rule figure3_4_pombe:
    input:
        pcc_df=f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_df.csv",
#        sensor_connectivity_df = f"{OUTPUT_PATH}/sensor_connectivity_df.csv",
        sensors_pcc = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_sensors_df_pombe.csv",
        effector_pcc = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_effectors_df_pombe.csv",
    params:
        save=False
    conda:
        "r_env.yml"
    output: "reports/03-Fig3abde_4acd_pombe-{nw_type}.html"
    script:
        "notebooks/03-Fig3abde_4acd-051821.Rmd"
rule figure3_4:
    input:
        pcc_df=f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_df.csv",
#        sensor_connectivity_df = f"{OUTPUT_PATH}/sensor_connectivity_df.csv",
        sensors_pcc = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_sensors_df.csv",
        effector_pcc = f"{OUTPUT_PATH}/{{nw_type}}/{{nw_type}}_effectors_df.csv",
    params:
        save=False
    output: "reports/03-Fig3abde_4acd-{nw_type}.html"
    script:
        "notebooks/03-Fig3abde_4acd-051821.Rmd"


# rule figure_networks:
#     input: 
#         pickle_file_name= PICKLE_FILE_NAME,
#         sensors_pcc = f"{OUTPUT_PATH}/sensors_df.csv",
#         effector_pcc = f"{OUTPUT_PATH}/effectors_df.csv"
#     params:
#         save=SAVE_FIGURES
#     log:
#         # optional path to the processed notebook
#         notebook="reports/01-Fig1bcd_3c_4b_5df-052421.ipynb"
#     conda:
#         "enm_snakemake.yml"
#     output:
#         notebook="reports/01-Fig1bcd_3c_4b_5df-052421.ipynb"
#     notebook: "notebooks/01-Fig1bcd_3c_4b_5df-052421.ipynb"

# rule figure_networks_html:
#     input: 
#         "reports/01-Fig1bcd_3c_4b_5df-052421.ipynb"
#     conda:
#         "enm_snakemake.yml"
#     output: "reports/01-Fig1bcd_3c_4b_5df-052421.html"
#     shell: "jupyter nbconvert {input} --to html"


# rule figs2:
#     input: 
#         strain_ids = 'data/interim/strain_ids.csv',
#         pcc_all = 'data/interim/costanzo_pcc_ALL',
#         pickle_file_name= PICKLE_FILE_NAME,
#         gaf= f"{RAW_INPUT_PATH}/ontology/sgd.gaf",
#         obo= f"{RAW_INPUT_PATH}/ontology/go-basic.obo",
#         background_file = f"{OUTPUT_PATH}/go_background_list",
#         sgd_info = f"{RAW_INPUT_PATH}/ontology/SGD_features.tab"
#     output: 
#         rewired_data_folder = directory('data/interim/rewired_data10test'),
#         notebook="reports/05-Figs2.ipynb"
#     conda:
#         "enm_snakemake.yml"
#     log:
#         notebook="reports/05-Figs2.ipynb"
#     params:
#         save=SAVE_FIGURES,
#         sim_num =10,
#         figure_folder = 'reports/figures/paper_figures_supp/'
#     notebook: "notebooks/05-Figs2.ipynb"


# rule figure_5_s5:
#     input: 
#         gaf= f"{RAW_INPUT_PATH}/ontology/sgd.gaf",
#         obo= f"{RAW_INPUT_PATH}/ontology/go-basic.obo",
#         background_file = f"{OUTPUT_PATH}/go_background_list",
#         sgd_info = f"{RAW_INPUT_PATH}/ontology/SGD_features.tab",
#         pickle_file_name= PICKLE_FILE_NAME,
#         sensors_pcc = f"{OUTPUT_PATH}/sensors_df.csv",
#         effector_pcc = f"{OUTPUT_PATH}/effectors_df.csv"
#     params:
#         save=SAVE_FIGURES
#     log:
#         # optional path to the processed notebook
#         notebook="reports/08-Fig5abc_figs5.ipynb"
#     conda:
#         "enm_snakemake.yml"
#     output:
#         notebook="reports/08-Fig5abc_figs5.ipynb",
#         ec1 = 'data/interim/eff_sens_path1.csv',
#         ec2 = 'data/interim/eff_sens_path2.csv',
#         ec3 = 'data/interim/eff_sens_path3.csv',
#         combined_data_for_colors = 'data/interim/eff_sens_combined_for_coloring.csv'

#     notebook: "notebooks/08-Fig5abc_figs5.ipynb"

# rule eff_sens_signaling:
#     input: 
#         gaf= f"{RAW_INPUT_PATH}/ontology/sgd.gaf",
#         obo= f"{RAW_INPUT_PATH}/ontology/go-basic.obo",
#         background_file = f"{OUTPUT_PATH}/go_background_list",
#         sgd_info = f"{RAW_INPUT_PATH}/ontology/SGD_features.tab",
#         sensors_pcc = f"{OUTPUT_PATH}/sensors_df.csv",
#         effector_pcc = f"{OUTPUT_PATH}/effectors_df.csv"
#     log:
#         # optional path to the processed notebook
#         notebook="reports/04-Signaling-related-effector-sensors.ipynb"
#     conda:
#         "enm_snakemake.yml"
#     output:
#         notebook="reports/04-Signaling-related-effector-sensors.ipynb",
#         sensors_signaling_df = 'data/interim/signaling_related_sensors.csv'
#     notebook: "notebooks/04-Signaling-related-effector-sensors.ipynb"


# rule signaling_related_sensors_html:
#     input: 
#         "reports/04-Signaling-related-effector-sensors.ipynb"
#     conda:
#         "enm_snakemake.yml"
#     output: "reports/04-Signaling-related-effector-sensors.html"
#     shell: "jupyter nbconvert {input} --to html"