enm
==============================

Elastic network model of the cell and identification of sensor and effector genes
doi: this will be updated

Project Organization
------------

    ├── LICENSE
    ├── Snakefile           <- Snakemake with rules to create intermediate data from raw data 
    ├── README.md          <- The top-level README for people using this project.
    ├── data
    │   ├── interim        <- Intermediate data that has been transformed and used for figure generation.
    │   └── raw            <- The original, immutable data dump.
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── notebooks          <- Jupyter and/or Rmarkdown notebooks for figure generation.
    │
    ├── reports            <- Folder to re-generate notebooks as HTML, PDF, LaTeX, etc.
    ├── reports_done       <- Generated notebooks as HTML, PDF, LaTeX, etc.
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so enm can be imported
    ├── enm                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes enm a Python module
    │   ├── Enm.py         <- Contains Enm class and related functions 
    │   ├── visualize.py    <- Makes enm a Python module
    │   ├── utils.py       <- Contains functions that have been used in this project 
    │
    ├── scripts            <- Contains scripts to generate intermediate data. All scripts could be run with Snakemake 
    └── tox.ini            <- tox file with settings for running tox; see tox.readthedocs.io


--------

This repo contains a Snakefile and thus all pipeline can be run using Snakemake. 

Snakemake will use the raw data provided under `data/raw` to generate all intermediate data and results for figure generation.

The Rmarkdown and Jupyter notebooks under `notebooks` directory can be used to create the figures in the paper. The html files for those notebooks are shared under `reports_done` folder. A rerun of snakemake pipeline will create `html` files to `reports` folder.

To re-create the figures/reports run snakemake:

```bash
snakemake -j10 --use-conda all
```

This will run rules based on the following DAG:
![Snakemake DAG](dag.png)

`--use-conda` directive will download and install necessary packages and run the scripts in a conda environment for both python/jupyter and r/rmarkdown files.

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
