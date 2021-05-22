enm
==============================

ENM applications on genetic networks

Project Organization
------------

    ├── LICENSE
    ├── Snakefile           <- Snakemake with rules to create intermediate data from raw data 
    ├── README.md          <- The top-level README for people using this project.
    ├── data
    │   ├── interim        <- Intermediate data that has been transformed.
    │   └── raw            <- The original, immutable data dump.
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── notebooks          <- Jupyter and/or Rmarkdown notebooks for figure generation.
    │
    ├── reports            <- Generated notebooks as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
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

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
