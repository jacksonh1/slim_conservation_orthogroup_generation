# least divergent, clustered ortholog groups from orthoDB
This repository contains tools to retrieve and process ortholog groups from a local copy of the orthoDB database (link) in preparation for downstream conservation analysis. <br>

## Pipeline:
1. Find a query protein in the orthoDB database (retrieve the corresponding orthoDB ID)
2. Retrieve the orthoDB-defined groups of homologous sequences (orthogroup IDs) containing the query protein
3. Select a group based on phylogenetic level
4. Filter out sequences in the group that are too short (relative to the length of the query sequence) or that contain non amino acid characters
5. Filter to least divergent orthologs (LDOs):
   - For each organism in the group, select the sequence that is most similar to the query sequence such that there is only one sequence per organism
6. Cluster the filtered LDOs using CD-HIT
7. Align the clustered sequences using MAFFT
8. output the alignment and the ortholog group information in a directory structure that is compatible with the conservation analysis pipeline (link)


# setup:

## Download orthoDB database files
- download the orthoDB database files from here: [link](https://data.orthodb.org/download/)
  - I used version 11.0 for my project
  - The files you need are:
    - `odb11v0_all_og_fasta.tab.gz`
    - `odb11v0_levels.tab.gz`
    - `odb11v0_species.tab.gz`
    - `odb11v0_level2species.tab.gz`
    - `odb11v0_genes.tab.gz`
    - `odb11v0_gene_xrefs.tab.gz`
    - `odb11v0_OGs.tab.gz`
    - `odb11v0_OG2genes.tab.gz`
- Unzip the files (todo: see if you can work with the compressed files instead)
- create a file called `.env` in this directory (where this README file is located)
  - or open the `.env` file here if it already exists 
- Edit the `.env` file with the following:
  - add a `ORTHODB_DATA_DIR` variable with the absolute path to the directory containing the downloaded files
  - example `.env` file: 
```
ORTHODB_DATA_DIR=/Users/username/project/data/orthodb/odb11v0/
```
- I have included example files in the `./data/orthoDB_sample_data/` directory as an example. These files are from orthoDB version 11.0 and are just subsets of the full files<br>

## conda environment

Create a conda environment with the necessary packages: <br>
`conda env create -f environment.yml` <br>


## install local tools in environment

I wrote this pipeline (code in `./src/`) for it to be installed in the environment as a local package. <br>
I took inspiration from this example: https://github.com/ericmjl/Network-Analysis-Made-Simple <br>
combined with information from here: https://setuptools.pypa.io/en/latest/userguide/package_discovery.html <br>

Install the modules/tools as a local package using setuptools. <br>
first activate the environment: <br>
`conda activate odb_groups_x86` <br>

Make sure that you are in this directory in terminal (where `setup.py` file is located) and run the following command to install: <br>
`pip install .` <br>

If you want to make modifications to the src code, you can install it as an editable install. Then you can make changes to the src files and have them be updated in the environment. <br>
If you want to do this, run the following command instead: <br>
`pip install -e .` <br><br>

## generate SQLite databases for orthoDB files
In this directory, run the following command to generate SQLite databases for the orthoDB files: <br>
```bash
bash ./prepare_data.sh
```

Warning - this will take a while to run and I don't think it can be parallelized (Let me know if this is possible because I would love to know how if so). <br>

