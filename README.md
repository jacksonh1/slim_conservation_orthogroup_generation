# Table of contents
- [Table of contents](#table-of-contents)
- [orthoDB groups for conservation analysis](#orthodb-groups-for-conservation-analysis)
  - [Pipeline overview:](#pipeline-overview)
  - [comments](#comments)
- [setup TLDR:](#setup-tldr)
- [setup:](#setup)
  - [Download orthoDB database files](#download-orthodb-database-files)
  - [conda environment](#conda-environment)
  - [install local tools in environment](#install-local-tools-in-environment)
  - [generate SQLite databases for orthoDB files](#generate-sqlite-databases-for-orthodb-files)
- [Usage](#usage)
- [Parameters](#parameters)
  - [alignment](#alignment)

# orthoDB groups for conservation analysis
This repository contains tools to retrieve and process ortholog groups from a local copy of the orthoDB database ([link](https://www.orthodb.org/)) in preparation for downstream conservation analysis. <br><br>

## Pipeline overview:
1. **Find a query protein in the orthoDB database** (retrieve the corresponding orthoDB ID)
   - you can find a protein by looking up its uniprot ID in the orthoDB database. You can also set the orthoDB ID manually if you already know it. 
   - *Note: If it can't find an orthoDB ID for a given UniprotID, it doesn't mean that the protein is absent from the orthoDB. It could still be present but was retrieved from a different database and a uniprot ID was not mapped to it. View OrthoDB documentation for more info on where the sequences come from. I have not solved this problem. In the future, it would be nice to develop a way to search for the actual full length sequence using blast or something, if it fails to find the uniprot ID*
2. **Retrieve the orthoDB-defined groups of homologous sequences** (orthogroup IDs) containing the query protein
3. **Select a group based on phylogenetic level**
4. **Filter** out sequences in the group that are too short (relative to the length of the query sequence) or that contain non amino acid characters
5. **Filter to least divergent orthologs (LDOs)**:
   - For each organism in the group, select the sequence that is most similar to the query sequence such that there is only one sequence per organism
6. **Cluster the filtered LDOs using CD-HIT**
7. **Align the clustered sequences** (uses MAFFT by default)
8. **output** the alignment and the ortholog group information in a directory structure that is compatible with the conservation analysis pipeline (link)
    - the group information is output in the form of a json file that can be imported as a python object (see `./examples/` for examples of how to use the object)

## comments
The main advantages of using these tools:
- the ability to query the orthoDB files quickly (using SQLite databases)
    - The orthoDB files are very large and it is impractical to access them using normal methods (e.g. loading table into python and filtering based on a column value). The key is to pre-index the columns that you are going to filter on, however, it can't be done using normal text files. So I used SQLite to create databases with indexed columns. It makes the lookup time faster by orders of magnitude. The scripts to create the databases here index the columns that are used for lookups in the pipeline.
- go straight from a query protein to a set of orthologs
    - The OrthoDB database can be used for many different analyses. My specific use case is to just retrieve orthologs for a query protein. 
    - Retrieving ortholog groups from the downloaded files requires cross-referencing multiple large files and linking different tables with ids. 
    - I've basically figured all of that out for you so that it's easy to go from a query protein to a set of orthologs.
- reducing ortholog groups to least divergent orthologs (LDOs). 
    - This simplifies a conservation analysis by removing paralogs. It crudely assumes that the most similar sequence to the query sequence is sequence in that organism to compare with.
    - You arguably lose information doing this, but it greatly simplifies the analysis when you are looking at conservation for a lot of proteins (too many to manually inspect all of the sequences/alignments)
- clustering
    - Clustering reduces sequence redundancy and results in a more even distribution of sequence diversity over the group, which is very helpful in a conservation analysis.
        - Imagine you had a group of 100 homologs: 60 from primates with 99 % identity, 20 from other mammals, and 20 distributed across more distant Eukaryotes. The conservation analysis would be dominated by the primates and would not be very informative. The clustering would collapse the 60 primate sequences into one sequence, which would make the analysis more informative.


# setup TLDR:
1. download this repository
2. download the orthoDB database files from here: [link](https://data.orthodb.org/download/)
3. navigate to this downloaded repository in terminal (where this README file is located)
4. edit the `.env` file with the location of the orthoDB downloads: `ORTHODB_DATA_DIR=/absolute/path/to/folder/with/orthodb_files/`
5. create a new environment with the dependencies: `conda env create -f environment.yml` <br>
6. activate the environment: `conda activate odb_groups_x86` <br>
7. install the local package: `pip install .` <br>
8. generate the SQLite databases: `bash ./prepare_data.sh` <br>

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
ORTHODB_DATA_DIR=/Users/username/project/data/orthodb/odb11v0aa/
```
- I have included example files in the `./data/orthoDB_sample_data/` directory. These files are from orthoDB version 11.0 and are just subsets of the full files<br>
- Note: *if you use a different version of orthoDB, you will need to change the file names in `./src/local_env_variables/env_variables.py`. File names are hard coded in the `orthoDB_files_object` class*

## conda environment

Create a conda environment with the necessary packages: <br>
`conda env create -f environment.yml` <br>
This will also include CD-HIT and MAFFT from bioconda. <br>


## install local tools in environment

**TLDR**: 
- activate the environment: `conda activate odb_groups_x86`
- run `pip install .` in this directory (where `setup.py` file is located) <br>

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

*Note: I am creating a separate database for each file. You could easily make one database with all of the tables. I tried this but it was significantly slower to query and I don't know why.* <br>

# Usage

There are two main ways to use this pipeline: <br>
1. generate ortholog groups for a small number of query proteins from a table
2. generate ortholog groups for every protein in a organism's proteome

example for choosing levels: https://www.orthodb.org/?query=pcare


# Parameters

## alignment
you can use your own alignment program by changing the `ALIGNER_EXECUTABLE` variable in the `.env` file
