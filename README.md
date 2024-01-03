# Table of contents
- [Table of contents](#table-of-contents)
- [orthoDB groups for conservation analysis](#orthodb-groups-for-conservation-analysis)
  - [Pipeline overview:](#pipeline-overview)
    - [outside resources used (and links):](#outside-resources-used-and-links)
- [setup TL;DR:](#setup-tldr)
- [Usage](#usage)
- [Advanced stuff](#advanced-stuff)
  - [brief explanation of the orthoDB data](#brief-explanation-of-the-orthodb-data)
  - [advanced configuration](#advanced-configuration)
  - [comments](#comments)
  - [source code](#source-code)
- [setup (more detailed):](#setup-more-detailed)
  - [Download this repository](#download-this-repository)
  - [Download orthoDB database files](#download-orthodb-database-files)
  - [conda environment](#conda-environment)
  - [install local tools in environment](#install-local-tools-in-environment)
  - [generate SQLite databases for orthoDB files](#generate-sqlite-databases-for-orthodb-files)

# orthoDB groups for conservation analysis
This repository contains tools to retrieve and process ortholog groups from a local copy of the orthoDB database files ([link](https://www.orthodb.org/)). The pipeline finds a protein of interest in the database, retrieves its homologous proteins as defined by the orthoDB, and processes the group of homologs in preparation for downstream conservation analysis. <br><br>
- Note: I refer  to these sequences as orthologs at many places throughout this repo but it is probably more accurate to refer to them as homologs. I use the term orthologs because that is what orthoDB calls them. see orthoDB [terminology](https://www.ezlab.org/orthodb_userguide.html#terminology)


## Pipeline overview:
1. **Find a query protein in the orthoDB database** (retrieve the corresponding orthoDB ID)
   - you can find a protein by looking up its uniprot ID in the orthoDB database. You can also set the orthoDB ID manually if you already know it. 
   - *Note: If it can't find an orthoDB ID for a given UniprotID, it doesn't mean that the protein is absent from the orthoDB. It could still be present but was retrieved from a different database and a uniprot ID was not mapped to it. View OrthoDB documentation for more info on where the sequences come from. I have not solved this problem. In the future, it would be nice to develop a way to search for the actual full length sequence using blast or something, if it fails to find the uniprot ID*
2. **Retrieve the orthoDB-defined groups of homologous sequences** (orthogroup IDs) containing the query protein
3. **Select a group based on phylogenetic level**
    - OrthoDB constructs groups of homologs at different phylogenetic levels (e.g. Eukaryota, Metazoa, Vertebrata, etc.). Thus, a protein will probably be part of multiple groups at different levels. Look at the orthoDB page for pcare as an example (https://www.orthodb.org/?query=pcare)
    - For short linear motifs, we have found that the Vertebrata level typically works well for conservation analysis
4. **Filter** out sequences in the group that are too short (relative to the length of the query sequence) or that contain non amino acid characters
5. **Filter to least divergent orthologs (LDOs)**:
   - For each organism in the group, select the sequence that is most similar to the query sequence such that there is only one sequence per organism
6. **Cluster the filtered LDOs using CD-HIT**
7. **Align the clustered sequences** (uses MAFFT by default)
8. **output** the alignment and the ortholog group information in a directory structure that is compatible with the conservation analysis pipeline (link)
    - the group information is output in the form of a json file that can be imported as a python object (see `./examples/` for examples of how to use the object)

### outside resources used (and links):
- [orthoDB database](https://www.orthodb.org/)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [CD-HIT](https://sites.google.com/view/cd-hit)
- [alfpy](https://github.com/aziele/alfpy)
- [BioPython](https://biopython.org/)
- and other python packages ... see `./environment.yml` for a full list of dependencies

# setup TL;DR:
1. download this repository
2. download the orthoDB database files from here: [link](https://data.orthodb.org/download/)
3. navigate to this downloaded repository in terminal (where this README file is located)
4. edit the `.env` file with the location of the orthoDB downloads: `ORTHODB_DATA_DIR=/absolute/path/to/folder/with/orthodb_files/`
5. create a new environment with the dependencies: `conda env create -f environment.yml` <br>
6. activate the environment: `conda activate odb_groups_x86` <br>
7. install the local package: `pip install .` <br>
8. generate the SQLite databases: `bash ./prepare_data.sh` <br>

Check out the longer [setup](#setup-more-detailed) instructions below for more details.

# Usage

**Look at the `./examples/` directory** where there are multiple different applications of the pipeline shown.

There are different ways to use this pipeline depending on your use case:
- You may want to run the pipeline on a small number of genes from a table
- Or you may want to run the pipeline on all of the proteins in an organism's proteome at different phylogenetic levels and use the output as a database for downstream conservation analysis. <br>
- examples of both of these use cases are shown in the `./examples/` directory

The main pipeline is executed via the script: `./src/scripts/odb_group_pipeline.py`
- This script runs the full pipeline for a single gene (see [pipeline overview](#pipeline-overview)) <br>
    - The gene can be specified using a uniprot ID or an odb_gene_id (see [brief explanation of the orthoDB data](#brief-explanation-of-the-orthodb-data) for more info on these ids)
- It can be run as a command line script or imported to be used in another script (like if you are running it on a lot of genes). 

How to use and configure the pipeline is explained in detail in `./examples/readme.md` <br>

# Advanced stuff
## brief explanation of the orthoDB data
You can view the readme file that comes with the orthoDB download for more information. <br>
The data is organized using some internal id numbers. <br>
Here is a brief explanation of the ids and how I've refered to them in the code:
- **odb_gene_id**: An internal orthoDB id. It defines a specific protein sequence in the database. The sequences in the fasta file have the orthoDB id as the sequence id. It is not consistent across versions of orthoDB.
  - example: `9606_0:001c7b`
  - odb_gene_id's are mapped to a variety of other ids corresponding to outside databases (e.g. uniprot, ensembl, etc.) or they were downloaded from some database and have a corresponding id. This information is stored in the `odb11v0_gene_xrefs.tab`/`odb11v0_genes.tab` files.
- **og_id**: An internal orthoDB id. It is probably not consistent across versions of orthoDB. It defines a group of homologous sequences. Each phylogenetic level of homologs has a unique og_id. So a single odb_gene_id probably belongs to multiple og_ids.
  - example og_id: `1567973at7742`
  - example - all of the og_id's that contain the odb_gene_id `9606_0:001c7b`:
    - `605262at9347`, `70995at314146`, `5821at9604`, `1742826at33208`, `5821at314295`, `4349at40674`, `1005199at32523`, `1567973at7742`, `5471876at2759`, `70995at9443`
  - Each of the og_id's is associated with a phylogenetic level at which it was constructed:
      | OG id          | level name       |   total non-redundant count of species underneath |
      |:---------------|:-----------------|--------------------------------------------------:|
      | 5471876at2759  | Eukaryota        |                                              1952 |
      | 1742826at33208 | Metazoa          |                                               817 |
      | 1567973at7742  | Vertebrata       |                                               470 |
      | 1005199at32523 | Tetrapoda        |                                               325 |
      | 4349at40674    | Mammalia         |                                               191 |
      | 605262at9347   | Eutheria         |                                               182 |
      | 70995at314146  | Euarchontoglires |                                                70 |
      | 70995at9443    | Primates         |                                                30 |
      | 5821at314295   | Hominoidea       |                                                 7 |
      | 5821at9604     | Hominidae        |                                                 5 |
      - These correspond with the groups on the website: https://www.orthodb.org/?query=9606_0%3A001c7b

## advanced configuration
you can use your own alignment program by changing the `ALIGNER_EXECUTABLE` variable in the `.env` file

## comments
The main advantages of using these tools:
- the ability to query the orthoDB files quickly (using SQLite databases)
    - The orthoDB files are very large and it is impractical to access them using normal methods (e.g. loading table into python and filtering based on a column value). The key is to pre-index the columns that you are going to filter on, however, it can't be done using normal text files. So I used SQLite to create databases with indexed columns. It makes the lookup time faster by orders of magnitude. The scripts to create the databases here index the columns that are used for lookups in the pipeline.
- go straight from a query protein to a set of orthologs
    - The OrthoDB database can be used for many different analyses. My specific use case is to just retrieve orthologs for a query protein. 
    - Retrieving ortholog groups from the downloaded files requires cross-referencing multiple large files and linking different tables with ids. 
    - I've basically figured all of that out for you so that it's easy to go from a query protein to a set of orthologs.
- reducing ortholog groups to least divergent orthologs (LDOs). 
    - This simplifies a conservation analysis by removing paralogs. It crudely assumes that the most similar sequence to the query sequence is the appropriate sequence in that organism to compare with.
    - You arguably lose information doing this (paralogs), but it greatly simplifies the analysis when you are looking at conservation for a lot of proteins (too many to manually inspect all of the sequences/alignments)
- clustering
    - Clustering reduces sequence redundancy and results in a more even distribution of sequence diversity over the group, which is very helpful in a conservation analysis.
        - Imagine you had a group of 100 homologous sequences: 60 from mammals, all with > 95 % identity and 40 distributed across more distant Vertebrates. The conservation analysis would be dominated by the mammals and would not be very informative. Preclustering would collapse highly similar sequences (>90% identity by default) into one sequence, which would make the analysis more informative.
Notes:
This pipeline could probably benefit from using a workflow manager like snakemake or nextflow. However, I'm not sure if it would be worth the extra effort required for people to understand if they are unfamiliar with those tools. 
I've considered using a config manager like hydra, but it has the same potential downfall. 

## source code
- see `./src/readme.md` for more info on the source code

# setup (more detailed):

## Download this repository
```bash
git clone https://github.com/jacksonh1/orthogroup_generation.git
```

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
- I have included example files in the `./data/orthoDB_sample_data/` directory. These files are from orthoDB version 11.0 and are just subsets of the full files<br>
- Note: *if you use a different version of orthoDB, you will need to change the file names in `./src/local_env_variables/env_variables.py`. File names are hard coded in the `orthoDB_files_object` class*

## conda environment

Create a conda environment with the necessary packages: <br>
`conda env create -f environment.yml` <br>
This will also include CD-HIT and MAFFT from bioconda. <br>


## install local tools in environment

**TL;DR**: 
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
