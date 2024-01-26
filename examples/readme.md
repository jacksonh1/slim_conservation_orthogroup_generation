# examples of how to use the pipeline

## examples shown here:

### `./uniprotid_mapping/`
A simple script to map uniprot ids to odb_gene_ids. Takes a table with a column of uniprot ids and tries to map each uniprot id to a odb_gene_id. outputs a copy of the table with an added column of odb_gene_ids<br>
run `./uniprotid_mapping/map_uniprotid.py --help` to see more.

### `./ex1_single_gene/`:

In this example, the pipeline is run for a single gene using the script as a command line script (via a one-line bash script `./ex1_single_gene/run.sh`).

### `./ex2_table_with_uniprot_ids/`:

In this example, a table with a column of uniprot ids is imported and ortholog groups are constructed for each uniprot id. The main pipeline is imported and used in a simple script `./ex2_table_with_uniprot_ids/get_orthologs_from_table.py`. The resulting file paths are mapped back to the original table and the table is written to a new file. <br>

### `./ex3_all_human_genes/`:
In this example, there is a script that runs the pipeline for all of the human genes in the database. It is more complex than the last one to show how I would implement the pipeline in a more complex scenario. <br>
- pipeline is imported and used in a script `./ex3_all_human_genes/pipeline_all_ids_in_species.py`
- uses multiprocessing
- accesses the `src/local_orthoDB_group_pipeline/sql_queries.py` tools to retrieve all of the human odb_gene_ids to run the pipeline on.
- dynamically modifies the pipeline parameters object in the script to change the ortholog level. In the `multiple_levels` function, the level is changed and the pipeline is run for each level via this line: `config.og_select_params.OG_level_name = og_level`
    - I want to show that the config file is imported as an object and its attributes can be changed programatically.
    - note, that if it is not explicitly defined, the default is used.
- I also included a check to make sure that the specified output directory doesn't already exist. If it does, the script will exit. This is to prevent overwriting data. <br>
- A filemap is also created after the pipeline is run using the tools in `../src/local_scripts/create_filemap.py`

- `./ex3_all_human_genes/pipeline_all_ids_in_species_CLI_version.py` is the same script as `./ex3_all_human_genes/pipeline_all_ids_in_species.py` but it is run as a command line script. I included this as a separate script because I am afraid that the argparse stuff could confuse new users<br>

Running on the full dataset using 62 cores, this pipeline took a few hours to run and generated around 10 Gb of data. When I ran it with align=True, it generated around 100 Gbs of data. <br>

### `./pipeline_walkthrough_for_single_gene/`:

This folder contains a notebook that walks through the pipeline step by step to show how the source code works.

# workflow example

1. database construction: construct groups for all human genes (using the script in `./ex3_all_human_genes/`)
    - generates a "database" of:
        - orthogroup info files `info_jsons/`
        - alignments: `alignments/`
        - file map: `filemap.json`
2. for a table of candidate motifs from human genes with a column of uniprot ids, map the uniprot ids to odb_gene_ids (using the script in `./uniprotid_mapping/`).
3. use the database filemap (e.g. `filemap.json`) to access the database files from the odb_gene_ids

