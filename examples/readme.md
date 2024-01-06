# examples of how to use the pipeline

## examples shown here:

### `./ex1_single_gene/`:

In this example, the pipeline is run for a single gene using the script as a command line script (via a one-line bash script `./ex1_single_gene/run.sh`).

### `./ex2_table_with_uniprot_ids/`:

In this example, a table with a column of uniprot ids is imported and ortholog groups are constructed for each uniprot id. The main pipeline is imported and used in a simple script `./ex2_table_with_uniprot_ids/get_orthologs_from_table.py`. The resulting file paths are mapped back to the original table and the table is written to a new file. <br>

### `./ex3_all_human_genes/`:
In this example, there is a script that runs the pipeline for all of the human genes in the database. It is more complex than the last one to show how I would implement the pipeline in a more complex scenario. <br>
- pipeline is imported and used in a script `./ex3_all_human_genes/pipeline_all_human_ids.py`
- uses multiprocessing
- accesses the `src/local_orthoDB_group_tools/sql_queries.py` tools to retrieve odb_gene_ids to run the pipeline on.
- dynamically modifies the configuration object in the script to change the ortholog level. In the `multiple_levels` function, the level is changed and the pipeline is run for each level via this line: `config.og_select_params.OG_level_name = og_level`
    - I want to show that the config file is imported as an object and its attributes can be changed programatically.
    - note, that if it is not explicitly defined, the default is used.

Running on the full dataset using 62 cores, this pipeline took a few hours to run and generated around 20 Gb of data


### `./pipeline_walkthrough_for_single_gene/`:

This folder contains a notebook that walks through the pipeline step by step to show how the source code works.

###  snakemake alignment? `./`



