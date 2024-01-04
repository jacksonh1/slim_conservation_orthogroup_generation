# examples of how to use the pipeline

## examples shown here:

### `./ex1_single_gene/`:

In this example, the pipeline is run for a single gene using the script as a command line script (via a one-line bash script `./ex1_single_gene/run.sh`).

### `./ex2_table_with_uniprot_ids/`:

In this example, a table with a column of uniprot ids is imported and ortholog groups are constructed for each uniprot id. The main pipeline is imported and used in a simple script `./ex2_table_with_uniprot_ids/get_orthologs_from_table.py`. The resulting file paths are mapped back to the original table and the table is written to a new file. <br>

### `./ex3_all_human_genes/`:
In this example, there is a script that runs the pipeline for all of the human genes in the database.
- The pipeline is imported and used in a script `./ex3_all_human_genes/pipeline_all_human_ids.py`
- This script uses multiprocessing
- It also accesses the `src/local_orthoDB_group_tools/sql_queries.py` tools to retrieve odb_gene_ids to run the pipeline on.

### `./pipeline_walkthrough_for_single_gene/`:

This folder contains a notebook that walks through the pipeline step by step to show how the source code works.

###  snakemake alignment? `./`



