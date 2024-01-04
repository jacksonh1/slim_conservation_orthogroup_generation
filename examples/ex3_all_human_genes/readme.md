# construct groups for all genes in an organism

This is the way that I've used this pipeline mostly. I generated the groups of filtered, clustered, least divergent orthologs (LDOs) for all of the human genes and then accessed the output files for further conservation analysis later.

You could run the script with the command `python ./pipeline_all_human_ids.py`

if you are running on a workstation, I would run it using nohup and save the output to a file like this:
`nohup python ./pipeline_all_human_ids.py > pipeline_all_human_ids.out &`

If using a cluster, I would specifiy the number of cores to use in the script by changing the `N_CORES` variable in the script

If you wanted to do this for another organism, you would just need to change the variable `SPECIES_ID` to the appropriate species code (e.g. `9606_0` for human).

You can find the species code for any organism in orthoDB in the `odb11v0_species.tab` file that you downloaded from the orthoDB website. The file in this repository is actually the full sized file and so you could just use that one as well (`../../data/orthoDB_sample_data/odb11v0_species.tab`).


Or you could do something like this:
```python
import local_env_variables.env_variables as env

for k,v in env.ODB_DATABASE.data_species_dict.items():
    if 'sapiens' in v:
        print(k,v)
```

## more details

If you look at the script in this directory, this line retrieves all of the odb_gene_id's for an organism code:
`odbgeneid_list = sql_queries.get_all_odbids_from_species_id('9606_0')`
The `odbgeneid_list` is then used to run the pipeline for each gene id in the list.
