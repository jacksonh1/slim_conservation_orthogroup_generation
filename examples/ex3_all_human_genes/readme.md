# construct groups for all genes in an organism

This is the way that I've used the pipeline mostly. I generated groups of filtered, clustered, least divergent orthologs (LDOs) for all of the human genes and then accessed the output files for further conservation analysis later. <br>
In this way I used this output as a database of alignments to annotate tables with conservation scores.

You could run the script with the command `python ./pipeline_all_human_ids.py` or use the command line interface version of the script `python ./pipeline_all_human_ids_CLI_version.py`. I'll assume that you're using the command line version.

to see the command line parameters:
`$python pipeline_all_human_ids_CLI_version.py --help`
```
run the pipeline for all genes in an organism

options:
  -h, --help            show this help message and exit
  -c <file>, --config <file>
                        path to config file (default: None)
  -n <int>, --n_cores <int>
                        number of cores to use (default: 6)
  -s <str>, --species_id <str>
                        species id to use (default: 9606_0)
  -o, --overwrite       if flag is provided and the main_output_folder exists, it will be removed and overwritten by the new files. Otherwise, an
                        error will be raised if the folder exists (default: False)
```

if you are running on a workstation, I would run it using nohup and save the output to a file like this:
`nohup python ./pipeline_all_human_ids_CLI_version.py -c "./params.yml" --species_id "9606_0" > pipeline_all_human_ids.out &`

If using a cluster, I would specifiy the number of cores to use in the script with the `--n_cores` parameter

If you wanted to do this for another organism, you would just need to change the variable `--species_id` to the appropriate species code (e.g. `9606_0` for human).

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
`odbgeneid_list = sql_queries.get_all_odbids_from_species_id(species_id)`
The `odbgeneid_list` is then used to run the pipeline for each gene id in the list.<br>
The function `create_filemap.create_filemap(` creates a json file that maps the id to the generated files.


