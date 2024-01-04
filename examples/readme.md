# examples of how to use the pipeline

## overview
The file `../src/local_scripts/odb_group_pipeline.py`, runs the pipeline for a single gene. This should be the only file that you need to access run the pipeline<br>
it outputs files with the results in json format and alignments in fasta format (if alignment is specified in the parameters). <br>
The output will look like this:
```bash
main_output_folder/
    ├── info_jsons/
    │   ├── 9606_0_001c7b_Eukaryota_5471876at2759_info.json
    │   ├── ...
    ├── ├── {odb_gene_id}_{og_level}_{og_id}_info.json
    │---alignments/
    │   ├── 9606_0_001c7b_Eukaryota_5471876at2759_clustered_ldos_aln.fasta
    │   ├── ...
    ├── ├── {odb_gene_id}_{og_level}_{og_id}_clustered_ldos_aln.fasta

```
The json files contain the information about the ortholog groups and are used for further conservation analysis. <br><br>

### using the tools as a command line script

The easiest way to explain how it works and how to use it is when it is used as a command line script. Here is the help message for the script:
```bash
$ python ../src/local_scripts/odb_group_pipeline.py --help
```
```
usage: odb_group_pipeline.py [-h] (-unid <str> | -odbid <str>) [-c <file>]

run main orthoDB group generation pipeline for a single gene
    processing parameters should be provided in a config file. (-c/--config)
    if no config file is provided, default parameters will be used
    The default parameters are:
- filter_params: {'min_fraction_shorter_than_query': 0.5}
- og_select_params: {'OG_selection_method': 'level_name', 'OG_level_name': 'Vertebrata'}
- ldo_select_params: {'LDO_selection_method': 'alfpy_google_distance', 'LDO_mafft_threads': 8}
- align_params: {'align': False, 'n_align_threads': 8}
- main_output_folder: ./processed_odb_groups_output
- write_files: True

options:
  -h, --help            show this help message and exit
  -unid <str>, --uniprot_id <str>
                        the uniprot id of the gene of interest
  -odbid <str>, --odb_gene_id <str>
                        the odb gene id of the gene of interest (e.g. "9606_0:001c7b")
  -c <file>, --config <file>
                        path to config file, default=None
```
- Either `-odbid` or `-unid` must be specified. <br>
  - This is the "query" sequence that the analysis is centered on (provided by either uniprot id or odb gene id). It is the sequence that you are retrieving orthologs for.<br>
- The `-c` argument is optional and specifies the location of a configuration file which can be used to specify the pipeline parameters. If it is not specified, the default parameters will be used. <br>
  - Any of the parameters can be specified in the config file and they will overwrite the defaults. Any parameters absent from the config file will take on the default value. More on the parameters below: <br>

You would probably very rarely want to run the pipeline as a command line script, but it is useful for understanding how the pipeline works. <br>
It would also be useful if you wanted to run the pipeline using a workflow manager like snakemake or nextflow. <br>

### using the tools as a module

importing the script as a module is very similar to using it as a command line script, see `./ex2_table_with_uniprot_ids/` and `./ex3_all_human_genes/` for examples.

In general you import the script as a module:
```python
import local_scripts.odb_group_pipeline as pipeline
```
This import should work from anywhere because you installed the src code as a package. <br>
The main function is `pipeline.main_pipeline` but it requires a configuration file to be loaded first. <br>
This is done with the `pipeline.load_config` function:
```python
config = pipeline.load_config("./params.yml")
```
Then you can run the pipeline with the `pipeline.main_pipeline` function, which takes the `config` object as an argument and either a uniprot id or odb gene id as an argument:
```python
pipeline.main_pipeline(config, odb_gene_id="9606_0:002f40")
```
or
```python
pipeline.main_pipeline(config, uniprot_id="Q8TC90")
```
In this way, you can run the pipeline for any number of genes in a script as is shown is ex2 and ex3 <br>

## pipeline parameters

The pipeline parameters are specified in a yaml file. <br>
If you are not familiar with yaml, it is fairly easy to understand by just looking at an example. Here is an example yaml file:<br>
```yaml 
filter_params:
  min_fraction_shorter_than_query: 0.5

og_select_params:
  OG_selection_method: level_name
  OG_level_name: Vertebrata

ldo_select_params:
  LDO_selection_method: alfpy_google_distance

align_params:
  align: false
  n_align_threads: 8

main_output_folder: processed_odb_groups_output
write_files: true
```
More configuration files are used in the examples here and you can also just copy and edit those for your own use (*e.g.* `./ex1_single_gene/params.yml`).

### pipeline parameters explained
Here is an explanation of the parameters and what they do:
- `filter_params`:
  - `min_fraction_shorter_than_query`: A number between 0 and 1. Sequences that are shorter than `min_fraction_shorter_than_query`*(`length of query sequence`) will be filtered out. <br>Default=0.5
- `og_select_params`:
  - `OG_selection_method`: the method for selecting the ortholog groups. Can be one of:
    - `level_name`: (Default) selects the ortholog groups at a specified taxonomic level
    - `most_species`: select the orthogroup level containing the most species 
- `ldo_select_params`:
  - `LDO_selection_method`: the method for selecting the LDOs. Can be one of:
    - `alfpy_google_distance`: (Default) uses the alfpy package to calculate the google distance between the query sequence and each ortholog sequence. The ortholog sequence with the smallest google distance is selected as the LDO.
    - `msa`: uses the mafft package to construct 1 multiple sequence alignment all of the ortholog sequences. Then the paralog in each organism with the highest percent identity (PID) to the query sequence is chosen as the LDO.
    - `msa_by_organism`: Uses the mafft package to construct a separate alignment for each organism in the group, composed of the paralogs in each organism and the query sequence. Then the paralog in each organism with the highest PID to the query sequence is chosen as the LDO.
    - `pairwise`: Performs a pairwise alignment between the query sequence and each individual ortholog in the group using BioPython. Selects the paralog in each organism with the highest PID to the query sequence is chosen as the LDO.
    - `LDO_msa_threads`: the number of threads to use for the msa. Default=8. Only used if `LDO_selection_method` is `msa` or `msa_by_organism`.
- `align_params`:
  - `align`: whether or not to align the clustered LDOs. Can be one of:
    - true: (Default) align the LDOs
    - false: do not align the LDOs
  - `n_align_threads`: the number of threads to use for alignment. Default=8
- `main_output_folder`: the folder to write the output files to. Default=`./processed_odb_groups_output`
- `write_files`: whether or not to write the output files. Can be one of:
  - true: (Default) write the output files
  - false: do not write the output files


## examples shown here:

### `./ex1_single_gene/`:

In this example, the pipeline is run for a single gene using the script as a command line script (via a one-line bash script `./ex1_single_gene/run.sh`).

### `./ex2_table_with_uniprot_ids/`:

In this example a table with a column of uniprot ids is imported and ortholog groups are constructed for each uniprot id. The resulting file paths are mapped back to the original table and the table is written to a new file. <br>

### `./ex3_all_human_genes/`:
In this example, there is a script that runs the pipeline for all of the human genes in the database.
The script is imported instead of run as a command line script. I've also written the script to run the pipeline using multiprocessing

### `./pipeline_walkthrough_for_single_gene/`:

This folder contains a notebook that walks through the pipeline step by step to show how the source code works.

###  snakemake alignment? `./`



