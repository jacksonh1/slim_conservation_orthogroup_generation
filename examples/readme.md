# examples of how to use the pipeline

The main pipeline (`../src/scripts/odb_group_pipeline.py`), runs the pipeline for a single gene id. <br>
it outputs files with the results in json format and alignments in fasta format (if alignment is specified in the parameters). <br>
The output will look like:
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
When used as a command line script, it takes a few arguments:
```bash
$ python ../src/scripts/odb_group_pipeline.py -h
```
```
usage: odb_group_pipeline.py [-h] (-unid <str> | -odbid <str>) [-c <file>]

run main orthoDB group generation pipeline for a single gene
    processing parameters should be provided in a config file. (-c/--config)
    if no config file is provided, default parameters will be used
    The default parameters are:
- filter_params: {'min_fraction_short_than_query': 0.5}
- og_select_params: {'OG_selection_method': 'level_name', 'OG_level_name': 'Vertebrata'}
- ldo_select_params: {'LDO_selection_method': 'alfpy_google_distance', 'LDO_msa_exe': 'mafft', 'LDO_msa_threads': 8}
- align_params: {'align': False, 'n_align_threads': 8, 'mafft_exe': 'mafft', 'mafft_additional_args': ''}
- cd_hit_exe: cd-hit
- cd_hit_additional_args: 
- main_output_folder: ./odb_group_construction_output
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
Either `-odbid` or `-unid` must be specified. <br>
This is the "query" sequence that the analysis is centered on (provided by either uniprot id or odb gene id). <br>
The `-c` argument is optional. If it is not specified, the default parameters will be used. <br>
Any of the parameters can be specified in the config file and they will overwrite the defaults. Any parameters absent from the config file will take on the default value. More on the parameters below: <br>


## pipeline parameters

The pipeline parameters are specified in a yaml file. <br>
The parameters are:
- `filter_params`: parameters for filtering the orthologs
  - `min_fraction_short_than_query`: the minimum fraction of the query sequence that the ortholog sequence must be in order to be included in the analysis
- `og_select_params:` parameters for selecting the ortholog groups
  - `OG_selection_method`: the method for selecting the ortholog groups. Can be one of:
    - `level_name`: (Default) selects the ortholog groups at a specified taxonomic level
    - `most_species`: select the orthogroup level containing the most species 
- `ldo_select_params`: parameters for selecting the least divergent orthologs (LDOs)
  - `LDO_selection_method`: the method for selecting the LDOs. Can be one of:
    - `alfpy_google_distance`: (Default) uses the alfpy package to calculate the google distance between the query sequence and each ortholog sequence. The ortholog sequence with the smallest google distance is selected as the LDO.
    - `msa`: uses the mafft package to construct 1 multiple sequence alignment all of the ortholog sequences. Then the paralog in each organism with the highest percent identity (PID) to the query sequence is chosen as the LDO.
    - `msa_by_organism`: Uses the mafft package to construct a separate alignment for each organism in the group, composed of the paralogs in each organism and the query sequence. Then the paralog in each organism with the highest PID to the query sequence is chosen as the LDO.
    - `pairwise`: Performs a pairwise alignment between the query sequence and each individual ortholog in the group using BioPython. Selects the paralog in each organism with the highest PID to the query sequence is chosen as the LDO.
    - `LDO_msa_threads`: the number of threads to use for the msa. Default=8. Only used if `LDO_selection_method` is `msa` or `msa_by_organism`.
- `align_params`: parameters for aligning the clustered LDOs
  - `align`: whether or not to align the LDOs. Can be one of:
    - `True`: (Default) align the LDOs
    - `False`: do not align the LDOs
  - `n_align_threads`: the number of threads to use for alignment. Default=8
```yaml 
filter_params:
  min_fraction_short_than_query: 0.5

og_select_params:
  OG_selection_method: level_name
  # OG_level_name: Vertebrata

ldo_select_params:
  LDO_selection_method: alfpy_google_distance
  # LDO_selection_method: msa
  # LDO_selection_method: msa_by_organism
  # LDO_selection_method: pairwise

align_params:
  align: false
  n_align_threads: 8

main_output_folder: odb_group_construction_output
write_files: true
```
```


## examples shown here:

### `./ex1_single_gene/`:

### `./ex2_table_with_uniprot_ids/`:

In this example a table with a column of uniprot ids is imported and ortholog groups are constructed for each uniprot id

### `./ex3_all_human_genes/`:

In this example, there is a script that runs the pipeline for all of the human genes in the database.
The script is imported instead of run as a command line script. I've also written the script to run the pipeline using multiprocessing

### `./pipeline_walkthrough_for_single_gene/`:

This folder contains a notebook that walks through the pipeline step by step to show how the source code works.



###  snakemake alignment? `./`



