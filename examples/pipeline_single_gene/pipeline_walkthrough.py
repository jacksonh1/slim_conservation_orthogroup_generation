# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: odb_groups_x86
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Pipeline overview:
# 1. **Find a query protein in the orthoDB database** (retrieve the corresponding orthoDB ID)
#    - you can find a protein by looking up its uniprot ID in the orthoDB database. You can also set the orthoDB ID manually if you already know it. 
#    - *Note: If it can't find an orthoDB ID for a given UniprotID, it doesn't mean that the protein is absent from the orthoDB. It could still be present but was retrieved from a different database and a uniprot ID was not mapped to it. View OrthoDB documentation for more info on where the sequences come from. I have not solved this problem. In the future, it would be nice to develop a way to search for the actual full length sequence using blast or something, if it fails to find the uniprot ID*
# 2. **Retrieve the orthoDB-defined groups of homologous sequences** (orthogroup IDs) containing the query protein
# 3. **Select a group based on phylogenetic level**
# 4. **Filter** out sequences in the group that are too short (relative to the length of the query sequence) or that contain non amino acid characters
# 5. **Filter to least divergent orthologs (LDOs)**:
#    - For each organism in the group, select the sequence that is most similar to the query sequence such that there is only one sequence per organism
# 6. **Cluster the filtered LDOs using CD-HIT**
# 7. **Align the clustered sequences** (uses MAFFT by default)
# 8. **output** the alignment and the ortholog group information in a directory structure that is compatible with the conservation analysis pipeline (link)
#     - the group information is output in the form of a json file that can be imported as a python object (see `./examples/` for examples of how to use the object)
#

# %% [markdown]
# ---

# %% [markdown]
# Here we will go through each step and explain what the pipeline does. <br>
# You don't have to worry about using any of this code directly, since you should be able to just run the script but it might be useful to understand what is going on under the hood if you want to write your own script.<br>
# There main script in `../orthoDB_construct_groups_for_all_human_genes/` combines all of the steps in a single script, using a configuration file to specify all of the parameters.<br>
# This notebook will explain what each of the parameters does?

# %% [markdown]
# # imports

# %%
import pandas as pd
from Bio import AlignIO, SeqIO
import local_env_variables.env_variables as env
from local_orthoDB_group_tools import (sql_queries, filters, uniprotid_search, og_selection, cluster, find_LDOs)

# %load_ext autoreload
# %autoreload 2

# %% [markdown]
# # Run parameters

# %%
min_fraction_of_query_length = 0.5
level='Vertebrata'
LDO_selection_method = 'alfpy_google_distance'

# %% [markdown]
# # 1. find a query protein in the orthoDB database

# %% [markdown]
# ## search for the uniprot id: Q8TC90

# %%
uniprot_id = 'Q8TC90'
print(uniprotid_search.uniprotid_2_odb_gene_id(uniprot_id))
odb_gene_id = uniprotid_search.uniprotid_2_odb_gene_id(uniprot_id)

# %% [markdown]
# # 2. get the ortholog groups (og_ids) that contain the query protein

# %%
og_ids = sql_queries.odb_gene_id_2_ogid_list(odb_gene_id)
print(og_ids)

# %%
# og_id information
og_selection.get_available_ogs(odb_gene_id)

# %% [markdown]
# # 3. select a group and retrieve the sequences

# %% [markdown]
# ## let's select the 'Vertebrata' group

# %%
ogid, oglevel = og_selection.select_OG_by_level_name(odb_gene_id=odb_gene_id, level_name='Vertebrata')
print(ogid)
print(oglevel)

# %%
# if we provide a level name that isn't found in the database, we get an error
og_selection.select_OG_by_level_name(odb_gene_id=odb_gene_id, level_name='Bacteria')

# %% [markdown]
# ## Get the sequences in the group

# %%
group_members = sql_queries.ogid_2_odb_gene_id_list(ogid)
sequence_dict = env.odb_database.get_sequences_from_list_of_seq_ids(group_members)
query_seqrecord = sequence_dict[odb_gene_id]
min_length = min_fraction_of_query_length * len(query_seqrecord)

# %%
print(query_seqrecord)

# %% [markdown]
# # 4. filter out sequences that are too short or contain non amino acid characters

# %%
import local_orthoDB_group_tools.filters as filters
print(len(sequence_dict))
filtered_sequence_dict = filters.filter_seqs_with_nonaa_chars(
    sequence_dict,
)
print(len(filtered_sequence_dict))
filtered_sequence_dict = filters.filter_shorter_sequences(
    filtered_sequence_dict,
    min_length=min_length,
)
if query_seqrecord.id not in filtered_sequence_dict:
    filtered_sequence_dict[query_seqrecord.id] = copy.deepcopy(query_seqrecord)

print(len(filtered_sequence_dict))

# %% [markdown]
# # 5. find least divergent orthologs (LDOs)

# %%
import local_orthoDB_group_tools.find_LDOs as find_LDOs

# %%
pid_df, ldos = find_LDOs.find_LDOs_main(
    seqrecord_dict=filtered_sequence_dict,
    query_seqrecord=query_seqrecord,
    pid_method = LDO_selection_method,
)

# %%
# sequence similarity matrix
pid_df.drop(columns=['sequence'])

# %%
ldo_seqrecord_dict = env.odb_database.get_sequences_from_list_of_seq_ids(ldos)

# %%
print(len(filtered_sequence_dict))
print(len(ldos))

# %% [markdown]
# # 6. cluster the LDOs using CD-HIT

# %%
import local_orthoDB_group_tools.cluster as cluster

# %%
clustered_ldo_seqrec_dict = cluster.cdhit_main(ldo_seqrecord_dict, query_seqrecord)

# %%
print(len(clustered_ldo_seqrec_dict))

# %% [markdown]
# # 7. align the clustered sequences

# %%
import local_seqtools.cli_wrappers as cli_wrappers

# %%
aln = cli_wrappers.mafft_align_wrapper(
    list(clustered_ldo_seqrec_dict.values()),
    output_type='alignment',
    n_align_threads=8,
)

# %%
print(aln[0:50, 175:220])

# %%

# %%

# %% [markdown]
# # putting it all together

# %%
import local_env_variables.env_variables as env
import json
from Bio import AlignIO, SeqIO
from pathlib import Path
from local_orthoDB_group_tools import (sql_queries, filters, uniprotid_search, og_selection, cluster, find_LDOs)
import local_seqtools.cli_wrappers as cli_wrappers
from local_config import conf
from attrs import asdict

# %%
processing_params = {
    "filter_params": {
        "min_fraction_short_than_query": 0.5,
    },
    "og_select_params": {
        "OG_selection_method": "level_name",
        "OG_level_name": 'Vertebrata',
    },
    "ldo_select_params": {
        "LDO_selection_method": "alfpy_google_distance",
    },
    "align_params": {
        "align": True,
        "n_align_threads": 8,
    },
    "main_output_folder": "./og_construction_output",
}

config = conf.pipeline_params.from_dict(processing_params)

# %%
uniprot_id = 'Q8TC90'
odb_gene_id = uniprotid_search.uniprotid_2_odb_gene_id(uniprot_id)


ogid, oglevel = og_selection.select_OG_by_level_name(
    odb_gene_id=odb_gene_id,
    level_name=config.og_select_params.OG_level_name
)
group_members = sql_queries.ogid_2_odb_gene_id_list(ogid)
sequence_dict = env.odb_database.get_sequences_from_list_of_seq_ids(group_members)
query_seqrecord = sequence_dict[odb_gene_id]


filtered_sequence_dict = filters.filter_seqs_with_nonaa_chars(
    sequence_dict,
)
min_length = config.filter_params.min_fraction_short_than_query* len(query_seqrecord)
filtered_sequence_dict = filters.filter_shorter_sequences(
    filtered_sequence_dict,
    min_length=min_length,
)
if query_seqrecord.id not in filtered_sequence_dict:
    filtered_sequence_dict[query_seqrecord.id] = copy.deepcopy(query_seqrecord)


pid_df, ldos = find_LDOs.find_LDOs_main(
    seqrecord_dict = filtered_sequence_dict,
    query_seqrecord = query_seqrecord,
    pid_method = config.ldo_select_params.LDO_selection_method,
)
ldo_seqrecord_dict = env.odb_database.get_sequences_from_list_of_seq_ids(ldos)


clustered_ldo_seqrec_dict = cluster.cdhit_main(ldo_seqrecord_dict, query_seqrecord)


if config.align_params.align:
    aln = cli_wrappers.mafft_align_wrapper(
        list(clustered_ldo_seqrec_dict.values()),
        output_type='alignment',
        n_align_threads=config.align_params.n_align_threads,
    )



# %%
alignment_output_file = Path(config.main_output_folder) / 'alignments' / f'{odb_gene_id.replace(":", "_")}_{oglevel}_{ogid}_clustered_ldos_aln.fasta'
alignment_output_file.parent.mkdir(parents=True, exist_ok=True)
with open(alignment_output_file, 'w') as f:
    AlignIO.write(aln, f, 'fasta')
print(alignment_output_file)

# %%
output_dict = {}
output_dict['query_uniprot_id'] = uniprot_id
output_dict['query_odb_gene_id'] = odb_gene_id
output_dict['query_sequence_str'] = str(query_seqrecord.seq)
output_dict['ogid'] = ogid
output_dict['oglevel'] = oglevel
output_dict['processing params'] = asdict(config)
output_dict['sequences'] = list(sequence_dict.keys())
output_dict['sequences_filtered'] = list(filtered_sequence_dict.keys())
output_dict['sequences_ldos'] = list(ldo_seqrecord_dict.keys())
output_dict['sequences_clustered_ldos'] = list(clustered_ldo_seqrec_dict.keys())
output_dict['alignment_clustered_ldos_file'] = str(alignment_output_file.resolve())
output_dict['alignment_clustered_ldos_file_relative'] = str(alignment_output_file.resolve().relative_to(Path.cwd()))

# %%
og_info_json_file = Path(config.main_output_folder) / 'info_jsons' / f'{odb_gene_id.replace(":", "_")}_{oglevel}_{ogid}_info.json'
og_info_json_file.parent.mkdir(parents=True, exist_ok=True)
with open(og_info_json_file, 'w') as f:
    json.dump(output_dict, f, indent=4)

# %%

# %%
output_dict = {}
output_dict['query_uniprot_id'] = uniprot_id
output_dict['query_odb_gene_id'] = odb_gene_id
output_dict['query_sequence_str'] = str(query_seqrecord.seq)
output_dict['ogs'] = {}
output_dict['ogs'][ogid] = {}
output_dict['ogs'][ogid]['oglevel'] = oglevel
output_dict['ogs'][ogid]['processing params'] = asdict(config)
output_dict['ogs'][ogid]['sequences'] = list(sequence_dict.keys())
output_dict['ogs'][ogid]['sequences_filtered'] = list(filtered_sequence_dict.keys())
output_dict['ogs'][ogid]['sequences_ldos'] = list(ldo_seqrecord_dict.keys())
output_dict['ogs'][ogid]['sequences_clustered_ldos'] = list(clustered_ldo_seqrec_dict.keys())
output_dict['ogs'][ogid]['alignment_clustered_ldos_file'] = str(alignment_output_file.resolve())
output_dict['ogs'][ogid]['alignment_clustered_ldos_file_relative'] = str(alignment_output_file.resolve().relative_to(Path.cwd()))

# %%
from pyprojroot import here

# %%
Path(config.main_output_folder).resolve().relative_to(here())

# %% [markdown]
# TODO:
# - [ ] figure out the file path thing. Probably have to convert to absolute paths
# - [ ] 

# %% [markdown]
# # 8. output files
#
# neccessary data to output:
# query protein:
# - uniprot id
# - odb gene id
# - species
#
# ortholog group:
# - og id
# - group sequence ids (odb gene ids):
#   - full group
#   - filtered group
#   - LDOs
#   - clustered LDOs
# - alignment file
#
# processing parameters:

# %%
print(asdict(config))

# %%
config

uniprot_id
odb_gene_id
query_seqrecord

ogid
sequence_dict
filtered_sequence_dict
ldo_seqrecord_dict
clustered_ldo_seqrec_dict

aln

# %%
