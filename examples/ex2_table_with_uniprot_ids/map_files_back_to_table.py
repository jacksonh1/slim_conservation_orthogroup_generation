# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.0
#   kernelspec:
#     display_name: odb_groups_x86
#     language: python
#     name: python3
# ---

# %%
import json
import re
from pathlib import Path

import pandas as pd

# %% [markdown]
# import table and specify pipeline output folder

# %%
table_path = './table.csv'
pipeline_output_path = './output/info_jsons/'
output_folder = Path(pipeline_output_path)

# %% [markdown]
# get list of json files in the folder

# %%
# .resolve() is used to get the absolute path
json_filelist = [f.resolve() for f in output_folder.glob('*.json')]


# %% [markdown]
# import each json file and extract the uniprot id. Use it to create a mapping between uniprot id and the json file path

# %%
def get_uniprot_id_from_json(json_file: Path) -> str:
    with open(json_file, 'r') as f:
        json_dict = json.load(f)
    return json_dict['query_uniprot_id']

json_file_map = {}
for file in json_filelist:
    print(get_uniprot_id_from_json(file), file)
    json_file_map[get_uniprot_id_from_json(file)] = file

# %% [markdown]
# map the file back to the table

# %%
table = pd.read_csv(table_path)
table.head()

# %%
table['ortholog group json'] = table['Uniprotid'].map(json_file_map)
table.head()

# %% [markdown]
# notice that it failed to find 2 of the uniprot ids in the database. That's because they are not in the sample dataset that I created for this repo

# %% [markdown]
# ---

# %% [markdown]
# you can easily add more information to the table from the json files<br>
# Let's try adding the number of sequences in the final ortholog group and the odb_gene_id for each protein

# %%
with open(json_filelist[0]) as f:
    test_json = json.load(f)
for i in test_json.keys(): print(i)


# %% [markdown]
# there are a lot of ways to do this but let's just repeat what we did above for the sake of time <br>
# in a real situation, you would probably want to think about how to do this in a more efficient way

# %%
def get_odb_gene_id_from_json(json_file: Path) -> str:
    with open(json_file, 'r') as f:
        json_dict = json.load(f)
    return json_dict['query_odb_gene_id']

def get_n_clustered_ldos_from_json(json_file: Path) -> str:
    with open(json_file, 'r') as f:
        json_dict = json.load(f)
    return len(json_dict['sequences_clustered_ldos'])


# %%
id_map = {k: get_odb_gene_id_from_json(v) for k, v in json_file_map.items()}
n_clustered_ldos_map = {k: get_n_clustered_ldos_from_json(v) for k, v in json_file_map.items()}
table['odb gene id'] = table['Uniprotid'].map(id_map)
table['n clustered ldos'] = table['Uniprotid'].map(n_clustered_ldos_map)

# %%
table.head()

# %%
table.to_csv('table_with_results.csv', index=False)
