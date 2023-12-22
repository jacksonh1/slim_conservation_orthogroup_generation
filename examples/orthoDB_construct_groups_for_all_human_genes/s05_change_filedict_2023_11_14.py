# %%
import json
from pathlib import Path

import local_env_variables.env_variables_and_filepaths as fp
from pyprojroot import here

ODBQUERY_RUN_ROOT = fp.orthodb_query_groups_root
# ODBQUERY_RUN_ROOT = Path('/Users/jackson/Dropbox (MIT)/work/07-SLiM_bioinformatics/01/pipeline_scripts/orthoDB_construct_groups_for_all_human_genes')
ORTHOLOG_GROUP_FOLDER = Path('./orthoDB_analysis_multiprocessed_minfrac-0.75/msa_by_organism/')
RN_DICT = {
    "fasta sequences - full OG": "file_sequences_full_OG",
    "fasta sequences - OG LDOs": "file_sequences_LDOs",
    "fasta sequences - OG LDOs cdhit": "file_sequences_orthologs_ldo_clustered",
    "fasta alignment - OG LDO cdhit": "file_alignment_orthologs_ldo_clustered",
    "json - query and ortho group info": "file_query_info_json",
}

def rename_dict_keys(d, rn_dict):
    return {rn_dict[i]: d[i] for i in d.keys() if i in rn_dict.keys()}

def convert_relative_to_root(file_dict_abs, root):
    return {i: str(Path(file_dict_abs[i]).relative_to(root)) for i in file_dict_abs.keys()}

def add_new_filedict(ogquery_info_dict, rn_dict, root=ODBQUERY_RUN_ROOT):
    file_dict_abs = ogquery_info_dict['output_file_dict_absolute']
    file_dict_rel = convert_relative_to_root(file_dict_abs, root)
    ogquery_info_dict['output_file_dict_rel'] = rename_dict_keys(file_dict_rel, rn_dict)
    return ogquery_info_dict


def main():
    json_files = list(ORTHOLOG_GROUP_FOLDER.glob('*/*/*.json'))
    for json_file in json_files:
        print(json_file)
        with open(json_file) as f:
            data = json.load(f)
        if 'output_file_dict_absolute' in data.keys():
            data = add_new_filedict(data, RN_DICT, root=ODBQUERY_RUN_ROOT)
            with open(json_file, 'w') as f:
                json.dump(data, f, indent=4)
        else:
            print(f'No output_file_dict_absolute key in json file: {json_file}')

if __name__ == '__main__':
    main()
