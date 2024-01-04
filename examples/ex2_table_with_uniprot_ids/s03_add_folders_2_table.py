import pandas as pd
import re
from pathlib import Path

# %%

table_path = './table.csv'
pipeline_output_path = './output/info_jsons/'
output_folder = Path(pipeline_output_path)

json_filelist = [f for f in output_folder.glob('*.json')]
json_filelist
# json_filelist = [f for f in ortho_group_folder.rglob('*_info___dict__.json')]
# json_filelist
# for json_file in json_filelist: print(json_file.parent.parent)
# for json_file in json_filelist: print(json_file.parent.parent.name)
# for json_file in json_filelist: print(json_file.name.replace('_info___dict__.json', ''))
# %%

p = re.compile(r'(.+)-.+_.+')

# for f in ortho_group_folder.glob('*/*/'): print(p.findall(f.name.split(':')[0])[0])
# for f in ortho_group_folder.glob('*/*/'): print(f.absolute())

output_folder_dict = {}
for f in ortho_group_folder.glob('*/'):
    uniprot_id = p.findall(f.name.split(':')[0])[0]
    if uniprot_id in output_folder_dict.keys():
        print(f'{uniprot_id} already in output_folder_dict. Duplicate uniprot_ids???')
    output_folder_dict[uniprot_id] = str(f.absolute())
output_folder_dict

# %%

# hits_df.columns
hits_df = pd.read_csv(hits_df_path)
hits_df['ortholog_group_folder'] = hits_df['UniprotID'].map(output_folder_dict)
hits_df_orthologs = hits_df.dropna(subset=['ortholog_group_folder']).copy()
# %%

# for f in ortho_group_folder.glob('*/params.json'): print(f)
# for json_file in json_filelist: print(json_file.name.replace('_info___dict__.json', ''))

# [f.name for f in Path('orthoDB_analysis_multiprocessed/msa_by_organism/Q9UPP2_9606_0:003212').glob('*/')]


def get_level_list(gene_id_folder: str) -> list:
    return [f.name for f in Path(gene_id_folder).glob('*/')]

# %%



# get_level_list(hits_df['ortholog_group_folder'].loc[0])



# %%

hits_df_orthologs['ortholog_levels'] = hits_df_orthologs['ortholog_group_folder'].apply(get_level_list)
hits_df_orthologs.to_csv(f'{Path(hits_df_path).stem}_with_orthologs.csv', index=False)

