# %%

import json

import local_orthoDB_analysis_tools_v3.odb_query_groups_data_classes as odb_query_groups_dataclasses

with open('./9606_0:000a0e_Eukaryota_info__dict__.json') as f:
    data = json.load(f)

querydata = odb_query_groups_dataclasses.odb_query_results.from_dict(data)
print(querydata.__repr__())
