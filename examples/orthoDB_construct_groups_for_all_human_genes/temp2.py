# %%

import json

import local_orthoDB_analysis_tools_v3.odb_query_groups_pydantic as odb_query_groups_pydantic

with open('./orthoDB_analysis_multiprocessed_minfrac-0.75/msa_by_organism/9606_0:000103/Vertebrata/9606_0:000103_Vertebrata_info__dict__.json') as f:
    data = json.load(f)

querydata = odb_query_groups_pydantic.odb_query_results(**data)
print(querydata.alignment_orthologs_ldo_clustered)
print(querydata.model_dump_json(indent=4))
# print(json.dumps(querydata.model_json_schema(), indent=4))

