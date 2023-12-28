
import local_env_variables.env_variables as env
import local_orthoDB_group_tools.database as database
import local_orthoDB_group_tools.sql_queries as sql_queries
import pandas as pd

print(env.load_data_species_dict())




'''

odbdb = database.orthoDB_database()
for k,v in odbdb.data_species_dict.items():
    if 'sapiens' in v:
        print(k,v)


l = sql_queries.odb_gene_id_2_ogid_list('9606_0:001c7b')
df = sql_queries.ogid_list_2_og_df(l)
d = odbdb.data_levels_df[['level NCBI tax id', 'level name', 'total non-redundant count of species underneath']]

df2 = pd.merge(df, d, on='level NCBI tax id', how='left')
df2 = df2.astype({'total non-redundant count of species underneath': int})
df2 = df2.sort_values('total non-redundant count of species underneath', ascending=False)
print(df2[['OG id', 'level name', 'total non-redundant count of species underneath']].to_markdown())
'''