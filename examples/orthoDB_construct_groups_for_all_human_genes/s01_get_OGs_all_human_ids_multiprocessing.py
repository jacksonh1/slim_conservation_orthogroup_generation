import multiprocessing

import numpy as np
import pandas as pd
from local_orthoDB_group_tools import main_pipeline
from local_orthoDB_group_tools import sql_queries

def multiple_levels(query_geneid):
    params = {
        "main_output_folder": "./orthoDB_analysis_multiprocessed_minfrac-0.75",
        "OG selection method": "level name",
        "min_fraction_short_than_query": "0.75",
        "LDO selection method": "alfpy_google_distance",
        "align": "True",
        "n_align_threads": "6",
        "write files": "True"
    }
    og_levels = ['Eukaryota', 'Mammalia', 'Metazoa', 'Tetrapoda', 'Vertebrata']
    assert params["OG selection method"] == "level name", "OG selection method must be level name"
    
    for og_level in og_levels:
        params['OG level name'] = og_level
        try:
            main_pipeline.pipeline_setID_directly(query_geneid, user_params=params, linux=False)
        except ValueError as err:
            # logger.error(f"{query_geneid} - {og_level} - {err}")
            print(f"{query_geneid} - {og_level} - {err}")


# for k,v in odbdb.data_species_dict.items():
#     if 'sapiens' in v:
#         print(k,v)

def main():
    geneid_list = sql_queries.get_all_odbids_from_species_id('9606_0')
    # geneid_list = ["9606_0:00194d", "9606_0:002f40"]
    # p = multiprocessing.Pool(3)
    # not sure if map is the best choice vs map_async or imap etc...
    # p.map(multiple_levels, gene_subset)
    # p.close()
    # p.join()
    for i in geneid_list:
        multiple_levels(i)


if __name__ == '__main__':
    main()
