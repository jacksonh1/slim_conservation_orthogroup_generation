import logging
import multiprocessing

import numpy as np
import pandas as pd
from local_orthoDB_tools import database_v6, main_pipeline

logging.basicConfig(filename='./s01_error.log', level=logging.DEBUG, 
                    format='%(asctime)s %(levelname)s %(name)s %(message)s')
logging.debug('Start of program')
logger=logging.getLogger(__name__)

def multiple_levels(query_geneid):
    params = {
        "main_output_folder": "./orthoDB_analysis_multiprocessed_minfrac-0.75",
        "OG selection method": "level name",
        "min_fraction_short_than_query": "0.75",
        "LDO selection method": "msa_by_organism",
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
            logger.error(f"{query_geneid} - {og_level} - {err}")

def main():
    geneid_list = list(database_v6.orthoDB_query.odb_database.data_geneid_2_og_list_dict.keys())
    # p = multiprocessing.Pool(3)
    # not sure if map is the best choice vs map_async or imap etc...
    # for this example, just run on a randomly selected 3 genes
    # p.map(multiple_levels, gene_subset)
    # p.close()
    # p.join()
    gene_subset = np.random.choice(geneid_list, 3)
    for i in gene_subset:
        multiple_levels(i)


if __name__ == '__main__':
    main()
