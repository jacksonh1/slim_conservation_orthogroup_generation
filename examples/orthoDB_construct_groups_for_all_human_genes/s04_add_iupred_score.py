import json
import multiprocessing
import logging
from pathlib import Path
import time
import local_seqtools.general_utils as tools
import numpy as np
import pandas as pd
from Bio import Align, AlignIO, Seq, SeqIO

logging.basicConfig(filename='./s04.log', level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(name)s %(message)s')
logging.debug('Start of program')
logger=logging.getLogger(__name__)


iupred_cutoff=0.4
gap_lower_limit=10
idr_length_cutoff=8
ORTHOLOG_GROUP_FOLDER = './orthoDB_analysis_multiprocessed_minfrac-0.75/msa_by_organism/'
PARAMS = {}

logging.debug(iupred_cutoff)
logging.debug(gap_lower_limit)
logging.debug(idr_length_cutoff)
logging.debug(ORTHOLOG_GROUP_FOLDER)

# %%


def merge_idr_regions(idr_regions, gap_lower_limit=10):
    """
    if number of residues between two consecutive IDRs is less than or equal to this, then merge them
    """
    small_gap_merged_idr_regions = []
    if len(idr_regions) <= 1:
        return idr_regions
    region = idr_regions[0]
    for i in range(len(idr_regions)-1):
        gapi = idr_regions[i+1][0]-region[1]
        remaining=len(idr_regions)-i-2
        if gapi <= gap_lower_limit:
            region = [region[0], idr_regions[i+1][1]]
            if remaining == 0:
                small_gap_merged_idr_regions.append(region)
        else:
            small_gap_merged_idr_regions.append(region)
            region=idr_regions[i+1]
            if remaining == 0:
                small_gap_merged_idr_regions.append(region)
    return small_gap_merged_idr_regions


def len_filter_idr_regions(idr_regions, idr_length_cutoff=10):
    new_idr_regions = []
    for indexes in idr_regions:
        # print(indexes[1]-indexes[0])
        if indexes[1]-indexes[0] > idr_length_cutoff:
            new_idr_regions.append(indexes)
    return new_idr_regions


def driver_filter_idr_regions(sequence_str, iupred_cutoff=0.4, gap_lower_limit=10, idr_length_cutoff=8):
    idr_regions = tools.get_idr_indexes_untouched_indexing(tools.iupred_scores_from_seqstring(sequence_str), cutoff=iupred_cutoff)
    idr_regions = merge_idr_regions(idr_regions, gap_lower_limit=gap_lower_limit)
    idr_regions = len_filter_idr_regions(idr_regions, idr_length_cutoff=idr_length_cutoff)
    return idr_regions


def write_idr_regions_to_json(idr_regions: list, json_file_path, **kwargs):
    with open(json_file_path) as f:
        info_dict = json.load(f)
    info_dict['idr_regions'] = idr_regions
    info_dict['idr_definition_parameters'] = kwargs
    # info_dict.update(kwargs)
    with open(json_file_path, 'w') as f:
        json.dump(info_dict, f, indent=4)


def write_idrs_to_jsons_in_folder(folder_path, iupred_cutoff=0.4, gap_lower_limit=10, idr_length_cutoff=8):
    logger.debug(f'{folder_path} started')
    json_file_list = list(folder_path.rglob('*.json'))
    with open(json_file_list[0]) as f:
        info_dict = json.load(f)
    idr_regions = driver_filter_idr_regions(
        info_dict['query_sequence_str'],
        iupred_cutoff=iupred_cutoff,
        gap_lower_limit=gap_lower_limit,
        idr_length_cutoff=idr_length_cutoff
    )
    logger.debug(f'{folder_path} idr regions: {idr_regions}')
    for json_file_name in json_file_list:
        write_idr_regions_to_json(
            idr_regions,
            json_file_name,
            iupred_cutoff=iupred_cutoff,
            gap_lower_limit=gap_lower_limit,
            idr_length_cutoff=idr_length_cutoff
        )

# %%

def main():
    folders = list(Path(ORTHOLOG_GROUP_FOLDER).glob('*/'))
    # for folder in folders:
    #     write_idrs_to_jsons_in_folder(
    #         folder,
    #         iupred_cutoff=iupred_cutoff,
    #         gap_lower_limit=gap_lower_limit,
    #         idr_length_cutoff=idr_length_cutoff
    #     )
    p = multiprocessing.Pool(
        multiprocessing.cpu_count() - 1
    )
    f_args = [(i, iupred_cutoff, gap_lower_limit, idr_length_cutoff) for i in folders]
    p.starmap(write_idrs_to_jsons_in_folder, f_args)
    p.close()
    p.join()


if __name__ == '__main__':
    main()
    logger.debug('End of program')
