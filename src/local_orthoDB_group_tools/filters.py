import copy

import pandas as pd


def filter_sequences_with_X(odbquery):
    '''
    filter sequences with X character
    '''
    xi = len(odbquery.sequences_full_OG_dict)
    filtered_og_seq_dict = {}
    for seq_id, seq in odbquery.sequences_full_OG_dict.items():
        if 'X' not in seq and '*' not in seq:
            filtered_og_seq_dict[seq_id] = seq
    if odbquery.query_sequence_id_str not in filtered_og_seq_dict:
        print(f'query sequence {odbquery.query_sequence_id_str} removed by filter_sequences_with_X')
        print(f'adding back query sequence {odbquery.query_sequence_id_str} to filtered_og_seq_dict')
        filtered_og_seq_dict[odbquery.query_sequence_id_str] = copy.deepcopy(odbquery.sequences_full_OG_dict[odbquery.query_sequence_id_str])
    xf = len(filtered_og_seq_dict)
    odbquery.filter_dict['filter - sequences with *s and Xs'] = {
        'number of sequences before filtering':xi,
        'number of sequences removed':xi-xf,
        'number of sequences after filtering':xf,
    }
    odbquery.sequences_full_OG_dict = filtered_og_seq_dict

def filter_sequences_by_length(odbquery, min_fraction_short_than_query=0.5):
    li=len(odbquery.sequences_full_OG_dict)
    filtered_og_seq_dict = {}
    for odb_gene_id, seq in odbquery.sequences_full_OG_dict.items():
        if len(seq) >= min_fraction_short_than_query * len(odbquery.query_sequence):
            filtered_og_seq_dict[odb_gene_id] = seq
    lf=len(filtered_og_seq_dict)
    odbquery.filter_dict['filter - sequence length relative to query length'] = {
        'min fraction of query sequence length':min_fraction_short_than_query,
        'number of sequences before filtering':li,
        'number of sequences removed':li-lf,
        'number of sequences after filtering':lf,
    }
    odbquery.sequences_full_OG_dict = filtered_og_seq_dict