import copy

import pandas as pd
from Bio import Seq, SeqIO


def filter_seqs_with_nonaa_chars(
    seqrecord_dict: dict[str, SeqIO.SeqRecord],
    prohibited_chars: list[str] = ["X", "x", "*"],
) -> dict[str, SeqIO.SeqRecord]:
    """
    filter sequences with non amino acid characters such as X and *.

    Returns a new dictionary with the filtered sequences
    """
    filtered_og_seq_dict = {}
    for seq_id, seq in seqrecord_dict.items():
        break_flag = False
        for char in prohibited_chars:
            if char in seq:
                break_flag = True
        if break_flag:
            continue
        filtered_og_seq_dict[seq_id] = copy.deepcopy(seq)
    return filtered_og_seq_dict


def filter_shorter_sequences(
    seqrecord_dict: dict[str, SeqIO.SeqRecord],
    min_length: int|float,
) -> dict[str, SeqIO.SeqRecord]:
    filtered_og_seq_dict = {}
    for seq_id, seq in seqrecord_dict.items():
        if len(seq) < min_length:
            continue        
        filtered_og_seq_dict[seq_id] = copy.deepcopy(seq)
    return filtered_og_seq_dict