import copy

import pandas as pd

import local_orthoDB_tools.database_v6 as odb_v6
import local_seqtools.alignment_tools as aln_tools
import local_seqtools.cli_wrappers as cli


def setup_df(seqrecord_dict_in):
    seqrecord_dict = copy.deepcopy(seqrecord_dict_in)
    seqrecord_list = [seq for seq in seqrecord_dict.values()]
    df = pd.DataFrame(columns=["id", "organism"], index=range(len(seqrecord_dict)))
    df["id"] = [seqrecord.id for seqrecord in seqrecord_list]
    df["organism"] = df["id"].str.split("|", expand=True)[0]
    df["sequence"] = df["id"].map(seqrecord_dict)
    return df, seqrecord_dict


def _alfpy_query_matrix(query_sequence, matrix):
    """
    get the row from the alfpy matrix that corresponds to the query gene
    """
    query_row = [c for c, i in enumerate(matrix.id_list) if query_sequence.id in i][0]
    query_row_distance = matrix.data[query_row]
    query_row_similarity = 1 - query_row_distance
    return matrix.id_list, query_row_similarity


def addpid_by_msa(
    df_in, query_sequence, seqrecord_dict, fast_msa=False, n_align_threads=8
):
    df = df_in.copy()
    seqrecord_list = [seq for seq in seqrecord_dict.values()]
    msa_seqrecord_list = cli.mafft_align_wrapper(
        seqrecord_list, fast=fast_msa, n_align_threads=n_align_threads
    )
    query_msa_seqrecord = [i for i in msa_seqrecord_list if i.id == query_sequence.id][
        0
    ]
    # could do this by applying a function to the `id` column but this seems a bit simpler
    df["PID"] = [
        aln_tools.compute_pairwise_percent_id_from_msa(
            query_msa_seqrecord,
            seqrecord,
        )
        for seqrecord in msa_seqrecord_list
    ]
    return df


def addpid_by_msa_by_organism(df_in, query_sequence, fast_msa=False, n_align_threads=8):
    df = df_in.copy()
    print("aligning sequences using mafft, one organism at a time")
    org_set = df["organism"].unique()
    pid_map_dict = {}
    # counter = 0
    for org_i in org_set:
        # counter += 1
        # get all the sequences for this organism
        seqs = list(df[df["organism"] == org_i]["sequence"].values)
        # add the query sequence to the list
        if query_sequence.id not in [seq.id for seq in seqs]:
            seqs.append(query_sequence)
        # align the sequences
        msa_i_dict = cli.mafft_align_wrapper(
            seqs, output_type="dict", fast=fast_msa, n_align_threads=n_align_threads
        )
        # compute the pairwise percent identity from the MSA
        for seq_id, seqrecord in msa_i_dict.items():
            pid_map_dict[seq_id] = aln_tools.compute_pairwise_percent_id_from_msa(
                seqrecord, msa_i_dict[query_sequence.id]
            )
        # if counter % 100 == 0:
        # print(f'finished {counter} of {len(org_set)} organisms')
    df["PID"] = df["id"].map(pid_map_dict)
    return df


def addpid_by_alfpy_google_distance(df_in, query_sequence):
    df = df_in.copy()
    print("comparing sequences using alignment free comparison (alfpy google distance)")
    org_set = df["organism"].unique()
    pid_map_dict = {}
    # counter = 0
    for org_i in org_set:
        # counter += 1
        # get all the sequences for this organism
        seqs = list(df[df["organism"] == org_i]["sequence"].values)
        # add the query sequence to the list
        if query_sequence.id not in [seq.id for seq in seqs]:
            seqs.append(query_sequence)
        matrix = aln_tools.alfpy_distance_matrix(seqs, word_size=2)
        id_list, query_row_similarity = _alfpy_query_matrix(query_sequence, matrix)
        for seq_id, similarity in zip(id_list, query_row_similarity):
            pid_map_dict[seq_id] = similarity
        # if counter % 100 == 0:
        # print(f'finished {counter} of {len(org_set)} organisms')
    df["PID"] = df["id"].map(pid_map_dict)
    return df


def addpid_by_pairwise(df_in, query_sequence, seqrecord_dict):
    df = df_in.copy()
    seqrecord_list = [seq for seq in seqrecord_dict.values()]
    print("aligning sequences using pairwise alignment")
    df["PID"] = [
        aln_tools.align_and_get_similarity(
            query_sequence,
            seqrecord,
            scoring_matrix_name=None,
            gap_opening_penalty=0,
            gap_extension_penalty=0,
        )
        for seqrecord in seqrecord_list
    ]
    return df


def _wrap_up(df, seqrecord_dict, odbquery: odb_v6.orthoDB_query):
    # remove sequences in the query organism that are not the query sequence
    df = df[(df["organism"] != odbquery.query_species_name.replace(" ", "_")) | (df["id"] == odbquery.query_sequence_id_str)]
    assert odbquery.query_sequence_id_str in df["id"].values, "query sequence not found in df"
    # groupby select the closest sequence for each organism
    ldo_df = df.loc[df.groupby("organism")["PID"].idxmax()].copy()
    ldo_df = ldo_df.sort_values("PID", ascending=False)
    ldo_seqrecord_list = [seqrecord_dict[i] for i in ldo_df["id"].values]
    df = df.drop("sequence", axis=1)
    odbquery.pid2query_df = df
    odbquery.sequences_LDO_list = ldo_seqrecord_list


def main(odbquery, pid_method="msa_by_organism", fast_msa=False, n_align_threads=8):
    assert pid_method in [
        "msa_by_organism",
        "alfpy_google_distance",
        "pairwise",
        "msa",
    ], "LDO selection method not recognized. must be one of: msa_by_organism, alfpy_google_distance, pairwise, msa"
    df, seqrecord_dict = setup_df(odbquery.sequences_full_OG_dict)
    if pid_method == "msa_by_organism":
        df = addpid_by_msa_by_organism(
            df, odbquery.query_sequence, n_align_threads=n_align_threads
        )
    elif pid_method == "alfpy_google_distance":
        df = addpid_by_alfpy_google_distance(df, odbquery.query_sequence)
    elif pid_method == "pairwise":
        df = addpid_by_pairwise(
            df, odbquery.query_sequence, seqrecord_dict
        )
    elif pid_method == "msa":
        df = addpid_by_msa(
            df,
            odbquery.query_sequence,
            seqrecord_dict,
            fast_msa=fast_msa,
            n_align_threads=n_align_threads,
        )
    _wrap_up(df, seqrecord_dict, odbquery)
