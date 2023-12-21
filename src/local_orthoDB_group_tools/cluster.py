import copy

import pandas as pd

import local_seqtools.cdhit_tools as cdhit_tools
import local_seqtools.cli_wrappers as cli

# ==============================================================================
# // helper functions for CD-HIT
# ==============================================================================


def add_PID_2_seqid(seqrecord_dict_in, df_in):
    """
    run on the final set of sequences after clustering
    """
    seqrecord_dict = copy.deepcopy(seqrecord_dict_in)
    df = df_in.copy()
    for i in seqrecord_dict.values():
        score = df.loc[df["id"] == i.id, "PID"].values[0]
        i.id = f"{i.id}-{score:.2f}"
    # filtering the dataframe and sorting to sort the sequences by PID
    df = df[df["id"].isin(seqrecord_dict.keys())]
    df = df.sort_values("PID", ascending=False)
    sorted_seqrecord_list = [seqrecord_dict[i] for i in df["id"].values]
    return sorted_seqrecord_list


def cdhit_clstr_retrieve_representative_sequences(clstr_dict, seqrecord_dict):
    """
    pull out representative seqs defined in cdhit clstr_dict from full seqrecord_dict
    """
    clustered_seq_dict = {}
    for cluster_id in clstr_dict.keys():
        id_i = clstr_dict[cluster_id]["representative_seq"]
        # id_i = re.findall(r'\d+\_\d\:.+$', rep_i)[0]
        clustered_seq_dict[id_i] = seqrecord_dict[id_i]
    return clustered_seq_dict


def cdhit_minidriver(odbquery, seqrecords_2_cluster_list, linux=True, repr_id_keywords=None):
    """
    in: list of seqrecords
    out: dict of clustered seqrecords
    for getting LDOs first and then clustering
    """
    print(f"clustering LDOs using cdhit")
    print(f"{len(seqrecords_2_cluster_list)} sequences before clustering")

    if repr_id_keywords is None:
        repr_id_keywords = []
    if hasattr(odbquery, "query_gene_id"):
        repr_id_keywords.append(odbquery.query_gene_id)

    _, cdhit_clstr_dict = cli.cd_hit_wrapper(seqrecords_2_cluster_list, output_type="dict", linux=linux)
    cdhit_clstr_dict=cdhit_tools.cd_hit_clstr_redefine_cluster_representative_by_keywords(
        cdhit_clstr_dict, repr_id_keywords
    )

    # should be fine to pull the representative sequences from the full OG dictionary instead of creating a LDO dictionary
    seqrecords_2_cluster_dict = {
        seqrecord.id: seqrecord for seqrecord in seqrecords_2_cluster_list
    }
    sequences_clustered_OG_dict = cdhit_clstr_retrieve_representative_sequences(
        cdhit_clstr_dict, seqrecords_2_cluster_dict
    )
    print(f"{len(sequences_clustered_OG_dict)} sequences after clustering")
    return sequences_clustered_OG_dict


# ==============================================================================
# // DRIVER - clustering LDOs with cdhit
# ==============================================================================
def driver_CDHIT(odbquery, align=True, n_align_threads=8, linux=True):
    sequences_OG_LDO_cdhit_dict = cdhit_minidriver(
        odbquery, odbquery.sequences_LDO_list, linux=linux
    )
    odbquery.sequences_OG_LDO_cdhit_list = add_PID_2_seqid(
        sequences_OG_LDO_cdhit_dict, odbquery.pid2query_df
    )
    if align:
        odbquery.alignment_OG_LDO_cdhit = cli.mafft_align_wrapper(
            odbquery.sequences_OG_LDO_cdhit_list, n_align_threads=n_align_threads
        )
    print(
        f"number of sequences in clustered LDO set: {len(odbquery.sequences_OG_LDO_cdhit_list)}"
    )
