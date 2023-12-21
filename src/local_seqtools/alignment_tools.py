import pandas as pd
from alfpy import word_distance, word_pattern, word_vector
from alfpy.utils import distmatrix
from alfpy.utils import seqrecords as alf_seqrecords
from Bio import Align, AlignIO, Seq, SeqIO

import local_env_variables.matrices as mats
import local_seqtools.general_utils as tools

# you could preload the matrices here
matrix_dict = {i:mats.load_precomputed_matrix_df(i) for i in mats.MATRIX_DF_DICT.keys()}

def score_alignment(seq1: str, seq2: str, subs_mat_df: pd.DataFrame, gap_open=10, gap_extend=0.5):
    """returns the score of an alignment between two sequences

    treats the following situation as opening 2 gaps:
        AAAA---
        ----BBB
    
    you can use this function to compute a similarity score between two sequences
    if you normalize the score by the max possible score which might be something like:
    >>> max_score = max(score_alignment(seq1, seq1), score_alignment(seq2, seq2))

    so the similarity score would be:
    >>> max_score = max(score_alignment(seq1, seq1), score_alignment(seq2, seq2))
    >>> similarity_score = score_alignment(seq1, seq2) / max_score

    Parameters
    ----------
    seq1 : str
        first sequence
    seq2 : str
        second sequence
    subs_mat_df : pd.DataFrame
        the scoring matrix as a pandas dataframe
    gap_open : int, optional
        penalty for opening a gap (will be subtracted from score), by default 10
    gap_extend : int, optional
        penalty for extending a gap (will be subtracted from score), by default 0.5

    Returns
    -------
    float
        the score of the alignment
    """    
    assert len(seq1) == len(seq2)
    score = 0
    gap1_open_flag = False
    gap2_open_flag = False
    for s1, s2 in zip(seq1, seq2):
        if s1 == '-':
            if gap1_open_flag:
                score += gap_extend
            else:
                score += gap_open
                gap1_open_flag = True
        elif s2 == '-':
            if gap2_open_flag:
                score += gap_extend
            else:
                score += gap_open
                gap2_open_flag = True
        else:
            score += subs_mat_df.loc[s1, s2]
            gap1_open_flag = False
            gap2_open_flag = False
    return score


def score_alignment_from_alignment_obj(alignment_obj, subs_mat_df, gap_open, gap_extend):
    seq1 = alignment_obj[0]
    seq2 = alignment_obj[1]
    return score_alignment(seq1, seq2, subs_mat_df, gap_open, gap_extend)


def slice_alignment_by_subseq(alignment, hit_subseq):
    '''
    ex:
    alignment:
        target          560 LPGPQITFPPPPPPPVDD-----SPPDFLPPPPPAANFGSHPPPP 600
                          0 ||------|||||||.||-----.||||..|||---.|...|||.  45
        query             0 LP------PPPPPPPLDDPELPPPPPDFMEPPP---DFVPPPPPS  36
    hit_subseq:
        LDDPELPPPPPDFMEPPP
    
    returns:
        LDDPELPPPPPDFMEPPP
        .||-----.||||..|||
        VDD-----SPPDFLPPPP
    '''
    _, index = tools.reindex_alignment_str(alignment[1])
    unaligned_hit = alignment[1].replace("-", "")
    st = unaligned_hit.find(hit_subseq)
    end = st + len(hit_subseq) - 1
    aln_st = index[st]
    aln_end = index[end]
    return alignment[aln_st:aln_end+1]


def pairwise_alignment(seq1, seq2, scoring_matrix_name = 'edssmat50', gap_opening_penalty = 10, gap_extension_penalty = 0.5):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    if scoring_matrix_name=='edssmat50':
        aligner.substitution_matrix = Align.substitution_matrices.read(mats.edssmat50_file)
    else:
        aligner.substitution_matrix = Align.substitution_matrices.load(scoring_matrix_name)
    aligner.extend_gap_score = -gap_extension_penalty
    aligner.open_gap_score = -gap_opening_penalty
    alignment=aligner.align(seq1, seq2)[0]
    return alignment


def compute_pairwise_percent_id_from_msa(seqrecord1, seqrecord2):
    """compute pairwise percent identity between two sequences
    - The sequences must be pre-aligned (i.e. they have the same length)
    - The returned percent identity is computed as the number of identical 
    residues divided by the LENGTH OF compared alignment sites. Alignment sites 
    where both sequences are gaps ('-' characters) are not compared.  I am hoping 
    that this will provide some sort of global measure of similarity between 
    two sequences as opposed to local similarity
    - Gaps are not counted as matches
    TODO: add a parameter to allow for substitution matrix
    TODO: do a pairwise alignment for each sequence individually and then compute percent identity
    

    Parameters
    ----------
    seqrecord1 : seqrecord
        first sequence
    seqrecord2 : seqrecord
        second sequence

    Returns
    -------
    float
        percent identity between two sequences.
    """    
    seq1 = str(seqrecord1.seq)
    seq2 = str(seqrecord2.seq)
    # print(seq1)
    # print(seq2)
    assert len(seq1) == len(seq2), 'sequences are not the same length'
    num_same = 0
    length = len(seq1)
    for i in range(len(seq1)):
        if seq1[i] == '-' and seq2[i] == '-':
            length -= 1
            continue
        if seq1[i] == seq2[i]:
            num_same += 1
    return num_same / length


def __normalized_pairwise_similarity(seq1, seq2, scoring_matrix_name = 'BLOSUM62', gap_opening_penalty = 10, gap_extension_penalty = 0.5):
    """
    Function is modified from the `_pairwise` function in Biopython (Bio.Phylo.TreeConstruction.DistanceCalculator)
    TODO: insert citation or copy right information here

    Calculate pairwise similarity from seq1 to seq2.
    I am specifically considering the distance between seq1 and seq2 and not between seq2 and seq1
    B/c it matters when you use matrices that do not have constant values across the diagonal
    (e.g. BLOSUM62), since we will be dividing the score for each pair in the alignment by the score
    of the residue in seq1 being paired with itself (i.e. conserved).

    See which matrices are available with:
    >>> from Bio import Align
    >>> print(Align.substitution_matrices.load())

    Returns a value between 0 (identical sequences) and 1 (completely
    different, or seq1 is an empty string.)
    
    Parameters
    ----------
    seq1 : str or seqrecord
        first sequence
    seq2 : str or seqrecord
        second sequence
    matrix : str, optional
        substitution matrix to use for the scoring, by default it will use BLOSUM62. If None, then it will calculate percent identity
    """
    # query_length = len(seq1)
    assert len(seq1) == len(seq2), 'sequences are not the same length. Are they aligned?'
    if scoring_matrix_name is None:
        # calculate percent identity
        num_same = 0
        length = len(seq1)
        for i in range(len(seq1)):
            if seq1[i] == '-' and seq2[i] == '-':
                length -= 1
                continue
            if seq1[i] == seq2[i]:
                num_same += 1
        norm_score = num_same / length
        # norm_score = compute_pairwise_percent_id_from_msa(seq1, seq2)
    else:
        score = 0
        max_score = 0
        gap1open = False
        gap2open = False
        scoring_matrix=Align.substitution_matrices.load(scoring_matrix_name)
        # calculate just the score of the aligned residues
        for i in range(0, len(seq1)):
            s1 = seq1[i]
            s2 = seq2[i]
            # if both residues are gaps, then skip it
            if seq1[i] == '-' and seq2[i] == '-':
                continue
            # if either residue is a gap, then add the gap penalty
            # I used 2 different gaps to deal with this situation:
            # AAAA---
            # ----BBB
            # i would say that this is 2 gaps?
            # Here's a version where it treats this situation as 2 gaps:
            # if s1 == '-':
            #     if gap1open:
            #         score -= gap_extension_penalty
            #     else:
            #         score -= gap_opening_penalty
            #         gap1open = True
            #     continue
            # # FYI - since I used a continue startement, I don't need an else statement
            # else:
            #     gap1open = False
            # if s2 == '-':
            #     if gap2open:
            #         score -= gap_extension_penalty
            #     else:
            #         score -= gap_opening_penalty
            #         gap2open = True
            #     continue
            # else:
            #     gap2open = False
            # here's a version where it treats this situation:
            # AAAA---
            # ----BBB 
            # as 1 gap:
            if s1 == '-' or s2 == '-':
                if gap1open:
                    score -= gap_extension_penalty
                else:
                    score -= gap_opening_penalty
                    gap1open = True
                continue
            else:
                gap1open = False
                # if neither residue is a gap, then add the score for the pair
                score += scoring_matrix[s1, s2]
        # calculate the max score for the alignment
        for a in seq1:
            if a == '-':
                continue
            max_score += scoring_matrix[a, a]
        # for b in seq2:
        #     if b == '-':
        #         continue
        #     max_score += scoring_matrix[b, b]
        # you could calculate the max score for each sequence and then take the higher one
        norm_score = score / max_score
    return norm_score


def __align_and_get_similarity(
    seqrecord1,
    seqrecord2,
    scoring_matrix_name=None,
    gap_opening_penalty=0,
    gap_extension_penalty=0,
):
    '''
    DEPRECATED
    aligns the two sequences with Biopython's PairwiseAligner using the BLOSUM62 scoring matrix and gap_opening_penalty = 10 and gap_extension_penalty = 0.5. These are hard coded for now.

    The `scoring_matix_name`, `gap_opening_penalty`, and `gap_extension_penalty` are passed to the `_normalized_pairwise_similarity` function to be used only in calculating similarity from the alignments. It just does percent identity by default.
    '''
    aln = pairwise_alignment(seqrecord1.seq, seqrecord2.seq)
    pid = _normalized_pairwise_similarity(
        aln[0],
        aln[1],
        scoring_matrix_name=scoring_matrix_name,
        gap_opening_penalty=gap_opening_penalty,
        gap_extension_penalty=gap_extension_penalty,
    )
    return pid


def __perc_ID_matrix_from_seqrecord_list(seqrecord_list, **kwargs):
    '''
    **kwargs are passed to the `align_and_get_similarity` function
    '''
    id_matrix = pd.DataFrame(
        columns = [i.id for i in seqrecord_list],
        index = [i.id for i in seqrecord_list]
    )
    for i in seqrecord_list:
        for j in seqrecord_list:
            id_matrix.loc[i.id, j.id] = align_and_get_similarity(i, j, **kwargs)
    id_matrix = id_matrix.infer_objects()
    return id_matrix


def alfpy_distance_matrix(seqrecord_list, word_size=2):
    id_list = [i.id for i in seqrecord_list]
    seq_str_list = [str(i.seq) for i in seqrecord_list]
    alf_seq_records = alf_seqrecords.SeqRecords(id_list=id_list, seq_list=seq_str_list)
    p = word_pattern.create(alf_seq_records.seq_list, word_size=word_size)
    # counts = word_vector.Counts(alf_seq_records.length_list, p)
    freqs = word_vector.Freqs(alf_seq_records.length_list, p) # why is that not used?
    dist = word_distance.Distance(freqs, 'google')
    matrix = distmatrix.create(alf_seq_records.id_list, dist)
    return matrix

