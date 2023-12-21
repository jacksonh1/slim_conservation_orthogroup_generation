from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Align, AlignIO, Seq, SeqIO

plt.style.use("custom_standard")
plt.style.use("custom_small")
import seaborn as sns

# from pyprojroot import here
# root = here()
# data = root / "data"
# matrices = root / "src" / "local_substitution_matrices"
# edssmat50_file = data / "2019-disorder-matrix/Matrices_and_Datasets/Matrices/EDSSMat50"

# BLOSUM62 = pd.read_csv(matrices / "BLOSUM62.csv", index_col=0)
# BLOSUM62_norm = pd.read_csv(matrices / "BLOSUM62_norm.csv", index_col=0)
# BLOSUM62_row_norm = pd.read_csv(matrices / "BLOSUM62_row_norm.csv", index_col=0)
# BLOSUM62_max_off_diagonal = pd.read_csv(
#     matrices / "BLOSUM62_max_off_diagonal.csv", index_col=0
# )
# BLOSUM62_max_off_diagonal_norm = pd.read_csv(
#     matrices / "BLOSUM62_max_off_diagonal_norm.csv", index_col=0
# )

# EDSSMat50 = pd.read_csv(matrices / "EDSSMat50.csv", index_col=0)
# EDSSMat50_norm = pd.read_csv(matrices / "EDSSMat50_norm.csv", index_col=0)
# EDSSMat50_row_norm = pd.read_csv(matrices / "EDSSMat50_row_norm.csv", index_col=0)
# EDSSMat50_max_off_diagonal = pd.read_csv(
#     matrices / "EDSSMat50_max_off_diagonal.csv", index_col=0
# )
# EDSSMat50_max_off_diagonal_norm = pd.read_csv(
#     matrices / "EDSSMat50_max_off_diagonal_norm.csv", index_col=0
# )

# MATRIX_DF_DICT = {
#     'BLOSUM62': matrices / "BLOSUM62.csv",
#     'BLOSUM62_norm': matrices / "BLOSUM62_norm.csv",
#     'BLOSUM62_row_norm': matrices / "BLOSUM62_row_norm.csv",
#     'BLOSUM62_max_off_diagonal': matrices / "BLOSUM62_max_off_diagonal.csv",
#     'BLOSUM62_max_off_diagonal_norm': matrices / "BLOSUM62_max_off_diagonal_norm.csv",
#     'EDSSMat50': matrices / "EDSSMat50.csv",
#     'EDSSMat50_norm': matrices / "EDSSMat50_norm.csv",
#     'EDSSMat50_row_norm': matrices / "EDSSMat50_row_norm.csv",
#     'EDSSMat50_max_off_diagonal': matrices / "EDSSMat50_max_off_diagonal.csv",
#     'EDSSMat50_max_off_diagonal_norm': matrices / "EDSSMat50_max_off_diagonal_norm.csv",
# }

# def load_precomputed_matrix_df(matrix_name=None):
#     '''
#     BLOSUM62
#     BLOSUM62_norm
#     BLOSUM62_row_norm
#     BLOSUM62_max_off_diagonal
#     BLOSUM62_max_off_diagonal_norm
#     EDSSMat50
#     EDSSMat50_norm
#     EDSSMat50_row_norm
#     EDSSMat50_max_off_diagonal
#     EDSSMat50_max_off_diagonal_norm
#     '''
#     if matrix_name is None:
#         print("Available matrices:")
#         for k in MATRIX_DF_DICT.keys():
#             print(k)
#     mat_df = pd.read_csv(MATRIX_DF_DICT[matrix_name], index_col=0)
#     return mat_df


def load_matrix_as_df(matrix_name):
    # convert to pandas dataframe
    mat = Align.substitution_matrices.load(matrix_name)
    AAs = sorted(
        [
            "V",
            "E",
            "K",
            "I",
            "H",
            "L",
            "G",
            "T",
            "M",
            "N",
            "S",
            "P",
            "A",
            "F",
            "W",
            "Y",
            "Q",
            "R",
            "C",
            "D",
        ]
    )
    AAs.extend(["B", "Z", "X", "*"])
    mat_df = pd.DataFrame(index=AAs, columns=AAs)
    for k in mat.keys():
        # skip = False
        # for i in ["B", "Z", "X", "*"]:
        #     if i in k:
        #         skip = True
        # if skip:
        #     continue
        mat_df.loc[k[0], k[1]] = mat[k]
    mat_df = mat_df.astype(float)
    return mat_df


def convert_matrix_array_2_df(matrix_array):
    # convert to pandas dataframe
    AAs = sorted(
        [
            "V",
            "E",
            "K",
            "I",
            "H",
            "L",
            "G",
            "T",
            "M",
            "N",
            "S",
            "P",
            "A",
            "F",
            "W",
            "Y",
            "Q",
            "R",
            "C",
            "D",
        ]
    )
    AAs.extend(["B", "Z", "X", "*"])
    mat_df = pd.DataFrame(index=AAs, columns=AAs)
    for k in matrix_array.keys():
        mat_df.loc[k[0], k[1]] = matrix_array[k]
    mat_df = mat_df.astype(float)
    return mat_df


def matrixdf_diagonal_2_max_off_diagonal(mat_df):
    AAs = sorted(
        [
            "V",
            "E",
            "K",
            "I",
            "H",
            "L",
            "G",
            "T",
            "M",
            "N",
            "S",
            "P",
            "A",
            "F",
            "W",
            "Y",
            "Q",
            "R",
            "C",
            "D",
        ]
    )
    AAs.extend(["B", "Z", "X", "*"])
    # set diagonal to min
    mat_df_newdiag = mat_df.copy()
    for i in AAs:
        mat_df_newdiag.loc[i, i] = mat_df.min().min()
    # set diagonal to max off-diagonal value
    max_off_diag = mat_df_newdiag.max().max()
    for i in AAs:
        mat_df_newdiag.loc[i, i] = max_off_diag
    return mat_df_newdiag


def normalize_matrix_df(mat_df):
    # normalize matrix
    mat_df_norm = (mat_df - mat_df.min().min()) / (
        mat_df.max().max() - mat_df.min().min()
    )
    return mat_df_norm


def rownorm_matrix_df(mat_df):
    mat_rownorm_df = mat_df.copy()
    for aa in mat_rownorm_df.index:
        row = mat_rownorm_df.loc[aa]
        mat_rownorm_df.loc[aa] = (row - row.min()) / (row.max() - row.min())
    return mat_rownorm_df


def plot_matrix(mat_df, title=None, ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 7))
    if title is not None:
        ax.set_title(title, fontsize=20)
    sns.heatmap(
        mat_df,
        annot=True,
        # cmap="coolwarm",
        cmap="vlag",
        fmt=".1g",
        # fmt=".2f",
        linewidths=0.5,
        ax=ax,
        annot_kws={"size": 8, "color": "black", "weight": "bold"},
        cbar=False,
        linecolor="black",
        clip_on=False,
        # vmin=0,
        # vmax=1,
    )
    return ax


def plot_and_save_matrix(filename, title=None):
    filename = Path(filename)
    mat_df = pd.read_csv(filename, index_col=0)
    if title is None:
        title = filename.stem
    plot_matrix(mat_df, title=title)
    plt.savefig(filename.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close()
