from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Align, AlignIO, Seq, SeqIO
from pyprojroot import here

import local_seqtools.substitution_matrices as submats

root = here()
data = root / "data"
output = root / "src" / "local_substitution_matrices"


def prep_matrix_df_rownorm(mat_df, mat_name, output_dir="./matrices"):
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    mat_rownorm_df = submats.rownorm_matrix_df(mat_df)

    mat_df.to_csv(output_dir / f"{mat_name}.csv")
    mat_rownorm_df.to_csv(output_dir / f"{mat_name}_row_norm.csv")
    return mat_rownorm_df


def prep_matrix_df_offdiagnorm(mat_df, mat_name, output_dir="./matrices"):
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    mat_df = mat_df.copy()
    mat_df_new_diag = submats.matrixdf_diagonal_2_max_off_diagonal(mat_df)
    mat_df_norm = submats.normalize_matrix_df(mat_df)
    mat_df_new_diag_norm = submats.normalize_matrix_df(mat_df_new_diag)

    mat_df.to_csv(output_dir / f"{mat_name}.csv")
    mat_df_new_diag.to_csv(output_dir / f"{mat_name}_max_off_diagonal.csv")
    mat_df_norm.to_csv(output_dir / f"{mat_name}_norm.csv")
    mat_df_new_diag_norm.to_csv(output_dir / f"{mat_name}_max_off_diagonal_norm.csv")
    return mat_df_new_diag_norm


matrix_name = "BLOSUM62"
matBLOSUM62_df = submats.load_matrix_as_df(matrix_name)
matBLOSUM62_df_rownorm = prep_matrix_df_rownorm(matBLOSUM62_df, matrix_name, output_dir=output)
matBLOSUM62_offdiagmaxnorm_df = prep_matrix_df_offdiagnorm(matBLOSUM62_df, matrix_name, output_dir=output)

matrix_name = "EDSSMat50"
mat = Align.substitution_matrices.read(
    data / "2019-disorder-matrix/Matrices_and_Datasets/Matrices/EDSSMat50"
)
matEDSS50_df = submats.convert_matrix_array_2_df(mat)
matEDSS50_df_rownorm = prep_matrix_df_rownorm(matEDSS50_df, matrix_name, output_dir=output)
matEDSS50_offdiagmaxnorm_df = prep_matrix_df_offdiagnorm(matEDSS50_df, matrix_name, output_dir=output)
