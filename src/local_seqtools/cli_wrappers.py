import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import Align, AlignIO, Seq, SeqIO

import local_env_variables.env_variables as env
import local_seqtools.cdhit_tools as cdhit_tools
import local_seqtools.general_utils as tools


def mafft_align_wrapper(
    input_seqrecord_list: list[SeqIO.SeqRecord],
    mafft_executable: str = env.MAFFT_EXECUTABLE,
    extra_args: str = env.MAFFT_ADDITIONAL_ARGUMENTS,
    output_type: str = "list",
    fast: bool = False,
    n_align_threads: int = 8,
) -> tuple[str, list[SeqIO.SeqRecord] | dict[str, SeqIO.SeqRecord] | AlignIO.MultipleSeqAlignment]:
    assert output_type in [
        "list",
        "dict",
        "alignment",
    ], f'`output_type` must be one of ["list", "alignment"]'
    # create temporary file
    temp_file = tempfile.NamedTemporaryFile(mode="w", delete=False)
    # write seqrecords to temporary file
    SeqIO.write(input_seqrecord_list, temp_file, "fasta")
    temp_file.close()
    # run mafft
    alignment_filename = f"{temp_file.name}-mafft.fa"
    # raise an error if the alignment file already exists. (it won't but just in case)
    if os.path.exists(alignment_filename):
        raise FileExistsError(f"{alignment_filename} already exists")
    if fast:
        mafft_command = f'{mafft_executable} --thread {n_align_threads} --quiet --retree 1 {extra_args} "{temp_file.name}" > "{alignment_filename}"'
    else:
        mafft_command = f'{mafft_executable} --thread {n_align_threads} --quiet --anysymbol {extra_args} "{temp_file.name}" > "{alignment_filename}"'
    # print(mafft_command)
    subprocess.run(mafft_command, shell=True, check=True)
    # read in mafft output
    if output_type == "list":
        mafft_output = tools.import_fasta(alignment_filename, output_format="list")
    elif output_type == "dict":
        mafft_output = tools.import_fasta(alignment_filename, output_format="dict")
    elif output_type == "alignment":
        mafft_output = AlignIO.read(alignment_filename, "fasta")
    # delete temporary file
    os.remove(alignment_filename)
    os.remove(temp_file.name)
    return mafft_command, mafft_output


def cd_hit_wrapper(
    input_seqrecord_list: list[SeqIO.SeqRecord],
    cd_hit_executable: str = env.CD_HIT_EXECUTABLE,
    extra_args: str = env.CD_HIT_ADDITIONAL_ARGUMENTS,
    output_type="list",
) -> tuple[str, list[SeqIO.SeqRecord] | dict[str, SeqIO.SeqRecord], dict[str, dict[str, str]]]:
    assert output_type in [
        "list",
        "dict",
    ], f'`output_type` must be one of ["list", "dict"]'

    # create temporary file
    temp_file = tempfile.NamedTemporaryFile(mode="w", delete=False)
    # write seqrecords to temporary file
    SeqIO.write(input_seqrecord_list, temp_file, "fasta")
    temp_file.close()
    clustered_seqs_filename = f"{temp_file.name}-cdhit.fa"
    # raise an error if the alignment file already exists. (it won't but just in case)
    if os.path.exists(clustered_seqs_filename):
        raise FileExistsError(f"{clustered_seqs_filename} already exists")

    clustered_seqs_clusters_filename = clustered_seqs_filename + ".clstr"
    command = f"{cd_hit_executable} -i {temp_file.name} -o {clustered_seqs_filename} -M 0 -d 0 {extra_args}"
    subprocess.run(command, shell=True, check=True)

    output_clstrs_dict = cdhit_tools.cd_hit_clstr_parser(
        clustered_seqs_clusters_filename
    )

    # read in clustal output
    if output_type == "list":
        output = tools.import_fasta(clustered_seqs_filename, output_format="list")
    elif output_type == "dict":
        output = tools.import_fasta(clustered_seqs_filename, output_format="dict")
    # delete temporary file
    os.remove(clustered_seqs_filename)
    os.remove(clustered_seqs_clusters_filename)
    os.remove(temp_file.name)
    return command, output, output_clstrs_dict



def clustal_align_wrapper(
    input_seqrecord_list, alignment_type="basic", output_type="list"
):
    assert output_type in [
        "list",
        "dict",
        "alignment",
    ], f'`output_type` must be one of ["list", "dict", "alignment"]'
    assert alignment_type in [
        "basic",
        "full",
    ], f'`output_type` must be one of ["basic", "full"]'
    # create temporary file
    temp_file = tempfile.NamedTemporaryFile(mode="w", delete=False)
    # write seqrecords to temporary file
    SeqIO.write(input_seqrecord_list, temp_file, "fasta")
    temp_file.close()
    alignment_filename = f"{temp_file.name}-clustal.fa"
    # raise an error if the alignment file already exists. (it won't but just in case)
    if os.path.exists(alignment_filename):
        raise FileExistsError(f"{alignment_filename} already exists")

    if alignment_type == "basic":
        clustal_command = f'clustalo -i "{temp_file.name}" -o "{alignment_filename}" -v --outfmt=fa --threads=6'
    elif alignment_type == "full":
        clustal_command = f'clustalo -i "{temp_file.name}" -o "{alignment_filename}" -v --outfmt=fa --full --threads=6'
    subprocess.run(clustal_command, shell=True, check=True)

    # read in clustal output
    if output_type == "list":
        clustal_output = tools.import_fasta(alignment_filename, output_format="list")
    elif output_type == "dict":
        clustal_output = tools.import_fasta(alignment_filename, output_format="dict")
    elif output_type == "alignment":
        clustal_output = AlignIO.read(alignment_filename, "fasta")
    # delete temporary file
    os.remove(alignment_filename)
    os.remove(temp_file.name)
    return clustal_output

        
def muscle_align_wrapper(
    input_seqrecord_list: list[SeqIO.SeqRecord],
    muscle_binary: str = "/Users/jackson/tools/muscle/muscle-5.1.0/src/Darwin/muscle",
    output_type:str = "list"
) -> list[SeqIO.SeqRecord] | dict[str, SeqIO.SeqRecord] | AlignIO.MultipleSeqAlignment:
    assert output_type in [
        "list",
        "dict",
        "alignment",
    ], f'`output_type` must be one of ["list", "dict", "alignment"]'

    temp_file = tempfile.NamedTemporaryFile(mode="w", delete=False)
    SeqIO.write(input_seqrecord_list, temp_file, "fasta")
    temp_file.close()
    alignment_filename = f"{temp_file.name}-muscle.fa"
    # raise an error if the alignment file already exists. (it won't but just in case)
    if os.path.exists(alignment_filename):
        raise FileExistsError(f"{alignment_filename} already exists")

    muscle_command = (
        f'{muscle_binary} -super5 "{temp_file.name}" -output "{alignment_filename}"'
    )
    subprocess.run(muscle_command, shell=True, check=True)

    if output_type == "list":
        muscle_output = tools.import_fasta(alignment_filename, output_format="list")
    elif output_type == "dict":
        muscle_output = tools.import_fasta(alignment_filename, output_format="dict")
    elif output_type == "alignment":
        muscle_output = AlignIO.read(alignment_filename, "fasta")

    # delete temporary file
    os.remove(alignment_filename)
    os.remove(temp_file.name)
    return muscle_output
