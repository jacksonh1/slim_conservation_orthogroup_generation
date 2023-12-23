import os
from dataclasses import dataclass
from pathlib import Path

import dotenv
import pandas as pd
# from attrs import define, field
from attrs import frozen
from Bio import Seq, SeqIO

# from pyprojroot import here
dotenv.load_dotenv()
# dotenv.find_dotenv()
orthodb_dir = Path(os.environ['ORTHODB_DATA_DIR'])

# ==============================================================================
# // getting odb filepaths
# ==============================================================================

@frozen
class orthoDB_files_object:
    all_seqs_fasta: str = str(orthodb_dir / "odb11v0_all_og.fasta")
    all_seqs_sqlite: str = str(orthodb_dir / "odb11v0_all_og.sqlite")
    gene_refs_tsv: str = str(orthodb_dir / "odb11v0_genes.tab")
    gene_refs_sqlite: str = str(orthodb_dir / "odb11v0_genes.sqlite")
    gene_xrefs_tsv = str(orthodb_dir / "odb11v0_gene_xrefs.tab")
    gene_xrefs_sqlite = str(orthodb_dir / "odb11v0_gene_xrefs.sqlite")
    ogs_tsv: str = str(orthodb_dir / "odb11v0_OGs.tab")
    ogs_sqlite: str = str(orthodb_dir / "odb11v0_OGs.sqlite")
    OG2genes_tsv: str = str(orthodb_dir / "odb11v0_OG2genes.tab")
    OG2genes_sqlite: str = str(orthodb_dir / "odb11v0_OG2genes.sqlite")
    levels_tsv: str = str(orthodb_dir / "odb11v0_levels.tab")
    levels2species_tsv: str = str(orthodb_dir / "odb11v0_level2species.tab")
    species_tsv: str = str(orthodb_dir / "odb11v0_species.tab")
orthoDB_files = orthoDB_files_object()

# ==============================================================================
# // data loading functions
# ==============================================================================
def load_data_all_odb_seqs():
    data_all_seqrecords_dict = SeqIO.index_db(
        str(orthoDB_files.all_seqs_sqlite),
        str(orthoDB_files.all_seqs_fasta),
        "fasta",
    )
    return data_all_seqrecords_dict

# Below is an attempt to deal with the fact that the orthoDB file names might change with different versions
# file_wildcards = {
#     "all_seqs_fasta": "*all_og.fasta",
#     "gene_refs_tsv": "*_genes.tab",
#     "gene_xrefs_tsv": "*gene_xrefs.tab",
#     "ogs_tsv": "*OGs.tab",
#     "OG2genes_tsv": "*OG2genes.tab",
#     "levels_tsv": "*levels.tab",
#     "levels2species_tsv": "*level2species.tab",
#     "species_tsv": "*_species.tab",
# }

# odb_data_files = {}
# for k, v in file_wildcards.items():
#     try:
#         odb_data_files[k] = next(orthodb_dir.glob(v))
#     except StopIteration:
#         raise FileNotFoundError(f"Could not find file for `{k}`. \nSearched with wildcard {orthodb_dir}/{v}")
#     if len(list(orthodb_dir.glob(v))) > 1:
#         raise FileNotFoundError(f"Found multiple files for `{k}`. \nSearched with wildcard {orthodb_dir}/{v}")

# # The sqlite files are created by the scripts in scripts-gen_SQLite_dbs and I 
# # want to import and use this library to retrieve those file names before they are created
# # so I am not checking for their existence here.
# odb_data_files["all_seqs_sqlite"] = odb_data_files["all_seqs_fasta"].with_suffix(".sqlite")
# odb_data_files["gene_refs_sqlite"] = odb_data_files["gene_refs_tsv"].with_suffix(".sqlite")
# odb_data_files["gene_xrefs_sqlite"] = odb_data_files["gene_xrefs_tsv"].with_suffix(".sqlite")
# odb_data_files["ogs_sqlite"] = odb_data_files["ogs_tsv"].with_suffix(".sqlite")
# odb_data_files["OG2genes_sqlite"] = odb_data_files["OG2genes_tsv"].with_suffix(".sqlite")
# # later - add a way to output the sqlite files to a different directory specified in the .env file
# odb_data_files = {k: str(v) for k, v in odb_data_files.items()}

# # I'm not really sure if this is an appropriate use of attrs, but I really like
# # having the autocomplete and type checking when working with these files throughout
# # the codebase
# @frozen
# class orthoDB_files_object:
#     all_seqs_fasta: str
#     all_seqs_sqlite: str
#     gene_refs_tsv: str
#     gene_refs_sqlite: str
#     gene_xrefs_tsv: str
#     gene_xrefs_sqlite: str
#     ogs_tsv: str
#     ogs_sqlite: str
#     OG2genes_tsv: str
#     OG2genes_sqlite: str
#     levels_tsv: str
#     levels2species_tsv: str
#     species_tsv: str

# orthoDB_files = orthoDB_files_object(**odb_data_files)

'''

# ==============================================================================
# // general env variables and filepaths
# ==============================================================================

root = here()
data = root / "data"
orthodb_query_groups_root = data / 'orthoDB_ortholog_group_runs'
conservation_analysis_root = root / "conservation_analysis_runs"

# ==============================================================================
# // important filepaths
# ==============================================================================

orthodb11v0 = data / "orthoDB"
all_human_seq_fasta = orthodb11v0 / "odb_proteome" / "all_human_proteins_in_odb.fasta"
all_human_seq_clustered_fasta = orthodb11v0 / "odb_proteome" / "all_human_proteins_in_odb_clustered_c1.fasta"
pipeline_script_example_folder = root / "pipeline_scripts" / "orthoDB_table_analysis"


# ==============================================================================
# // functions to load data
# ==============================================================================

def load_data_levels_df():
    data_levels_df = pd.read_csv(
        DATABASE_FILES["levels"],
        sep="\t",
        header=None,
        names=[
            "level NCBI tax id",
            "level name",
            "total non-redundant count of genes in all underneath clustered species",
            "total count of OGs built on it",
            "total non-redundant count of species underneath",
        ],
    )
    return data_levels_df


def load_data_uniprot_gene_key_dfs():
    """
    Relevant files:
        gene xrefs uniprot human - tsv
        gene refs human - tsv
        gene refs full - tsv
    """
    data_human_gene_key_df = pd.read_csv(
        DATABASE_FILES["gene refs human - tsv"],
        sep="\t",
        header=None,
        names=[
            "gene ID",
            "species ID",
            "source ID",
            "synonyms",
            "UniprotID",
            "Ensemble",
            "NCBI id/name",
            "description",
        ],
    )
    data_human_xref_gene_key_df = pd.read_csv(
        DATABASE_FILES["gene xrefs uniprot human - tsv"],
        sep="\t",
        header=None,
        names=["gene ID", "UniprotID", "DB name"],
    ).drop("DB name", axis=1)
    # return {"gene_key": data_human_gene_key_df, "xref_gene_key": data_human_xref_gene_key_df}
    return data_human_gene_key_df, data_human_xref_gene_key_df


def load_data_all_seqs():
    data_all_seqrecords_dict = SeqIO.index_db(
        str(DATABASE_FILES["all sequences - sqlite"]),
        str(DATABASE_FILES["all sequences - fasta"]),
        "fasta",
    )
    return data_all_seqrecords_dict


def load_data_og_df():
    data_og_df = pd.read_csv(
        DATABASE_FILES["ortholog groups (OGs) - tsv"],
        sep="\t",
        header=None,
        names=["OG id", "level NCBI tax id", "OG name"],
    )
    return data_og_df


def load_species_df():
    odb_species_df = pd.read_csv(
        DATABASE_FILES["species"],
        sep="\t",
        header=None,
        names=[
            "NCBI id",
            "species ID",
            "species name",
            "assembly ID",
            "n clustered genes",
            "n OGs",
            "mapping type",
        ],
    )
    return odb_species_df

def load_xref_uniprot():
    data_xref_uniprot_gene_key_df = pd.read_csv(
        DATABASE_FILES["gene xrefs uniprot - tsv"],
        sep="\t",
        header=None,
        names=["gene ID", "UniprotID", "DB name"],
    ).drop("DB name", axis=1)
    return data_xref_uniprot_gene_key_df

'''
