import copy
import json
import os
import pathlib
import re
from pathlib import Path

import local_orthoDB_tools_v3.sql_queries as sql_queries
import pandas as pd
from Bio import Seq, SeqIO
from pyprojroot import here

import local_env_variables.env_variables as env
import local_orthoDB_group_tools.uniprotid_search as uniprotid_search


class orthoDB_database:

    database_files = env.orthoDB_files

    def __init__(self):
        self._load_data_all_seqs()
        self._load_data_species_dict()
        self._load_data_levels_df()

    # ==============================================================================
    # // DATA LOADING METHODS
    # ==============================================================================
    def _load_data_all_seqs(self):
        self.data_all_seqrecords_dict = SeqIO.index_db(
            str(self.database_files.all_seqs_sqlite),
            str(self.database_files.all_seqs_fasta),
            "fasta",
        )

    def _load_data_species_dict(self):
        self.data_species_dict = (
            pd.read_csv(
                self.database_files.species_tsv,
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
            )[["species ID", "species name"]]
            .set_index("species ID")
            .to_dict()["species name"]
        )

    def _load_data_levels_df(self):
        self.data_levels_df = pd.read_csv(
            self.database_files.levels_tsv,
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


class orthoDB_query:
    """
    This object's role should be to store the information about the selected ortholog group and the query protein, and to provide methods to query the orthoDB database for information about the ortholog group.
    """

    odb_database = orthoDB_database()

    def __init__(self, output_dir_base=None):
        if output_dir_base is None:
            output_dir_base = "./orthoDB_analysis"
        self.output_dir_base = Path(output_dir_base)
        self.output_dir_base.mkdir(parents=True, exist_ok=True)
        self.output_file_dict = {}
        self.filter_dict = {}

    # ==============================================================================
    # // helper functions for finding the protein in the database using uniprot id
    # ==============================================================================
    def _uniprotid_2_geneid_human(self, uniprotid, duplicate_action="longest"):
        """
        Given a uniprot id, return the gene id (orthodb id)
        If multiple gene ids are found, return the `duplicate_action` one (first or longest)
        """
        # TODO: switch lookup to SQL query
        if duplicate_action not in ["first", "longest"]:
            raise ValueError(
                f"duplicate_action must be 'first' or 'longest', not {duplicate_action}"
            )
        matches = self._search_df(
            self.odb_database.data_human_gene_key_df, "UniprotID", uniprotid
        )
        if matches is None:
            print(f"{uniprotid} not found in gene key table, searching in xref table")
            matches = self._search_df(
                self.odb_database.data_human_xref_gene_key_df, "UniprotID", uniprotid
            )
        if matches is None:
            print("not found in xref key table")
            raise ValueError(
                f"Uniprot id `{uniprotid}` not found in human gene key or xref tables"
            )
        self.query_uniprot_id = uniprotid
        self._matches = matches
        if len(matches) > 1:
            print(f"Multiple matches for `{uniprotid}` in table")
            print(matches)
            print('view matches in "self._matches"')
            print(f'choosing a single id from multiple matches by "{duplicate_action}"')
            if duplicate_action == "first":
                query_gene_id = matches["gene ID"].values[0]
                self.query_gene_id = query_gene_id
                return query_gene_id
            if duplicate_action == "longest":
                seq_list = []
                for gene_id in matches["gene ID"].values:
                    seq = self.odb_database.data_all_seqrecords_dict[gene_id]
                    seq_list.append(seq)
                sorted_seq_list = sorted(
                    seq_list, key=lambda x: len(x.seq), reverse=True
                )
                self.duplicate_uniprot_match_seq_str_dict = {i.id: str(i.seq) for i in sorted_seq_list}
                query_gene_id = sorted_seq_list[0].id
                self.query_gene_id = query_gene_id
                return query_gene_id
        else:
            query_gene_id = matches["gene ID"].values[0]
            self.query_gene_id = query_gene_id
            return query_gene_id

    def _get_geneid_info(self, query_gene_id=None):
        if query_gene_id is None:
            query_gene_id = self.query_gene_id
        gene_info = self.odb_database.data_human_gene_key_df[
            self.odb_database.data_human_gene_key_df["gene ID"] == query_gene_id
        ]

        # storing info in class variables
        self.query_description = gene_info.description.values[0]
        self.query_species_id = gene_info["species ID"].values[0]
        self.query_species_name = self.odb_database.data_species_dict[
            self.query_species_id
        ]
        self.query_ogid_list = self.odb_database.data_geneid_2_og_list_dict[
            query_gene_id
        ]
        # if there are no OGs for this gene, raise a ValueError
        if len(self.query_ogid_list) == 0:
            raise ValueError(f"no ortholog groups found for gene id {query_gene_id}")
        qseq = self.odb_database.data_all_seqrecords_dict[query_gene_id]
        self.query_sequence = qseq
        # self.query_sequence_id = str(qseq.id)

    def _available_OGs_info(self):
        """get info about the list of OGs available for the selected geneid"""
        query_og_df = sql_queries.ogid_list_2_og_df(self.query_ogid_list)
        query_level_df = self.odb_database.data_levels_df[
            self.odb_database.data_levels_df["level NCBI tax id"].isin(
                query_og_df["level NCBI tax id"].unique()
            )
        ].copy()
        query_available_OGs_info_df = pd.merge(
            query_og_df, query_level_df, on="level NCBI tax id", how="inner"
        )
        query_available_OGs_info_df = query_available_OGs_info_df[
            [
                "OG id",
                "level NCBI tax id",
                "level name",
                "total non-redundant count of species underneath",
                "OG name",
            ]
        ]
        query_available_OGs_info_df[
            "total non-redundant count of species underneath"
        ] = query_available_OGs_info_df[
            "total non-redundant count of species underneath"
        ].astype(
            float
        )
        # query_OG_info_df = query_OG_info_df.infer_objects()
        self.query_available_OGs_info_df = query_available_OGs_info_df.sort_values(
            by="total non-redundant count of species underneath"
        )

    # ==============================================================================
    # // DRIVER - finding and setting the geneid in the database
    # ==============================================================================
    def driver_set_geneid_by_uniprotid_search(
        self, uniprotid, duplicate_action="longest"
    ):
        self.query_gene_id = uniprotid_search.uniprotid_2_geneid(uniprotid, duplicate_action=duplicate_action)
        self.query_uniprot_id = uniprotid
        _ = self._get_geneid_info()
        self._available_OGs_info()

    def driver_set_geneid_directly(self, geneid):
        """
        TODO: get UniprotID from gene or xref table. Should probably use SQL query for this
        """
        self.query_gene_id = geneid
        _ = self._get_geneid_info()
        self._available_OGs_info()

    # ==============================================================================
    # // writing output to files
    # ==============================================================================
    def set_query_output_directory(self, output_folder_name=None):
        if output_folder_name is None:
            if hasattr(self, "query_uniprot_id"):
                output_dir_path = (
                    self.output_dir_base
                    / f"{self.query_uniprot_id}-{self.query_gene_id}"
                    / f"{self.selected_query_og_level_name}"
                )
            else:
                output_dir_path = (
                    self.output_dir_base
                    / f"{self.query_gene_id}"
                    / f"{self.selected_query_og_level_name}"
                )
        else:
            output_dir_path = self.output_dir_base / output_folder_name
        self.query_output_dir_path = output_dir_path
        self.query_output_dir_path.mkdir(parents=True, exist_ok=True)

    def write_all_sequences(self, file_prefix=None):
        """
        TODO: put the filenames into in the set_query_output_directory function
        """
        if file_prefix is None:
            file_prefix = self.file_prefix
        output_filename1 = self.query_output_dir_path / f"{file_prefix}_full_OG.fasta"
        self._write_sequence_list_2_fasta(
            list(self.sequences_full_OG_dict.values()), output_filename1
        )
        self.output_file_dict["fasta sequences - full OG"] = output_filename1

        output_filename2 = self.query_output_dir_path / f"{file_prefix}_OG_LDOs.fasta"
        self._write_sequence_list_2_fasta(self.sequences_LDO_list, output_filename2)
        self.output_file_dict["fasta sequences - OG LDOs"] = output_filename2

        output_filename3 = (
            self.query_output_dir_path / f"{file_prefix}_OG_LDOs_cdhit.fasta"
        )
        self._write_sequence_list_2_fasta(
            self.sequences_OG_LDO_cdhit_list, output_filename3
        )
        self.output_file_dict["fasta sequences - OG LDOs cdhit"] = output_filename3

        if hasattr(self, "alignment_OG_LDO_cdhit"):
            output_filename4 = (
                self.query_output_dir_path
                / f"{file_prefix}_OG_LDOs_cdhit_mafftaln.fasta"
            )
            self._write_sequence_list_2_fasta(
                self.alignment_OG_LDO_cdhit, output_filename4
            )
            self.output_file_dict["fasta alignment - OG LDO cdhit"] = output_filename4

    def write_out_params(self):
        if hasattr(self, "query_uniprot_id"):
            info_json_filename = (
                self.query_output_dir_path
                / f"{self.query_uniprot_id}_{self.selected_query_og_level_name}_info__dict__.json"
            )
        else:
            info_json_filename = (
                self.query_output_dir_path
                / f"{self.query_gene_id}_{self.selected_query_og_level_name}_info__dict__.json"
            )
        self.output_file_dict["json - query and ortho group info"] = info_json_filename
        self.output_file_dict_absolute = {
            k: str(v.absolute()) for k, v in self.output_file_dict.items()
        }
        self.output_file_dict = {k: str(v) for k, v in self.output_file_dict.items()}
        self.output_file_dict_absolute = {
            k: str(v) for k, v in self.output_file_dict_absolute.items()
        }
        output_dict = {}
        output_dict["json_info_file"] = str(info_json_filename)
        for k, v in self.__dict__.items():
            if k.startswith("sequence") or k.startswith("alignment"):
                continue
            if k == "query_sequence":
                continue
            if type(v) in [pd.DataFrame, pd.Series]:
                continue
            if type(v) == pathlib.PosixPath:
                output_dict[k] = str(v)
                continue
            output_dict[k] = v

        output_dict["num_sequences_full_OG"] = len(self.sequences_full_OG_dict)
        output_dict["num_sequences_LDO"] = len(self.sequences_LDO_list)
        output_dict["num_sequences_LDO_cdhit"] = len(self.sequences_OG_LDO_cdhit_list)
        self.output_dict = output_dict
        with open(info_json_filename, "w") as f:
            json.dump(output_dict, f, indent=4)
        # with open(self.query_output_dir_path / 'filter_results.json', 'w') as f:
        # json.dump(self.filter_dict, f, indent=4)

    # ==============================================================================
    # // DRIVER - write files
    # ==============================================================================
    def write_files(self, output_folder_name=None):
        self.set_query_output_directory(output_folder_name=output_folder_name)
        self.file_prefix = self.query_output_dir_path.stem
        self.write_all_sequences()
        self.write_out_params()

    # ==============================================================================
    # // functions to print info to screen
    # ==============================================================================
    def print_database_filepaths(self):
        for k, v in self.odb_database.database_files.items():
            print(f"- {k}\n    {v}")

    def print_orthoDB_readme(self):
        with open(self.odb_database.database_files["readme"], "r") as f:
            print(f.read())

    # ==============================================================================
    # // general utilities (static methods)
    # ==============================================================================
    @staticmethod
    def _search_df(df, column, query):
        matches = df[df[column] == query].copy()
        if len(matches) == 0:
            print(f"`{query}` not found in {column} column")
            return None
        return matches

    @staticmethod
    def _convert_sequence_dict_to_list(seq_dict):
        return [seq for seq in seq_dict.values()]

    @staticmethod
    def _write_sequence_list_2_fasta(seqrecord_list, filename):
        with open(filename, "w") as f:
            SeqIO.write(seqrecord_list, f, "fasta")

    # for i in ldo_seqrecord_list:
    #     score = ldo_df.loc[ldo_df["id"] == i.id, "PID"].values[0]
    #     i.id = f"{i.id}-{score:.2f}"
