import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Align, AlignIO, Seq, SeqIO

from local_env_variables import env_variables as env


def uniprotid_2_odb_gene_id_refs(uniprotid, db_path: str|Path = env.orthoDB_files.gene_refs_sqlite) -> list[str]:
    """return the odb_gene_id from a uniprot ID"""
    connection = sqlite3.connect(db_path)
    cursor = connection.cursor()
    res = cursor.execute(f"SELECT odb_gene_id FROM gene_refs WHERE Uniprotid='{uniprotid}'")
    odb_gene_ids = res.fetchall()
    odb_gene_ids = [x[0] for x in odb_gene_ids]
    connection.close()
    return odb_gene_ids


def uniprotid_2_odb_gene_id_xrefs(uniprotid, db_path: str|Path = env.orthoDB_files.gene_xrefs_sqlite) -> list[str]:
    """return the odb_gene_id from a uniprot ID"""
    connection = sqlite3.connect(db_path)
    cursor = connection.cursor()
    # res = cursor.execute(f"SELECT odb_gene_id FROM gene_xrefs WHERE xref_id='{uniprotid}' AND DB_name='UniProt'")
    # This is much faster if you don't specify the DB_name and then filter the results
    res = cursor.execute(f"SELECT * FROM gene_xrefs WHERE xref_id='{uniprotid}'")
    odb_gene_ids = res.fetchall()
    odb_gene_ids = [x[1] for x in odb_gene_ids if x[3]=='UniProt']
    connection.close()
    return odb_gene_ids


def odb_gene_id_2_species_id(odb_gene_id, db_path: str|Path = env.orthoDB_files.gene_refs_sqlite) -> str:
    """return the species ID from an orthodb ID"""
    connection = sqlite3.connect(db_path)
    cursor = connection.cursor()
    res = cursor.execute(f"SELECT species_id FROM gene_refs WHERE odb_gene_id='{odb_gene_id}'")
    species_id = res.fetchall()[0][0]
    connection.close()
    return species_id



def ogid_2_odb_gene_id_list(ogid, db_path: str|Path = env.orthoDB_files.OG2genes_sqlite) -> list[str]:
    connection = sqlite3.connect(db_path)
    cursor = connection.cursor()
    res = cursor.execute(f"SELECT odb_gene_id FROM OG2genes WHERE OG_id='{ogid}'")
    odb_gene_ids = res.fetchall()
    odb_gene_ids = [x[0] for x in odb_gene_ids]
    connection.close()
    return odb_gene_ids


def odb_gene_id_2_ogid_list(odb_gene_id, db_path: str|Path = env.orthoDB_files.OG2genes_sqlite) -> list[str]:
    connection = sqlite3.connect(db_path)
    cursor = connection.cursor()
    res = cursor.execute(f"SELECT OG_id FROM OG2genes WHERE odb_gene_id='{odb_gene_id}'")
    og_ids = res.fetchall()
    og_ids = [og_id[0] for og_id in og_ids]
    og_ids = list(set(og_ids))
    connection.close()
    return og_ids


def ogid_list_2_og_df(ogid_list: list[str], db_path: str|Path = env.orthoDB_files.ogs_sqlite) -> pd.DataFrame:
    db_path = env.orthoDB_files.ogs_sqlite
    connection = sqlite3.connect(db_path)
    cursor = connection.cursor()
    # og_df = pd.DataFrame.f
    og_query_results = []
    for og_id in ogid_list:
        # res = cursor.execute(f"SELECT species_id FROM genes WHERE odb_gene_id='{seqrecord.id.split('|')[1]}'")
        res = cursor.execute(f"SELECT * FROM OGs WHERE OG_id='{og_id}'")
        og_query_results.append(res.fetchall()[0])

    connection.close()
    og_df=pd.DataFrame.from_records(og_query_results, columns=['index', 'OG id', "level NCBI tax id", "OG name"])
    og_df = og_df.drop('index', axis=1)
    og_df['level NCBI tax id'] = og_df['level NCBI tax id'].astype(int)
    return og_df


def get_all_odbids_from_species_id(species_id: str, db_path: str|Path = env.orthoDB_files.gene_refs_sqlite) -> list[str]:
    connection = sqlite3.connect(db_path)
    cursor = connection.cursor()
    res = cursor.execute(f"SELECT * FROM gene_refs WHERE species_id='{species_id}'")
    results = res.fetchall()
    gene_list = list(set([i[1] for i in results]))
    connection.close()
    return gene_list