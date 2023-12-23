import copy
import sqlite3

from local_orthoDB_group_tools import database, sql_queries


def load_all_sequences_in_OG_ID(odbquery: database.orthoDB_query, ogid=None):
    """
    load all sequences into a dictionary using a OG ID
    by default it will use the self.selected_query_ogid as the OG ID
    """
    if ogid is None:
        ogid = odbquery.selected_query_ogid
    id_list = sql_queries.ogid_2_geneid_list(ogid)
    og_seq_dict = {}
    for gene_id in id_list:
        og_seq_dict[gene_id] = odbquery.odb_database.data_all_seqrecords_dict[gene_id]
    odbquery.sequences_full_OG_dict = og_seq_dict
    # return og_seq_dict


def get_OG_info(odbquery: database.orthoDB_query):
    """set instance variables with info about the ortholog group (OG) after selecting it"""
    query_OG_info_dict = odbquery.query_available_OGs_info_df[
        odbquery.query_available_OGs_info_df["OG id"] == odbquery.selected_query_ogid
    ].to_dict("records")[0]
    odbquery.selected_query_og_level_name = query_OG_info_dict["level name"]
    odbquery.selected_query_og_level_ncbi_tax_id = query_OG_info_dict[
        "level NCBI tax id"
    ]
    odbquery.selected_query_og_name = query_OG_info_dict["OG name"]


def add_species_2_seqid(odbquery: database.orthoDB_query):
    """
    Run on the full OG to avoid confusion
    """
    seqrecord_dict = copy.deepcopy(odbquery.sequences_full_OG_dict)
    for seqrecord in seqrecord_dict.values():
        # grab species number from orthoDB ID with split(':') (This might not be the species id for all orthoDB IDs, but it seems to be in this release of the database as least)
        # could try to construct a mapping dictionary from the orthoDB ID to the species name. Would have to load the giant Gene table though
        species_name = odbquery.odb_database.data_species_dict[
            seqrecord.id.split(":")[0]
        ].replace(" ", "_")

        # you could use a regex to grab the orthodbID to make sure that you don't mess up the renaming if you run this function multiple times
        # re.findall(r'\d+\_\d\:.+$', seqrecord.id)[0])
        seqrecord.id = f"{species_name}|{seqrecord.id}"
        # seqrecord.id = f'{species_name}'
        seqrecord.description = ""
    # don't think you need to return the dictionary since I think the id's are changed globally or "in place"
    # even though this is in a function
    # ids don't change unless I include the next 2 lines so maybe it's not "in place" this time but I dom't know why
    new_seqrecord_dict = {
        seqrecord.id: seqrecord for seqrecord in seqrecord_dict.values()
    }
    odbquery.sequences_full_OG_dict = new_seqrecord_dict
    # re-set the query sequence after changing the sequence ids
    odbquery.query_sequence = odbquery.sequences_full_OG_dict[
        f"{odbquery.query_species_name.replace(' ','_')}|{odbquery.query_gene_id}"
    ]
    odbquery.query_sequence_str = str(odbquery.query_sequence.seq)
    odbquery.query_sequence_id_str = str(odbquery.query_sequence.id)


def add_species_2_seqid_from_genes_SQL(odbquery: database.orthoDB_query):
    '''
    Run on the full OG to avoid confusion
    '''
    connection = sqlite3.connect(odbquery.odb_database.database_files.gene_refs_sqlite)
    cursor = connection.cursor()
    seqrecord_dict = copy.deepcopy(odbquery.sequences_full_OG_dict)
    for seqrecord in seqrecord_dict.values():
        try:
            species_id = sql_queries.odb_id_2_species_id(seqrecord.id)
        except IndexError:
            breakpoint()
        assert len(species_id) != 0, f"no species name found for gene id {seqrecord.id}"
        
        species_name = odbquery.odb_database.data_species_dict[species_id].replace(" ", "_")

        seqrecord.id = f'{species_name}|{seqrecord.id}'
        seqrecord.description = ''
    connection.close()
    new_seqrecord_dict = {seqrecord.id: seqrecord for seqrecord in seqrecord_dict.values()}
    odbquery.sequences_full_OG_dict = new_seqrecord_dict
    # re-set the query sequence after changing the sequence ids
    odbquery.query_sequence = odbquery.sequences_full_OG_dict[f"{odbquery.query_species_name.replace(' ','_')}|{odbquery.query_gene_id}"]
    odbquery.query_sequence_str = str(odbquery.query_sequence.seq)
    odbquery.query_sequence_id_str = str(odbquery.query_sequence.id)


def wrap_up_OG_selection(odbquery: database.orthoDB_query):
    """
    Run after selecting an OG
    """
    get_OG_info(odbquery)
    load_all_sequences_in_OG_ID(odbquery)
    add_species_2_seqid_from_genes_SQL(odbquery)


def select_OG_by_level_name(odbquery: database.orthoDB_query, level_name):
    filtered_OG_info_df = odbquery._search_df(
        odbquery.query_available_OGs_info_df, "level name", level_name
    )
    if filtered_OG_info_df is None:
        print(f"No OGs found for level name `{level_name}`")
        print("view available levels in `odbquery.query_available_OGs_info_df`")
        raise ValueError
        return None
    if len(filtered_OG_info_df) > 1:
        print(f"Multiple OGs found for level name `{level_name}`")
        print(odbquery.query_available_OGs_info_df)
        raise ValueError
        return None
    odbquery.selected_query_ogid = filtered_OG_info_df["OG id"].values[0]
    wrap_up_OG_selection(odbquery)


def select_OG_by_level_with_most_species(odbquery: database.orthoDB_query):
    temp = odbquery.query_available_OGs_info_df.sort_values(
        by="total non-redundant count of species underneath", ascending=False
    ).reset_index(drop=True)
    odbquery.selected_query_ogid = temp.loc[0, "OG id"]
    wrap_up_OG_selection(odbquery)


def select_OG_by_target_number_of_species(odbquery: database.orthoDB_query, target_number_of_species):
    temp = odbquery.query_available_OGs_info_df.copy()
    temp["diff"] = abs(
        target_number_of_species
        - temp["total non-redundant count of species underneath"]
    )
    odbquery.selected_query_ogid = temp.loc[temp["diff"].idxmin(), "OG id"]
    wrap_up_OG_selection(odbquery)
