import pandas as pd

import local_orthoDB_group_tools.sql_queries as sql_queries

query_ogid_list = sql_queries.odb_gene_id_2_ogid_list(odb_gene_id)


def get_og_info_from_og_list(ogid_list: list[str], ) -> pd.DataFrame:
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