import multiprocessing

import local_config.conf as conf
import local_orthoDB_group_tools.sql_queries as sql_queries
import local_scripts.odb_group_pipeline as pipeline

SPECIES_ID = "9606_0"
N_CORES = multiprocessing.cpu_count()-2

def multiple_levels(
    config: conf.PipelineParams, query_odb_gene_ids: str, og_levels: list
):
    for og_level in og_levels:
        config.og_select_params.OG_level_name = og_level
        try:
            pipeline.main_pipeline(config, odb_gene_id=query_odb_gene_ids)
        except ValueError as err:
            # logger.error(f"{query_geneid} - {og_level} - {err}")
            print(f"{query_odb_gene_ids} - {og_level} - {err}")


def main(config: conf.PipelineParams, og_levels: list, multiprocess=True):
    odbgeneid_list = sql_queries.get_all_odb_gene_ids_from_species_id("9606_0")
    if multiprocess:
        p = multiprocessing.Pool(N_CORES)
        f_args = [(config, i, og_levels) for i in odbgeneid_list]
        p.starmap(multiple_levels, f_args)
        p.close()
        p.join()
    else:
        for i in odbgeneid_list:
            multiple_levels(config, i, og_levels)


if __name__ == "__main__":
    og_levels = ["Eukaryota", "Mammalia", "Metazoa", "Tetrapoda", "Vertebrata"]
    config = pipeline.load_config("./params.yml")
    main(config, og_levels, multiprocess=True)
