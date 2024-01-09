import multiprocessing
import shutil
from pathlib import Path

import local_config.conf as conf
import local_orthoDB_group_tools.sql_queries as sql_queries
import local_scripts.odb_group_pipeline as pipeline
import local_scripts.create_filemap as create_filemap

SPECIES_ID = "9606_0"
N_CORES = multiprocessing.cpu_count()-2
OVERWRITE = True


def multiple_levels(
    config: conf.PipelineParams, query_odb_gene_id: str, og_levels: list
):
    '''
    run the pipeline for a single odb_gene_id for multiple og_levels
    '''
    for og_level in og_levels:
        config.og_select_params.OG_level_name = og_level
        try:
            pipeline.main_pipeline(config, odb_gene_id=query_odb_gene_id)
        except ValueError as err:
            # logger.error(f"{query_geneid} - {og_level} - {err}")
            print(f"{query_odb_gene_id} - {og_level} - {err}")


def main(
    config: conf.PipelineParams,
    og_levels: list,
    multiprocess=True,
    species_id=SPECIES_ID,
    n_cores=N_CORES,
    overwrite=False,
):
    odbgeneid_list = sql_queries.get_all_odb_gene_ids_from_species_id(species_id)
    if Path(config.main_output_folder).exists():
        if overwrite:
            shutil.rmtree(config.main_output_folder)
        else:
            raise FileExistsError(
                f"main_output_folder already exists: {config.main_output_folder}. Use -o flag to overwrite"
            )
    if multiprocess:
        p = multiprocessing.Pool(n_cores)
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
    main(
        config,
        og_levels,
        multiprocess=True,
        species_id=SPECIES_ID,
        n_cores=N_CORES,
        overwrite=OVERWRITE,
    )
    create_filemap.create_filemap(config.main_output_folder, output_file=Path(config.main_output_folder) / "filemap.json")
