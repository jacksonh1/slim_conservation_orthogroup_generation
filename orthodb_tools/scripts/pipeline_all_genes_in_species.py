#!/usr/bin/env python

import argparse
import multiprocessing
import shutil
import traceback
from pathlib import Path
from attrs import asdict

import orthodb_tools.config.orthodb_pipeline_parameters as conf
from orthodb_tools.config import orthodb_pipeline_parameters
import orthodb_tools.sql_queries as sql_queries

# import orthodb_tools.scripts.create_filemap as create_filemap
import orthodb_tools.orthogroup_processing.pipeline as pipeline

SPECIES_ID = "9606_0"
N_CORES = multiprocessing.cpu_count() - 2


def multiple_levels(
    config: conf.PipelineParams, query_odb_gene_id: str, og_levels: list
):
    """
    run the pipeline for a single odb_gene_id for multiple og_levels
    """
    for og_level in og_levels:
        config.og_select_params.OG_level_name = og_level
        try:
            pipeline.orthogroup_pipeline(config, odb_gene_id=query_odb_gene_id)
        except ValueError as err:
            traceback.print_exc()
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


def main_cli():
    OG_LEVELS = ["Eukaryota", "Mammalia", "Metazoa", "Tetrapoda", "Vertebrata"]
    d_params = ""
    for k, v in asdict(
        orthodb_pipeline_parameters.PipelineParams(),
        filter=lambda attr, value: not str(attr.name).startswith("_"),
    ).items():
        d_params += f"- {k}: {v}\n"

    parser = argparse.ArgumentParser(
        description=f"""run the pipeline for all genes in an organism. 
    processing parameters should be provided in a config file. (-c/--config)
    if no config file is provided, default parameters will be used
    The default parameters are:
{d_params}""",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-s",
        "--species_id",
        type=str,
        metavar="<str>",
        required=True,
        help=f"""species id to use. For example, '{SPECIES_ID}' for human""",
    )
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        metavar="<file>",
        default=None,
        help="""path to config file""",
    )
    parser.add_argument(
        "-n",
        "--n_cores",
        type=int,
        metavar="<int>",
        default=N_CORES,
        help=f"""number of cores to use. Default is the number available on your machine minus 2 ({N_CORES})""",
    )
    parser.add_argument(
        "-l",
        "--og_levels",
        nargs="*",
        metavar="<list>",
        default=OG_LEVELS,
        help=f"""list of phylogenetic levels at which to construct ortholog groups. Default is {OG_LEVELS}""",
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        action="store_true",
        help="""if flag is provided and the main_output_folder exists, it will be removed and overwritten by the new files. Otherwise, an error will be raised if the folder exists""",
    )
    args = parser.parse_args()
    config = pipeline.load_config(args.config)
    main(
        config,
        og_levels=args.og_levels,
        multiprocess=True,
        species_id=args.species_id,
        n_cores=args.n_cores,
        overwrite=args.overwrite,
    )
    # create_filemap.create_filemap(
    #     config.main_output_folder,
    #     output_file=Path(config.main_output_folder) / "filemap.json",
    # )


if __name__ == "__main__":
    main_cli()
