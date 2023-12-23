import json
import os

from pyprojroot import here

from local_orthoDB_group_tools import (cluster, database, filters, find_LDOs,
                                       og_select)

root = here()
default_params_json = root / "src" / "local_orthoDB_tools" / "params.json"
DEFAULT_PARAM_DICT = {
    "main_output_folder": "./orthoDB_analysis",
    "OG selection method": "level name",
    "min_fraction_short_than_query": 0.5,
    "OG level name": "Eukaryota",
    "LDO selection method": "msa_by_organism",
    "align": True,
    "n_align_threads": 32,
    "write files": True
}

def str2bool(s: str):
    if s == 'True':
        return True
    elif s == 'False':
        return False
    else:
        raise ValueError("str2bool only accepts 'True' or 'False'")


def load_params(user_parameters, default_parameters=DEFAULT_PARAM_DICT):
    if isinstance(user_parameters, str):
        with open(user_parameters, "r") as f:
            user_parameters = json.load(f)
    
    for key in user_parameters.keys():
        if key not in default_parameters.keys():
            raise ValueError(f"parameter {key} not recognized")
    for k,v in default_parameters.items():
        if k not in user_parameters.keys():
            print(f"`{k}` parameter not provided")
            print(f"\t- using default value: `{v}`\n")

    # update default params with params provided by the user
    params = default_parameters.copy()
    params.update(user_parameters)
    # if params['align'] is a boolean, use it as is
    # otherwise, try to convert it to a bool
    if not isinstance(params["align"], bool):
        params["align"] = str2bool(params["align"])
    if not isinstance(params["write files"], bool):
        params["write files"] = str2bool(params["write files"])
    params["n_align_threads"] = int(params["n_align_threads"])
    params["min_fraction_short_than_query"] = float(params["min_fraction_short_than_query"])
    return params


def pipeline(query_uniprotid, user_params=None, linux=True):
    """run the orthoDB pipeline
    steps (arguments provided in params.json file):
    1. Find uniprotid in orthoDB database
    2. Select orthologous group (OG) based on criteria
    3. Filter sequences based on criteria
    4. select least divergent orthologs (LDOs) based on criteria
    5. cluster LDOs using cdhit
    6. write output files

    example params or parameter json file:
    {
        "main_output_folder": "./orthoDB_analysis",
        "OG selection method": "level name",
        "min_fraction_short_than_query": "0.5",
        "OG level name": "Eukaryota",
        "LDO selection method": "msa_by_organism",
        "align": "True",
        "n_align_threads": "32",
        "write files": "True"
    }

    Parameters
    ----------
    query_uniprotid : str
        uniprot id of sequence you want to find orthologs for
    params : dict
        the processing parameters in a dict

    Returns
    -------
    local_orthoDB_tools.database.orthoDB_query
        local_orthoDB_tools.database.orthoDB_query object
    """
    # TODO: all this default params stuff makes me think I should maybe use argparse or something else
    # could put the params handling in a separate function and call it from here and from the command line

    if user_params is None:
        user_params = {}
    params = load_params(user_params)

    output_folder = os.path.join(params["main_output_folder"], params["LDO selection method"])

    odbquery = database.orthoDB_query(
        output_folder
    )

    # select geneid in odb database
    try: 
        odbquery.driver_set_geneid_by_uniprotid_search(query_uniprotid)
    except ValueError:
        # print(e)
        print("Uniprotid not found in orthoDB database")
        print("or protein has no ortholog groups in the database")
        raise

    # Ortholog group selection
    print(f"selecting ortholog group for {query_uniprotid}")
    assert params["OG selection method"] in [
        "most_species",
        "target_number_of_species",
        "level name",
    ], "OG selection method not recognized. must be one of: most_species, target_number_of_species, level name"
    if params["OG selection method"] == "most_species":
        og_select.select_OG_by_level_with_most_species(odbquery)
    if params["OG selection method"] == "target_number_of_species":
        og_select.select_OG_by_target_number_of_species(
            odbquery, int(params["target number of species"])
        )
    if params["OG selection method"] == "level name":
        og_select.select_OG_by_level_name(odbquery, params["OG level name"])

    # filters
    print(f"filtering sequences for {query_uniprotid}")
    filters.filter_sequences_with_X(odbquery)
    filters.filter_sequences_by_length(odbquery, params["min_fraction_short_than_query"])

    print(f"finding LDOs for {query_uniprotid}")
    # Least divergent orthologs
    assert params["LDO selection method"] in [
        "msa_by_organism",
        "alfpy_google_distance",
        "pairwise",
        "msa",
    ], "LDO selection method not recognized. must be one of: msa_by_organism, alfpy_google_distance, pairwise, msa"
    find_LDOs.main(odbquery, params["LDO selection method"], n_align_threads=params['n_align_threads'])

    # clustering
    print(f"running CDHIT for {query_uniprotid}")
    cluster.driver_CDHIT(odbquery, align=params['align'], n_align_threads=params['n_align_threads'], linux=linux)
    
    ## TODO: add filepath of folder 

    if params["write files"]:
        print(f"saving files for {query_uniprotid}")
        odbquery.write_files()
        with open(os.path.join(output_folder, 'params.json'), 'w') as f:
            json.dump(params, f, indent=4)

    return odbquery


def pipeline_json_params(query_uniprotid, parameter_file_json=None, linux=True):
    """run the orthoDB pipeline
    steps (arguments provided in params.json file):
    1. Find uniprotid in orthoDB database
    2. Select orthologous group (OG) based on criteria
    3. Filter sequences containing X/* values and sequences shorter than query (by fraction of query length)
    4. select least divergent orthologs (LDOs) based on criteria
    5. cluster LDOs using cdhit
    6. write output files

    example params.json file:
    {
        "main_output_folder": "./orthoDB_analysis",
        "OG selection method": "level name",
        "min_fraction_short_than_query": "0.5",
        "OG level name": "Eukaryota",
        "LDO selection method": "msa_by_organism",
        "align": "True",
        "n_align_threads": "46",
        "write files": "True"
    }

    Parameters
    ----------
    query_uniprotid : str
        uniprot id of sequence you want to find orthologs for
    parameter_file_json : str, optional
        the processing parameters in a json file, by default "./params.json"

    Returns
    -------
    local_orthoDB_tools.database.orthoDB_query
        local_orthoDB_tools.database.orthoDB_query object
    """    
    if parameter_file_json is None:
        print("No parameter file provided, using default parameters")
        for k, v in DEFAULT_PARAM_DICT.items():
            print(f"{k}: {v}")
        param_dict = DEFAULT_PARAM_DICT
    else:
        with open(parameter_file_json) as f:
            param_dict = json.load(f)
    odbquery = pipeline(query_uniprotid, param_dict, linux=linux)
    return odbquery


    # if isinstance(user_params, str):
    #     with open(user_params, "r") as f:
    #         user_params = json.load(f)
    
    # for key in user_params.keys():
    #     if key not in DEFAULT_PARAM_DICT.keys():
    #         raise ValueError(f"parameter {key} not recognized")
    # for k,v in DEFAULT_PARAM_DICT.items():
    #     if k not in params.keys():
    #         print(f"parameter {key} not provided")
    #         print(f"using default value: {v}, for parameter {key}")

    # # update default params with params provided by the user
    # params = DEFAULT_PARAM_DICT.copy()
    # params.update(user_params)


def pipeline_setID_directly(odb_gene_id, user_params={}, linux=True):
    """run the orthoDB pipeline
    steps (arguments provided in params.json file):
    1. Find odb_gene_id in orthoDB database
    2. Select orthologous group (OG) based on criteria
    3. Filter sequences based on criteria
    4. select least divergent orthologs (LDOs) based on criteria
    5. cluster LDOs using cdhit
    6. write output files

    example params or parameter json file:
    {
        "main_output_folder": "./orthoDB_analysis",
        "OG selection method": "level name",
        "min_fraction_short_than_query": "0.5",
        "OG level name": "Eukaryota",
        "LDO selection method": "msa_by_organism",
        "align": "True",
        "n_align_threads": "32",
        "write files": "True"
    }

    Parameters
    ----------
    odb_gene_id : str
        ortho db gene id of the sequence you want to find orthologs for
    params : dict
        the processing parameters in a dict

    Returns
    -------
    local_orthoDB_tools.database.orthoDB_query
        local_orthoDB_tools.database.orthoDB_query object
    """
    # TODO: all this default params stuff makes me think I should maybe use argparse or something else
    # could put the params handling in a separate function and call it from here and from the command line

    params = load_params(user_params)

    output_folder = os.path.join(params["main_output_folder"], params["LDO selection method"])

    odbquery = database.orthoDB_query(
        output_folder
    )

    # select geneid in odb database
    try: 
        odbquery.driver_set_geneid_directly(odb_gene_id)
    except ValueError:
        # print(e)
        print("geneid not found in orthoDB database")
        print("or protein has no ortholog groups in the database")
        raise

    # Ortholog group selection
    print(f"selecting ortholog group for {odb_gene_id}")
    assert params["OG selection method"] in [
        "most_species",
        "target_number_of_species",
        "level name",
    ], "OG selection method not recognized. must be one of: most_species, target_number_of_species, level name"
    if params["OG selection method"] == "most_species":
        og_select.select_OG_by_level_with_most_species(odbquery)
    if params["OG selection method"] == "target_number_of_species":
        og_select.select_OG_by_target_number_of_species(
            odbquery, int(params["target number of species"])
        )
    if params["OG selection method"] == "level name":
        og_select.select_OG_by_level_name(odbquery, params["OG level name"])

    # filters
    print(f"filtering sequences for {odb_gene_id}")
    filters.filter_sequences_with_X(odbquery)
    filters.filter_sequences_by_length(odbquery, params["min_fraction_short_than_query"])

    print(f"finding LDOs for {odb_gene_id}")
    # Least divergent orthologs
    assert params["LDO selection method"] in [
        "msa_by_organism",
        "alfpy_google_distance",
        "pairwise",
        "msa",
    ], "LDO selection method not recognized. must be one of: msa_by_organism, alfpy_google_distance, pairwise, msa"
    find_LDOs.main(odbquery, params["LDO selection method"], n_align_threads=params['n_align_threads'])

    # clustering
    print(f"running CDHIT for {odb_gene_id}")
    cluster.driver_CDHIT(odbquery, align=params['align'], n_align_threads=params['n_align_threads'], linux=linux)
    
    ## TODO: add filepath of folder 

    if params["write files"]:
        print(f"saving files for {odb_gene_id}")
        odbquery.write_files()
        with open(os.path.join(output_folder, 'params.json'), 'w') as f:
            json.dump(params, f, indent=4)

    return odbquery
