# source code structure

## `local_config/`
- contains code for importing and parsing the pipeline parameters from a yml file

## `local_env_variables/`
`env_variables.py`:
- reads the environment variables from the `.env` file
- has some data importing functions and stores the main data files as objects
- This is imported in many of the other src files

## `local_orthoDB_group_tools/`
- contains the main code for each individual step of the pipeline

## `local_seqtools/`
- contains some basic tools for working with sequences

## `scripts/`
- contains the main script for running the pipeline (`odb_group_pipeline.py`)



