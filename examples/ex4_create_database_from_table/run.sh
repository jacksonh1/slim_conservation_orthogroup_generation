python "../../src/local_scripts/map_uniprotid.py" -i "./table.csv" --uni_column "Uniprotid"
python "../../src/local_scripts/pipeline_input_table.py" -c "./params.yaml" -o --table "./table_mapped_odbgeneid.csv" --odb_gene_id_column gene_id
python "../../src/local_scripts/create_filemap.py" --main_output_folder "./output/" --output_file "./output/database_key.json"