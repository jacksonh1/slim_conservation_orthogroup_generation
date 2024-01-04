import local_config.conf as conf
import local_orthoDB_group_tools.sql_queries as sql_queries
import local_scripts.odb_group_pipeline as pipeline
import pandas as pd


def main(config, uniprot_id_list):
    for uniprot_id in uniprot_id_list:
        try:
            pipeline.main_pipeline(config, uniprot_id=uniprot_id)
        except ValueError as err:
            print(f"{uniprot_id} - {err}")

if __name__ == "__main__":
    config = pipeline.load_config("./params.yml")
    df = pd.read_csv('./table.csv')
    uniprot_id_list = df['Uniprotid'].unique().tolist()
    main(config, uniprot_id_list)
