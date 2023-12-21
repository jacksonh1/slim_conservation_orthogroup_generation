echo "building SQLite databases from orthoDB tables"
echo "This will probably take a while"
python ./scripts-gen_SQLite_dbs/make_SQLite_database_fasta.py
python ./scripts-gen_SQLite_dbs/make_SQLite_database_genes.py
python ./scripts-gen_SQLite_dbs/make_SQLite_database_gene_xrefs.py
python ./scripts-gen_SQLite_dbs/make_SQLite_database_OGs.py
python ./scripts-gen_SQLite_dbs/make_SQLite_database_OG2genes.py
