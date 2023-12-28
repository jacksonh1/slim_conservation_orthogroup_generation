# construct groups for all genes in an organism

This is the way that I've used this pipeline mostly. I generated the groups of filtered, clustered, least divergent orthologs (LDOs) for all of the human genes and then accessed the output files for further conservation analysis.

If you wanted to do this for another organism, you would just need to change the `--species` flag to the appropriate species code (e.g. `--species=9606_0` for human).

You can find the species code for any organism in orthoDB in the `odb11v0_species.tab` file that you downloaded from the orthoDB website. The file in this repository is actually the full sized file and so you could just use that one (`../../data/orthoDB_sample_data/odb11v0_species.tab`).


Or you could do something like this:
```python
import local_orthoDB_group_tools.database as database

odbdb = database.orthoDB_database()
for k,v in odbdb.data_species_dict.items():
    if 'sapiens' in v:
        print(k,v)
```
