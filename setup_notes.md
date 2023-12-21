Installed my own module (`./src/`) using setuptools. <br>
I took inspiration from this example: https://github.com/ericmjl/Network-Analysis-Made-Simple <br>
combined with information from here: https://setuptools.pypa.io/en/latest/userguide/package_discovery.html

first I activated a virtual environment to install my module into: <br> 
`conda activate `

then ran the following command to install: <br>
`pip install -e .` <br>
output: <br>
```
Obtaining file:///Users/jackson/Dropbox%20%28MIT%29/work/07-SLiM_bioinformatics/01
  Preparing metadata (setup.py) ... done
Installing collected packages: local-tools
  Running setup.py develop for local-tools
Successfully installed local-tools-0.1
```


This is an editable install. So as files are added to `./src/local_seqtools/` or the files are changed, the changes will be updated when you use the tools
