# clusterDBLalpha
A python script to generate an otu table from cleaned DBLa sequences.

##Installation
The script depends on python2.7. The location of the usearch binary needs to be hardcoded into the top of the python script.

##Usage
The fasta file is assumed to be either in the format used by Thomas Rask's pipeline or the updated pipeline which uses the Usearch format to assign reads to isolates.
```
usage: clusterDBLa.py [-h] -o OUTPUTDIR -r READ [--perID PERID] [--cpu CPU]
                      [--verbose]

Cluster cleaned DBLalpha sequence tags.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUTDIR, --outputDir OUTPUTDIR
                        location of output directory. Will be created if it
                        doesn't exist
  -r READ, --read READ  location of fasta file containing sequences.
  --perID PERID         percentage ID threshold. (default=0.96)
  --cpu CPU             number of cpus to use. (default=1)
  --verbose             print verbose output (default=False)
```
