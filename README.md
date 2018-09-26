# Enterobase Typer

This is a script to type assemblies according to schemes retrieved from Enterobase.

Enterobase databases can be manually retrieved with [EnterobasePull](https://github.com/bfssi-forest-dussault/EnterobasePull).
Database files retrieved from Enterobase can be automatically formatted via makeblastdb by the typer script.

### Usage
```bash
Usage: enterobase_typer.py [OPTIONS]

Options:
  -i, --input_assembly PATH  Path to input assembly in FASTA format
                             [required]
  -db, --database PATH       Path to your MLST database  [required]
  -o, --outdir PATH          Root directory to store all output files
                             [required]
  --create_db                Set this flag to re-create the blastDB files
                             using makeblastdb in the specified database
                             directory.
  --version                  Specify this flag to print the version and exit.
  --help                     Show this message and exit.
  ```