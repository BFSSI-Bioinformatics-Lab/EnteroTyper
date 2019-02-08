# EnterobaseTyper

This is a set of scripts to type assemblies according to schemes retrieved from Enterobase.

Enterobase databases can be manually retrieved with [EnterobasePull](https://github.com/bfssi-forest-dussault/EnterobasePull). Database files retrieved from Enterobase can be automatically formatted via makeblastdb by the typer script.

### Single-sample Usage (enterobase_typer.py)
```
Usage: enterobase_typer.py [OPTIONS]

Options:
  -i, --input_assembly PATH  Path to input assembly in FASTA format
                             [required]
  -db, --database PATH       Path to your MLST database  [required]
  -o, --out_dir PATH         Root directory to store all output files
                             [required]
  --create_db                Set this flag to create the blastDB files using
                             makeblastdb in the specified database
                             directory.Will re-create the database files if
                             they are already present.
  -v, --verbose              Set this flag to enable more verbose logging.
  --version                  Specify this flag to print the version and exit.
  --help                     Show this message and exit.
  ```

### Multi-sample Usage (enterobase_typer_multi.py)
```
Usage: enterobase_typer_multi.py [OPTIONS]

Options:
  -i, --input_dir PATH  Path to directory containing FASTA assemblies
                        [required]
  -db, --database PATH  Path to your MLST database  [required]
  -o, --out_dir PATH    Root directory to store all output files  [required]
  --create_db           Set this flag to create the blastDB files using
                        makeblastdb in the specified database directory.Will
                        re-create the database files if they are already
                        present.
  -v, --verbose         Set this flag to enable more verbose logging.
  --version             Specify this flag to print the version and exit.
  --help                Show this message and exit.
```

### Sequence Concatenation Pipeline (sequence_concatenation_pipeline.py)

This helper script will take the output reports from the
enterobase_typer.py script as input to generate a concatenated sequence
file that can be fed to tree generating software.

```
Usage: sequence_concatenation_pipeline.py [OPTIONS] [TARGETS]...

  Takes a list of target *.BLASTn_Detailed_Report.tsv files, followed by
  several options. Extracts sequences from the BLASTn report files as FASTA
  files, aligns them all with MUSCLE, and then concatenates all sequences
  into a single FASTA.

Options:
  -o, --outdir PATH     Root directory to store all output files  [required]
  -db, --database PATH  Path to your MLST database  [required]
  -v, --verbose         Set this flag to enable more verbose logging.
  --help                Show this message and exit.
```
