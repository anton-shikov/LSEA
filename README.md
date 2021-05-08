# LSEA v0.2.1

*This is a new version of the tool that is now being extensively benchmarked. Please note that this is a development repository. Please use the contents of this repository with caution*

## Prerequisites
- plink v1.9
- bedtools
- scipy for Python3
- rpy2 for Python3
- LDSC (https://github.com/bulik/ldsc)

## Usage

To run full pipeline:  
```
python3 LSEA_2.1.py -input ./data/in.tsv \
                    -plink_dir <plink_folder_path> \
                    -bfile <.bed file for PLINK> \
                    -universe ./data/universe.json \
                    -n 100000 -m 100 -h2 0.22
```

To generate universe from BED and GMT files, run:  
```
python3 universe_generator.py -variants ./data/positions.tsv \
                              -features ./data/gencode_formatted.bed \
                                        ./data/c2.all.v7.0.symbols.gmt \
                              -interval 100000 -o universe.json
```

To generate universe using a directory with per-set BED files, run:
```
python3 universe_generator.py -variants ./data/positions.tsv \
                              -feature_files_dir <dir with BED files> \
                              -interval 100000 -o universe_dir.json
