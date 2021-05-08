# LSEA v0.2.1

*This is a new version of the tool that is now being extensively benchmarked. Please note that this is a development repository. Please use the contents of this repository with caution*

## Prerequisites
- plink v1.9
- bedtools
- scipy for Python3
- rpy2 for Python3
- LDSC (https://github.com/bulik/ldsc)

## Usage

To see description for every option run `python3 LSEA_2.0.py --help`.

## Data format

Input file (tsv) should contain headers. You should pass as parameter 4 column names - for chromosome, position, id and pval. See  `data/in.tsv`.
Output file contains pathway, p-value and q-value. See `data/result.tsv` for example.

## Examples
To run full pipeline with *plink* run:  
`
python3 LSEA_2.1.py -input ./data/in.tsv -plink_dir <plink_folder_path> -bfile <.bed file for PLINK> -universe ./data/universe.json -n 100000 -m 100 -h2 0.22
`


To use pregenerated file from *plink* run:  
`
python3 LSEA_2.0.py -tsv ./in.tsv -use_clumped ./<clumped_file>.clumped -interval 100000 -gene_file ./data/gencode_formatted.tsv -json ./data/universe.json -column_names CHR BP SNP P
`  

To generate universe run:  
`
python3 universe_generator.py -path ./in.tsv -gene_file ./data/gencode_formatted.tsv -msig_path ./data/c2.all.v7.0.symbols.gmt  -interval 100000
`
