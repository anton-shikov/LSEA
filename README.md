## Prerequisites
- plink v1.9
- bedtools
- scipy for Python3
- rpy2 for Python3

## Usage

To see description for every option run `python3 LSEA_2.0.py --help`.

## Data format

Input file (tsv) should contain headers. You should pass as parameter 4 column names - for chromosome, position, id and pval. See  `data/in.tsv`.
Output file contains pathway, p-value and q-value. See `data/result.tsv` for example.

## Examples
To run full pipeline with *plink* run:  
`
python3 LSEA_2.0.py -tsv ./in.tsv -pldir <plink_folder_path> -b_file <files_for_plink> -gene_file ./data/gencode_formatted.tsv -json ./data/universe.json -column_names CHR BP SNP P
`


To use pregenerated file from *plink* run:  
`
python3 LSEA_2.0.py -tsv ./in.tsv -use_clumped ./<clumped_file>.clumped -interval 100000 -gene_file ./data/gencode_formatted.tsv -json ./data/universe.json -column_names CHR BP SNP P
`  

To generate universe run:  
`
python3 universe_generator.py -path ./in.tsv -gene_file ./data/gencode_formatted.tsv -msig_path ./data/c2.all.v7.0.symbols.gmt  -interval 100000
`
## Methods
Work of the program can be described as folows:
1. Create "clumps" from input SNPs using *plink*
2. Create intervals around their centers
3. Merge intersecting intervals
4. For each gene find, how many intervals it is in
5. For every pathway count, how many intervals correspond to it
6. Use FDR correction.

## Citations

- Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ & Sham PC (2007) PLINK: a toolset for whole-genome association and population-based linkage analysis. American Journal of Human Genetics, 81.
- Aaron R. Quinlan, Ira M. Hall, BEDTools: a flexible suite of utilities for comparing genomic features, Bioinformatics, Volume 26, Issue 6, 15 March 2010, Pages 841–842
