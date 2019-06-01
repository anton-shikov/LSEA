# LSEA
LSEA (locus set enrichment analysis) is a tool for performing gene set enrichment analysis on independant loci (clumps) regarding LD (linkage disequilibrium).

## Getting Started

LSEA could be applied for gene set enrichment analysis for data obtained from GWAS-summary statistics files in tsv-format. It is based on simple hypergeometric test, however it transforms genes and gene sets into clumps and clum sets to eliminate multiple signals from genes in LD to enchanse analysis precision.

### Prerequisites

-python 3.7 or higher
-scipy
```
~$ pip3 install scipy
```
-pandas
```
~$ pip3 install pandas
```
-numpy
```
~$ pip3 install numpy
```
-plink (link)
-SNPeff (link)

### Installing

To install this tool clone this repository to your PC.

```
~$ git clone https://github.com/LSEA
```

## Running and using tool

Firstly, you need to prepare tsv-file with GWAS summary statistics with the following structure: <br> 
CHR	COORDINATE	RSID	REF	ALT	PVAL

To launch this tool you also need to specify path to plink and snpeff directories.

## Example of usage
```
~$ python3 LSEA.py -af <input tsv-file> -sn <path to SNPeff> -pld <path to plink> -bf <bfile for plink> -p
```
This command will apply LSEA algorytm to input file and will generate tsv file with the following structure: <br>
BIOCARTA_INTRINSIC_PATHWAY	2.0446237642438122e-14	2.2517441515617103e-10	(17776, 11, 36, 6, 'F11;FGB;FGA;F5;FGG;KLKB1') <br>
The first column contains the name of the set, the second and the third represent p-value and corrected q-value of hypergeometric test, the last coloumn includes information about total number of independant loci, number of loci in quiery, number of loci in gene set, number of loci common for quiery and gene set and, finally, the genes list.<br>
Note that the genes list could be smaller then number of genes, because only indepedant loci are counted for analysis.
-p (--precompiled flag) points that precompiled universe of clumpes based of UK Biobank data is used.

##Information about flags:  
  -af File              Input file in tsv-format
  -vf File              Annotated vcf-file if it is already prepared
  -pl File              Plink result of clumping if it is already prepared
  --precompiled, -p     use precompiled loci
  -sn File              Path to SNPeff
  -g Str                Flag for specifying genome for SNPeff annotattion
  -pld Dir              Path to plink
  -bf Str               Bfile for plink

##Creating your own universe:  
If you don't want to use precompiled clumped universe you can use options for creating your own clumping universe based on GWAS summary statictics. Use -cu (--create_universe) option to create clump universe for yor data. For this function you have to prepare results of clumping for GWAS data that you need to get results of clumping (.clumped files). If you have multivle files (e.g. for different phenotypes use -cld <directory> flag for specifying directory with clumped files)


## Author

* **Anton Shikov** - *Initial work* - [anton-shikov](https://github.com/anton-shikov)


## License

This project is free and available for everyone.

