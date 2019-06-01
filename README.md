# LSEA
LSEA (locus set enrichment analysis) is a tool for performing gene set enrichment analysis on independent loci, taking into account LD (Linkage disequilibrium).

## Getting Started

LSEA could be applied for gene set enrichment analysis for data obtained from GWAS-summary statistics files in tsv-format. It is based on simple hypergeometric test, however it transforms genes and gene sets into independant loci and sets of independant loci to eliminate multiple signals from genes in LD to enhance analysis precision.

### Prerequisites

-python (3.7 or higher) <br>
-scipy (1.0.0 or higher)
```
~$ pip3 install scipy
```
-pandas (0.17.1 or higher)
```
~$ pip3 install pandas
```
-numpy (1.14.1 or higher)
```
~$ pip3 install numpy
```
-plink (1.07 or higher) - http://zzz.bwh.harvard.edu/plink/ <br>
-SNPeff (4.3T or higher) - http://snpeff.sourceforge.net/ <br>

### Installing

To install this tool clone this repository to your PC.

```
~$ git clone https://github.com/LSEA
```

## Running and using tool

Firstly, you need to prepare tsv-file from GWAS summary statistics with the following structure: <br> 
<table>
  <tr>
    <td>CHR</td>
    <td>COORDINATE</td>
    <td>RSID</td>
    <td>REF</td>
    <td>ALT</td>
    <td>PVAL</td>
  </tr>
</table>

To launch this tool you will also need to specify path to plink and snpeff directories.

## Example of usage
```
~$ python3 LSEA.py -af <input tsv-file> -sn <path to SNPeff> -pld <path to plink> -bf <bfile for plink> -p
```
This command will apply LSEA algorithm to the input file and will generate tsv-file with the following structure: 
<table>
  <tr>
    <td>BIOCARTA_INTRINSIC_PATHWAY</td>
    <td>2.0446237642438122e-14</td>
    <td>2.2517441515617103e-10</td>
    <td>(17776, 11, 36, 6, 'F11;FGB;FGA;F5;FGG;KLKB1')</td>
  </tr>
</table>
The first column contains the name of the set, the second and the third represent p-value and corrected q-value of hypergeometric test, the last coloumn includes information about total number of independant loci, number of loci in quiery, number of loci in gene set, number of loci common for quiery and gene set and, finally, the genes list.<br> 
Note that the genes list could be smaller then the number of common loci, because only indepedant loci are counted for analysis. <br>
-p (--precompiled flag) points that precompiled universe of independant loci based on UK Biobank data is used.<br>
Information about HLA-locus is excluded from analisys due to high ambiguity of LD-scores within the HLA-locus.

## Information about flags: 
  -af <input.tsv> Input file in tsv-format <br>
  -vf <imput.vcf> Annotated vcf-file if it is already prepared <br>
  -pl <imput.clumped> Plink result of clumping (.clumped file) if it is already prepared <br>
  --precompiled, -p Use precompiled loci <br>
  -sn <path to SNPeff directory> Path to SNPeff <br>
  -g <genome> Flag for specifying genome for SNPeff annotattion <br>
  -pld <path to plink directory> Path to plink <br>
  -bf <bfile> Bfile for plink <br>

## Creating your own universe:  
If you don't want to use precompiled universe of independant loci you can use options for creating your own universe based on GWAS summary statictics files. Use -cu (--create_universe) option to create universe of independant from your data. For this function you have to prepare results of clumping for your GWAS data (obtain .clumped file). If you have multiple files (e.g. for different phenotypes use -cld <directory> flag for specifying directory with clumped files)


## Author

* **Anton Shikov** - *Initial work* - [anton-shikov](https://github.com/anton-shikov)


## License

This project is free and available for everyone.

