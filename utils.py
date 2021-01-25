from collections import defaultdict
import csv
import subprocess

def count_intervals(dict, genes):
    res = defaultdict(int)
    for name in dict.keys():
        gene_list = dict[name]  # List of genes for a trait
        interval_set = set()  # Count every interval once
        for gene in genes:
            if gene in gene_list:  # If gene belongs to the trait
                # Add every corresponding interval to set
                for interval in genes[gene]:
                    interval_set.add(int(interval))
        res[name] = len(interval_set)
    return res

# Extract gene names that are in intersection
def get_intersected_genes(path_to_bed, path_to_gene_file, out_file, skip_intersect):
    genes = defaultdict(list)  # Gene -> interval id
    if not skip_intersect:
        subprocess.call("bedtools intersect -a {0} -b {1} -loj > {2}".format(
            path_to_bed, path_to_gene_file, out_file), shell=True)
    with open(out_file, 'r', newline='') as inter:  # Our result of clumping (SNPs sets)
        my_reader = csv.reader(inter, delimiter='\t')
        for row in my_reader:
            genes[row[-1]].append(row[3])
    return genes

def create_snp_dict(tsv_file, names):
    input_dict = defaultdict(list)
    names = list(map(lambda x: x.lower(), names))
    with open(tsv_file, 'r', newline='') as csvfile:  # Convert tsv to dictionary
        my_reader = csv.reader(csvfile, delimiter='\t')
        csv_headings = list(map(lambda x: x.lower(), next(my_reader)))
        try:  # We find necessary indices from header
            chr_index = csv_headings.index(names[0])
            pos_index = csv_headings.index(names[1])
            id_index = csv_headings.index(names[2])
            p_index = csv_headings.index(names[3])
        except ValueError:
            print("Check that your tsv file has headers and they are correct!")
            exit(1)
        for row in my_reader:  # Start from second row
            input_dict[row[id_index]] = [row[chr_index]] + \
                [row[pos_index]] + [row[p_index]]  # Name -> Chromosome, pos, pval
    return input_dict

def get_msig_dict(path):
    d = defaultdict(list)
    with open(path, 'r', newline='') as db:
        my_reader = csv.reader(db, delimiter='\t')
        for row in my_reader:
            for gene in row[2:]:
                d[row[0]].append(gene)
    return d
