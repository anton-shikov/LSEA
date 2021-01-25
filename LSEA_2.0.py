import argparse
import subprocess
import sys
import os
import csv
import json
from collections import defaultdict
from scipy.stats import hypergeom
from rpy2 import robjects
from utils import create_snp_dict, get_intersected_genes, count_intervals


def run_plink(plink_path, bfile_path, tsv_file, input_dict):
    print('Calculating independent loci with PLINK')
    tsv_plink = os.getcwd() + "/" + \
        tsv_file.split('/')[-1].split('.')[0] + '_for_plink.tsv'
    out_plink = os.getcwd() + "/" + \
        tsv_file.split('/')[-1].split('.')[0]
    with open(tsv_plink, 'w', newline='') as csvfile:
        my_writer = csv.writer(csvfile, delimiter='\t')
        init_row = ["SNP", "Chr", "Pos", "P"]
        my_writer.writerow(init_row)
        for snp in input_dict:
            # Only extract necessary info from dictionary corresponding to csv
            row = [snp] + input_dict[snp][0:3]
            my_writer.writerow(row)
    subprocess.call('{0}/plink --bfile {1} --clump {2} --clump-field P --clump-p1 1e-05 --clump-r2 0.1 --clump-snp-field SNP --clump-kb 500 --out {3} --allow-no-sex \
      2> {4}'.format(plink_path, bfile_path, tsv_plink,
                     out_plink, 'PLINK_clumping.log'), shell=True)
    return out_plink + ".clumped"


def normalize_path(path):
    if path[-1] == "/":
        path = path[:-1]
    return path


def make_bed_file(clumped_file):
    with open("./clumps.bed", 'w', newline='') as bed_file:  # Here we write to new file
        my_writer = csv.writer(bed_file, delimiter='\t')
        with open(clumped_file, 'r', newline='') as cl_file:  # Our result of clumping (SNPs sets)
            my_reader = csv.reader(cl_file, delimiter='\t')
            for row in my_reader:
                if len(row) != 0:
                    # What to do with NONE?
                    row = list(filter(lambda x: len(x) !=
                                      0 and x != " ", row[0].split(" ")))
                    if row[0] == "CHR":  # Skip first row
                        continue
                    bed_row = [row[0], int(row[3]) -
                               interval, int(row[3]) + interval]
                    my_writer.writerow(bed_row)
    # Sort file
    subprocess.call(
        "bedtools sort -i clumps.bed > clumps_sorted.bed", shell=True)
    # Merge regions
    subprocess.call(
        "bedtools merge -i clumps_sorted.bed > merged.bed", shell=True)
    with open("./merged_fixed_size.bed", 'w', newline='') as bed_file:  # Here we write to new file
        my_writer = csv.writer(bed_file, delimiter='\t', lineterminator='\n')
        with open('merged.bed', 'r', newline='') as inter:
            my_reader = csv.reader(inter, delimiter='\t')
            for row in my_reader:
                middle_point = int(row[1]) + (int(row[2]) - int(row[1])) // 2
                new_row = [row[0], middle_point -
                           interval, middle_point + interval]
                my_writer.writerow(new_row)
    # Add numeration to intervals
    subprocess.call(
        "awk {'print $0\"\t\"FNR'} merged_fixed_size.bed > merged_with_line_numbers.bed", shell=True)


def p_val_for_gene_set(n_big, k_big, n, k):
    return hypergeom.sf(k, n_big, k_big, n)


def calcuate_qvals(pvals):
    x = robjects.r('''
            f <- function(pvals) {
                p.adjust(pvals, method = "fdr")
            }
            ''')
    r_f = robjects.globalenv['f']
    v = robjects.FloatVector(pvals)
    return list(r_f(v))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='LSEA')
    parser.add_argument('-pldir', help='Path to PLINK directory', metavar='dir',
                        type=str, required=False)
    parser.add_argument('-b_file', help='Bfile for PLINK', metavar='path',
                        type=str, required=False)
    parser.add_argument('-tsv', help='Input file in tsv-format', metavar='file',
                        type=str, required=True)
    parser.add_argument('-interval', help='Size of interval taken from each clumping center (Default: 500000)',
                        metavar='int', type=int, required=False, default=500000)
    parser.add_argument('-use_clumped', help='Path to .clumped file if already exists',
                        metavar='path', type=str, required=False)
    parser.add_argument('-gene_file', help='Path to file with genes',
                        metavar='path', type=str, required=True)
    parser.add_argument('-qval_threshold', help='Q-value threshold for output (Default: 0.1)',
                        metavar='float', type=float, required=False, default="0.1")
    parser.add_argument('-column_names', help='Column names for input tsv. Should be 4 values - for chromosome, position, id and pval (Default: chr, pos, id, p)',
                        metavar='name', nargs=4, type=str, required=False)
    parser.add_argument('-json', help='Json file with universe',
                        metavar='path', type=str, required=True)
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    tsv_file = args.tsv
    path_to_plink_dir = args.pldir
    if path_to_plink_dir is not None:
        path_to_plink_dir = normalize_path(args.pldir)
    path_to_bfile = args.b_file
    interval = args.interval
    clumped_file = args.use_clumped
    path_to_gene_file = args.gene_file  # TODO: add default file
    qval_thresh = args.qval_threshold
    col_names = args.column_names
    json_file = args.json
    if col_names is None:
        col_names = ["chr", "pos", "id", "p"]
    input_dict = create_snp_dict(tsv_file, col_names)
    if clumped_file is None:
        clumped_file = run_plink(
            path_to_plink_dir, path_to_bfile, tsv_file, input_dict)

    base = json.load(open(json_file, "r"))
    interval_counts_for_universe = base["interval_counts"]

    make_bed_file(clumped_file)
    n_intervals = sum(1 for line in open('merged_with_line_numbers.bed'))
    genes = get_intersected_genes(
        "./merged_with_line_numbers.bed", path_to_gene_file, "inter.tsv", False)
    msig_dict = base["msig_dict"]
    interval_counts = count_intervals(msig_dict, genes)


    pvals = []
    for w in sorted(interval_counts, key=interval_counts.get, reverse=True):
        if interval_counts[w] != 0:
            pvals.append(p_val_for_gene_set(base["universe_intervals_number"], interval_counts_for_universe[w], n_intervals, interval_counts[w]))

    qvals = calcuate_qvals(pvals)
    with open("./result.tsv", 'w', newline='') as file:
        my_writer = csv.writer(file, delimiter='\t')
        my_writer.writerow(["Gene_set", "p-value", "q-value"])
        i = 0
        for w in sorted(interval_counts, key=interval_counts.get, reverse=True):
            if interval_counts[w] != 0:
                if qvals[i] <= qval_thresh:
                    row = [w, pvals[i], qvals[i]]
                    my_writer.writerow(row)
                i += 1
