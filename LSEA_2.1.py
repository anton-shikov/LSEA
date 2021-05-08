import argparse
import subprocess
import sys
import os
import csv
import json
import re
from collections import defaultdict
from scipy.stats import hypergeom
from rpy2 import robjects
from utils import get_overlapping_features, count_intervals


def tsv_len(tsv_file):
    with open(tsv_file) as f:
        for i, l in enumerate(f):
            pass
    return i

def get_h2(tsv_file, ldsc_dir, column_names, ldsc_id):
    global h2
    subprocess.call(
        '{0}/munge_sumstats.py --N {1} --snp {2} --p {3} --out for_h2 --sumstats {4}'.
            format(ldsc_dir, N, ldsc_id, column_names[3], tsv_file), shell=True
    )
    subprocess.call(
        '{0}/ldsc.py --h2 for_h2.sumstats.gz --ref-ld {0}/eur_w_ld_chr_custom/genotypes_genome_hapgen.MAF_higher_0.05.with_cms \
        --w-ld {0}/eur_w_ld_chr_custom/genotypes_genome_hapgen.MAF_higher_0.05.with_cms --out h2'.format(ldsc_dir),
        shell=True
    )
    pattern = re.compile('Total Observed scale h2:')
    separator = ":"
    for line in open('h2.log'):
        for match in re.finditer(pattern, line):
            informative = line.partition(separator)[2]
            h2 = float(re.findall(r"[-+]?\d*\.\d+|\d+", informative)[0])
    return h2


def get_cutoff(N, M, h2, tsv_file):
    K = tsv_len(tsv_file)
    lin = 1.547e-01
    b1 = 1.244e-01
    b2 = 8.973e-02
    b3 = 1.899e-01
    b4 = -3.372e-04
    log_p_cutoff = lin * (N ** b1 + K ** b2 + h2 ** b3 - M * b4 +
    (N ** b1) * (K ** b2) + (N ** b1) * (h2 ** b3) +
    (N ** b1) * (M * b4) + (K ** b2) * (h2 ** b3) + (h2 ** b3) * (M * b4) +
    (N ** b1) * (K ** b2) * (h2 ** b3) + (N ** b1) * (K ** b2) * (M * b4) + (K ** b2) * (h2 ** b3) * (M * b4) +
    (N ** b1) * (K ** b2) * (h2 ** b3) * (M * b4))
    return 10**-(log_p_cutoff)


def run_plink(plink_path, bfile_path, tsv_file, input_dict, n, m, h2):
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
    subprocess.call('{0}/plink --bfile {1} --clump {2} --clump-field P --clump-p1 {3} --clump-r2 0.1 --clump-snp-field SNP --clump-kb 500 --out {4} --allow-no-sex --allow-extra-chr \
      2> {5}'.format(plink_path, bfile_path, tsv_plink, get_cutoff(n, m, h2, tsv_file),
                     out_plink, 'PLINK_clumping.log'), shell=True)
    return out_plink + ".clumped"


def normalize_path(path):
    if path[-1] == "/":
        path = path[:-1]
    return path


def get_snp_info(tsv_file, names):
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


def make_bed_file(clumped_file, interval):
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
    parser.add_argument('-input', help='Input file in TSV format. The file should contain at least four columns: chromosome name, variant position, variant ID, and GWAS p-value. The file MUST have a header row.', metavar='file',
                        type=str, required=True)
    parser.add_argument('-universe', help='Json file with universe',
                        metavar='path', type=str, required=True)
    parser.add_argument('-out', help='Output file name (Default: result.tsv)',
                        metavar='name', type=str, required=False, default='result.tsv')
    parser.add_argument('-n', help='Number of individuals in the GWAS analysis (sample size)',
                        metavar='int', type=int, required=True)
    parser.add_argument('-m', help='Expected number of causal SNPs in the dataset (Default: 30)',
                        metavar='int', type=int, required=False, default=30)
    parser.add_argument('-h2', help='Expected number of total variance explained by SNPs (SNP-based heritability). If not specified, will be estimated using LDSC',
                        metavar='float', type=float, required=False)
    parser.add_argument('-plink_dir', help='Path to a directory with PLINK executable', metavar='dir',
                        type=str, required=False)
    parser.add_argument('-bfile', help='Genotypes in PLINK .bed format that will be used for LD-based clumping. Only file prefix should be given.', metavar='prefix',
                        type=str, required=False)
    parser.add_argument('-ldsc_dir', help='Path to LDSC directory (provide only if h2 is unknown)', metavar='dir',
                        type=str, required=False)
#    parser.add_argument('-interval', help='Size of interval taken from each clumping center (Default: 500000)',
#                        metavar='int', type=int, required=False, default=500000)
#    parser.add_argument('-use_clumped', help='Path to .clumped file if already exists',
#                        metavar='path', type=str, required=False)
#    parser.add_argument('-gene_file', help='Path to file with genes',
#                        metavar='path', type=str, required=True)
    parser.add_argument('-qval_threshold', help='Q-value threshold for output (Default: 0.1)',
                        metavar='float', type=float, required=False, default="0.1")
    parser.add_argument('-column_names', help='Column names for input TSV. These names will be used if the column names do not match the default: chr, pos, id, p. Names should be given in the same order (chromosome, position, ID, p-value)',
                        metavar='name', nargs=4, type=str, required=False)
    parser.add_argument('-alt_id', help='Column name for alternative variant ID matching your LD scores estimations (use only if variant IDs for LD scores do not match variant IDs in your association data (Default: id)',
                        metavar='name', type=str, required=False)
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    tsv_file = args.input
    path_to_plink_dir = args.plink_dir
    if path_to_plink_dir is not None:
        path_to_plink_dir = normalize_path(args.plink_dir)
    path_to_bfile = args.bfile
    qval_thresh = args.qval_threshold
    col_names = args.column_names
    ldsc_name = args.alt_id
    json_file = args.universe
    ldsc_path = args.ldsc_dir
    samples_number = args.n
    causal_number = args.m
    heritability = args.h2
    out_name = args.out
    if col_names is None:
        col_names = ["chr", "pos", "id", "p"]
    if ldsc_name is None:
        ldsc_name = col_names[2]
    if heritability is None:
        heritability = get_h2(tsv_file, ldsc_path, col_names, ldsc_name)
    input_dict = get_snp_info(tsv_file, col_names)
    clumped_file = run_plink(path_to_plink_dir, path_to_bfile, tsv_file, input_dict, samples_number, causal_number, heritability)

    universe= json.load(open(json_file, "r"))
    interval = universe["interval"]
    with open('./features.bed', 'w') as feature_file:
        for feature in universe["features"]:
            bed_line = '\t'.join(universe["features"][feature])
            feature_file.write(f'{bed_line}\n')
    interval_counts_for_universe = universe["interval_counts"]

    make_bed_file(clumped_file, interval)
    n_intervals = sum(1 for line in open('merged_with_line_numbers.bed'))
    target_features = get_overlapping_features("./merged_with_line_numbers.bed", './features.bed', "inter.tsv")
    feature_set = universe["gene_set_dict"]
    interval_counts = count_intervals(feature_set, target_features)

    pvals = []
    for w in sorted(interval_counts, key=interval_counts.get, reverse=True):
        if interval_counts[w] != 0:
            pvals.append(p_val_for_gene_set(universe["universe_intervals_number"], 
                interval_counts_for_universe[w], n_intervals, interval_counts[w]))

    qvals = calcuate_qvals(pvals)
    with open(f"./{out_name}", 'w', newline='') as file:
        my_writer = csv.writer(file, delimiter='\t')
        my_writer.writerow(["Gene_set", "p-value", "q-value"])
        for i, w in enumerate(sorted(interval_counts, key=interval_counts.get, reverse=True)):
            if interval_counts[w] != 0:
                if qvals[i] <= qval_thresh:
                    row = [w, pvals[i], qvals[i]]
                    my_writer.writerow(row)
