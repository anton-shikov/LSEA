import argparse
import sys
import csv
import json
from utils import create_snp_dict, get_intersected_genes, count_intervals, get_msig_dict

def create_universe(dict, interval):
    with open("./universe.bed", 'w', newline='') as bed_file:  # Here we write to new file
        my_writer = csv.writer(bed_file, delimiter='\t')
        id = 1
        for snp in dict:
            info = dict[snp]
            # Chromosome, start and end coordinates, id of an interval
            start = max(0, int(info[1]) - interval)
            bed_row = [info[0], start,
                       int(info[1]) + interval, id]
            id += 1
            my_writer.writerow(bed_row)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Universe generator for LSEA')
    parser.add_argument('-column_names', help='Column names for input tsv. Should be 4 values - for chromosome, position, id and pval (Default: chr, pos, id, p)',
                    metavar='name', nargs=4, type=str, required=False)
    parser.add_argument('-path', help='Path to file from which universe will be created',
                        metavar='path', type=str, required=True)
    parser.add_argument('-interval', help='Size of interval taken from each clumping center (Default: 500000)',
                        metavar='int', type=int, required=False, default=500000)
    parser.add_argument('-gene_file', help='Path to file with genes',
                        metavar='path', type=str, required=True)
    parser.add_argument('-msig_path', help='Path to file with set/gene relationships',
                        metavar='path', type=str, required=True)
    parser.add_argument('-skip_intersect', help='Use if inter2.tsv is already created',
                        action = "store_true", required=False)
    parser.add_argument('-o', help='Output path for json',
                        metavar='path', type=str, required=True)
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    path = args.path
    interval = args.interval
    col_names = args.column_names
    path_to_gene_file = args.gene_file  # TODO: add default file
    msig_path = args.msig_path
    skip_intersect = args.skip_intersect
    out_path = args.o
    col_names = args.column_names
    
    if col_names is None:
        col_names = ["chr", "pos", "id", "p"]
    input_dict = create_snp_dict(path, col_names)
    msig_dict = get_msig_dict(msig_path)
    create_universe(input_dict, interval)
    genes_in_universe = get_intersected_genes(
        "./universe.bed", path_to_gene_file, "inter2.tsv", skip_intersect)
    interval_counts_for_universe = count_intervals(
        msig_dict, genes_in_universe)
    out_dict = {"universe_intervals_number" : len(input_dict), "interval_counts" : interval_counts_for_universe, "msig_dict" : msig_dict}

    base = open(out_path, "w")
    json.dump(out_dict, base)
    base.close()
