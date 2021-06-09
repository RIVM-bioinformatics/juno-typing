#!/usr/bin/env python3
"""Extract information from serotypefinder.

Usage:
  extract_alleles_seotypefinder.py <result_serotype1> <result_serotype2> ... <result_serotypeN>

<result_serotypeX> is the full name (with path) to the result_serotype.csv file generated
while using the script extract_alleles_serotypefinder.py. Output will be a
data frame with the information of the two O-alleles and the H-allele for each sample.
For example:

    wzx   wzy wzm wzt fli
0  O157  O157          H7
1   O63                H6

  ...
"""

import argparse
import re
import pandas as pd
import os
from itertools import chain


def find_allele(serotype_df, allele_string):
    allele = []
    for column_name in serotype_df.columns.tolist():
        result_allele = serotype_df.loc['serotype', column_name]
        if allele_string in column_name:
            if result_allele not in allele:
                allele.append(result_allele)
    return '/'.join(allele)


def get_sample_serotype(serotype_csv):
    # Read csv file
    result_csv = pd.read_csv(serotype_csv, index_col = 0)
    # Extract information from individual samples
    genes = ['wzx', 'wzy', 'wzm', 'wzt', 'fl']
    result_summary = []
    for gene in genes:
        result_summary.append(find_allele(result_csv, gene))
    return result_summary

def report_o_type(row_df):
    row_df = row_df[['wzx', 'wzy', 'wzm', 'wzt']]
    # Split string in case any locus already consists of two alleles
    reported_alleles = [str(item).split('/') for item in row_df if item is not '']
    if len(reported_alleles) == 0:
        return 'Error! No O type found'
    else:
        # Keep only unique alleles and report them all together (separated by a /)
        final_o_type = list(chain(*reported_alleles))
        # By request of Kim van der Zwaluw, only the numbers of the serotype are reported not the letters after it (e.g. O128ac reported as O128)
        final_o_type = [re.findall("\d+", sample) for sample in final_o_type]
        final_o_type = list(set(chain(*final_o_type)))
        # If there is more than one serotype, sort them from smaller to larger (requested by Kim)
        if len(final_o_type) > 1:
            ranking_o_types = sorted(range(len(final_o_type)), key=lambda k: int(final_o_type[k]))
            final_o_type = [final_o_type[ind] for ind in ranking_o_types]
        # Give back the serotype name
        final_o_type = ["O{}".format(final_o_type[i]) for i in range(0, len(final_o_type))]
        return '/'.join(final_o_type)

def merge_multiple_reports(list_files, sample_names):
    results = list(map(get_sample_serotype, list_files))
    results_df = pd.DataFrame(results)
    results_df.columns = ['wzx', 'wzy', 'wzm', 'wzt', 'fli']
    results_df.index = [sample_names]
    results_df['O type'] = [report_o_type(row) for index, row in results_df.iterrows()]
    results_df['H type'] = results_df['fli']
    results_df.loc[results_df['H type']=='', 'H type'] = 'Error! No H type found'
    return results_df

def main(list_files, output_file):
    sample_names = list(map(os.path.dirname, list_files))
    sample_names = list(map(os.path.basename, sample_names))
    results_df = merge_multiple_reports(list_files, sample_names)
    results_df.to_csv(output_file, index = True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a summary report for the serotyping result of one or more E. coli samples")
    parser.add_argument('input_paths', type=str, nargs='+',
                        help="One or more input files (with path) containing results of serotyper finder (result_serotype.csv)")
    parser.add_argument('out_file', type=str,
                        help="Name (with path) of desired output csv file")
    main(parser.parse_args().input_paths, parser.parse_args().out_file)



