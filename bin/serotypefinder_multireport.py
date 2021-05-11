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
import pandas as pd
import os
import numpy as np


def find_allele(serotype_df, allele_string):
    allele = []
    for column_name in serotype_df.columns.tolist():
        result_allele = serotype_df.loc['serotype', column_name]
        if allele_string in column_name:
            if result_allele not in allele:
                allele.append(result_allele)
    return '/'.join(allele)


def get_sample_serotype(serotype_csv):
    #Read csv file
    result_csv = pd.read_csv(serotype_csv, index_col = 0)
    # Extract information from individual samples
    genes = ['wzx', 'wzy', 'wzm', 'wzt', 'fl']
    result_summary = []
    for gene in genes:
        result_summary.append(find_allele(result_csv, gene))
    return result_summary

def main(list_files, output_file):
    sample_names = list(map(os.path.dirname, list_files))
    sample_names = list(map(os.path.basename, sample_names))
    results = list(map(get_sample_serotype, list_files))
    results_df = pd.DataFrame(results)
    results_df.columns = ['wzx', 'wzy', 'wzm', 'wzt', 'fli']
    results_df.index = [sample_names]
    # Conditions for finding O type
    conditions = [
        (results_df['wzx'] != '') & (results_df['wzx'] == results_df['wzy']),
        (results_df['wzx'] != '') & (results_df['wzy'] == ''),
        (results_df['wzx'] == '') & (results_df['wzy'] != ''),
        (results_df['wzx'] == 'O50/O2') & (results_df['wzy'] == 'O50'),
        (results_df['wzx'] == 'O50/O2') & (results_df['wzy'] == 'O2'),
        (results_df['wzx'] == 'O169') & (results_df['wzy'] == 'O169/O183'),
        (results_df['wzx'] == 'O183') & (results_df['wzy'] == 'O169/O183'),
        (results_df['wzx'] == 'O123') & (results_df['wzy'] == 'O123/O186'),
        (results_df['wzx'] == 'O186') & (results_df['wzy'] == 'O123/O186'),
        (results_df['wzx'] == 'O17/O77') & (results_df['wzy'] == 'O17/O77'),
        (results_df['wzx'] == 'O44') & (results_df['wzy'] == 'O17/O44'),
        (results_df['wzx'] == 'O118/0151') & (results_df['wzy'] == 'O118/0151'),
        (results_df['wzx'] == 'O118/0151') & (results_df['wzy'] == 'O118/0151'),
        (results_df['wzx'] == 'O13') & (results_df['wzy'] == 'O13/O135'),
        (results_df['wzx'] == 'O164') & (results_df['wzy'] == 'O124'),
        (results_df['wzx'] == 'O90') & (results_df['wzy'] == 'O127'),
        (results_df['wzx'] == 'O134') & (results_df['wzy'] == 'O46'),
        (results_df['wzx'] == 'O128ab') & (results_df['wzy'] == 'O128ac'),
        (results_df['wzm'] == 'O162') & (results_df['wzy'] == 'O101'),
        (results_df['wzx'] != '') & (results_df['wzx'] != results_df['wzy']),
        (results_df['wzx'] == '') & (results_df['wzy'] == '') & (results_df['wzm'] == '') & (results_df['wzt'] == ''),
        (results_df['wzx'] == '') & (results_df['wzy'] == '') & (results_df['wzm'] == results_df['wzt']),
        (results_df['wzx'] == '') & (results_df['wzy'] == '') & (results_df['wzm'] != results_df['wzt'])
        ]
    values  = [results_df['wzx'], results_df['wzx'], results_df['wzy'], 
        'O50', 'O2', 'O169', 'O183', 'O123', 'O186', 'O17', 'O44', 'O118', 
        'O151', 'O135', 'O164', 'O127', 'O134', 'O128 (a,ab, abc)', 'O101', 
        'Undetermined', 'Error! No O type found', results_df['wzm'], 'Undetermined'] 
    results_df['O type'] = np.select(conditions, values)
    results_df['H type'] = results_df['fli']
    results_df.loc[results_df['H type']=='', 'H type'] = 'Error! No H type found'
    results_df.to_csv(output_file, index = True)
    print(results_df)
    return results_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a summary report for the serotyping result of one or more E. coli samples")
    parser.add_argument('input_paths', type=str, nargs='+',
                        help="One or more input files (with path) containing results of serotyper finder (result_serotype.csv)")
    parser.add_argument('out_file', type=str,
                       help="Name (with path) of desired output csv file")
    print("Input file(s): {}".format(parser.parse_args().input_paths))
    print("Output file: {}".format(parser.parse_args().out_file))
    main(parser.parse_args().input_paths, parser.parse_args().out_file)



