#!/usr/bin/env python3
"""Multireport from seroba
"""

import argparse
import pandas as pd
import os


def main(list_files, output_file):
    assert all([os.path.isfile(file) for file in list_files]), "One or more of the specified input files do not exist"
    names = ['Sample', 'Serotype', 'Contamination']
    results = [pd.read_csv(file, sep = '\t', header = None, names = names) for file in list_files]
    results_df = results[0].append(results[1:])
    results_df.to_csv(output_file, index = False)
    print(results_df)
    return results_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a summary report for the serotyping result of one or more E. coli samples")
    parser.add_argument('input_paths', type=str, nargs='+',
                        help="One or more input files (with path) containing results of serotyper finder (result_serotype.csv)")
    parser.add_argument('-o', '--out', 
                        type=str,
                        help="Name (with path) of desired output csv file")
    print("Input file(s): {}".format(parser.parse_args().input_paths))
    print("Output file: {}".format(parser.parse_args().out))
    main(parser.parse_args().input_paths, parser.parse_args().out)



