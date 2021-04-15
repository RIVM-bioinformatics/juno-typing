import pandas as pd
import argparse
import os
import re

# Extract sample
def extract_sample_name(input):
    assert type(input) is str, "The input directory must be a string"
    search_name = re.match('(.*)/mlst7/(.*)/results_tab.tsv', input)
    if(input):
        return str(search_name.group(2))

# Extract info from one report
def extract_from_mlst7(input):
    assert os.path.exists(input), "One or more input files do not exist. Please make sure that you provided the right paths to your input files."
    sample_name = extract_sample_name(input)
    mlst7_report = pd.read_csv(input, sep='\t')[['Locus', 'Allele']]
    mlst7_report.rename(columns = {'Allele':sample_name}, inplace = True)
    mlst7_report = mlst7_report.set_index('Locus')
    return mlst7_report
    
# Create and save multi_report
def main(args):
    multi_report = list( map(extract_from_mlst7, args.input) )
    for i in range(1, len(multi_report)):
        if set(multi_report[0].index) == set(multi_report[i].index):
            multi_report[0] = multi_report[0].join(multi_report[i])
        else:
            multi_report[0] = multi_report[0].append(multi_report[i])
    multi_report[0].to_csv(args.out_report, index = True)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', nargs='+', 
                        help="Input files (tsv) produced by CGE-MLST7 that need to be combined into one multireport")
    parser.add_argument('-o', '--out_report', type=str, 
                       help="Path (relative or absolute) and name of the output file that contains the multireport")
    main(parser.parse_args())
    