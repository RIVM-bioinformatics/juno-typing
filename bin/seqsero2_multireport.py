import pandas as pd
import argparse
import os

# Extract info from one report
def extract_from_seqsero(input):
    assert os.path.exists(input), "One or more input files do not exist. Please make sure that you provided the right paths to your input files."
    seqsero_report = pd.read_csv(input, sep='\t')
    # seqsero_report = seqsero_report.iloc[:,[0,3,4,5,6,7,8,9,10]]
    return seqsero_report

# Create and save multi_report
def main(args):
    multi_report = list( map(extract_from_seqsero, args.input) )
    multi_report = multi_report[0].append([multi_report[i] for i in range(1, len(multi_report))], ignore_index=True)
    multi_report["Sample name"] = [ file_n.split("_")[0] for file_n in multi_report["Sample name"].tolist() ] # In file names, remove everything after the first underscore
    multi_report.to_csv(args.out_report, index = False)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', nargs='+', 
                        help="Input files (tsv) produced by seqsero2 that need to be combined into one multireport")
    parser.add_argument('-o', '--out_report', type=str, 
                       help="Path (relative or absolute) and name of the output file that contains the multireport")
    main(parser.parse_args())