import pandas as pd
import argparse
import os
import re
import json

# Extract sample
def extract_sample_name(input):
    assert type(input) is str, "The input directory must be a string"
    search_name = re.match('(.*)/mlst7/(.*)/data.json', input)
    if(input):
        return str(search_name.group(2))

# Extract info from one report
def extract_from_mlst7(input):
    assert os.path.exists(input), "One or more input files do not exist. Please make sure that you provided the right paths to your input files."
    sample_res = [extract_sample_name(input)]
    with open(input) as json_file:
        data = json.load(json_file)
    sample_res.append(data['mlst']['results']['sequence_type'])
    sample_res.append(data['mlst']['user_input']['organism'])
    allele_list = [allele for allele in data['mlst']['results']['allele_profile']]
    sample_res.append('-'.join(allele_list))
    alleles_res = []
    for allele_name in allele_list:
        alleles_res.append(data['mlst']['results']['allele_profile'][allele_name]['allele'])
    sample_res.append('-'.join(alleles_res))
    return sample_res
    
# Create and save multi_report
def main(args):
    multi_report = list( map(extract_from_mlst7, args.input) )
    colnames = ['Sample', 'ST_type', 'Scheme_used', 'genes_in_scheme', 'alleles']#'locus_1','locus_2','locus_3','locus_4','locus_5','locus_6','locus_7']
    multi_report = pd.DataFrame(multi_report, 
                                columns = colnames)
    multi_report.to_csv(args.out_report, index = False)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', nargs='+', 
                        help="Input files (json) produced by CGE-MLST7 that need to be combined into one multireport")
    parser.add_argument('-o', '--out_report', type=str, 
                       help="Path (relative or absolute) and name of the output file that contains the multireport")
    main(parser.parse_args())
    