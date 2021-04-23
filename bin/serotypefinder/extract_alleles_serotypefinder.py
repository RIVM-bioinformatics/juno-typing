#!/usr/bin/env python3
"""Extract information from serotypefinder.

Usage:
  extract_alleles_seotypefinder.py <file_name>

<file_name> is the full name (with path) to the data.json file generated
while using SerotypeFinder (CGE) for E. coli samples. Output will be a
data frame with the information of the two O-alleles and the H-allele.
For example:

           wzx_186_AB627352_O183 wzy_134_AB812069_O169/O183 fliC_308_AY250001_H18
gene                       wzx                        wzy                  fliC
serotype                  O183                  O169/O183                   H18
identity                   100                      97.25                 98.32
coverage                   100                        100                   100

  ...
"""

import argparse
import json
import pandas as pd

def main(data_json_path, output_csv):
    result_file = open(str(data_json_path)).read()
    result = json.loads(result_file)
    h_type = result['serotypefinder']['results']['H_type']
    o_type = result['serotypefinder']['results']['O_type']
    if h_type != 'No hit found':
        h_type = pd.DataFrame.from_dict(h_type)
    else:
        h_type = pd.DataFrame.from_dict({"no_H_hit": {"gene": "No hit found", "serotype": "No hit found", "identity": 0, "HSP_length": 0, "template_length": 0, "position_in_ref": "NA", "contig_name": "NA", "positions_in_contig": "NA", "accession": "NA", "coverage": 0, "hit_id": "NA"}})
    if o_type != 'No hit found':
        o_type = pd.DataFrame.from_dict(o_type)
    else:
        o_type = pd.DataFrame.from_dict({"no_O_hit": {"gene": "No hit found", "serotype": "No hit found", "identity": 0, "HSP_length": 0, "template_length": 0, "position_in_ref": "NA", "contig_name": "NA", "positions_in_contig": "NA", "accession": "NA", "coverage": 0, "hit_id": "NA"}})
    alleles = pd.concat([h_type, o_type], axis = 1)
    alleles.to_csv(output_csv)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract information from serotypefinder results")
    parser.add_argument("file_path", type=str, 
                       help="Name (and path) of the json file containing the results from SerotypeFinder")
    parser.add_argument("output_path", type=str, 
                       help="Name (and path) of the output csv file where the results will be written")
    main(parser.parse_args().file_path, parser.parse_args().output_path)
