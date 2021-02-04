#############################################################################
##### Title: Generate list of species used in this run                  #####
##### Author: Alejandra Hernández Segura                                #####
#############################################################################

import argparse
from os import path
import yaml

def main(args):
    assert path.exists(args.sample_sheet), "Argument must be an existing file."
    
    samples = {}
    with open("sample_sheet.yaml") as sample_sheet_file:
        SAMPLES = yaml.safe_load(sample_sheet_file) 

    # Species that are analyzed in this run
    SPECIES = set()
    for sample in SAMPLES:
        sp = SAMPLES[sample]["species"]
        SPECIES.add(sp)

    for spp in SPECIES:
        print(spp)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("sample_sheet", 
                       help="Sample sheet with all the input files are located")
    main(parser.parse_args())
