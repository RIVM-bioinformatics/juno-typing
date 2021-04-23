#!/bin/python

"""Get species name (formated for CGE-MLST).

Usage:
    get_mlst7_db_name.py </path/to/kmerfinder/results/data.json>
"""

import argparse
import pathlib
import json
import yaml
import pandas as pd

with open("files/dictionary_correct_species.yaml") as translation_yaml:
    translation_tbl = yaml.safe_load(translation_yaml)

def parse_species_result(text_file):
    species_file = open(text_file)
    species = species_file.read()
    species = species.split()[0][0] + species.split()[1]
    species_file.close()
    return species.lower()

def choose_right_db(species, translation_tbl):
    if species in translation_tbl.keys():
        species = translation_tbl[species]
    return species

def main(args):
    assert args.species_result.is_file(), "kmerfinder_res must be an existing file."
    species = parse_species_result(args.species_result)
    print(choose_right_db(species, translation_tbl))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("species_result", type=pathlib.Path, 
                       help="text file containing the top hit of KmerFinder.")
    main(parser.parse_args())

