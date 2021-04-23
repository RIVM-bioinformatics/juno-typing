"""Generate sample sheet for Juno pipeline.

Usage:
  generate_sample_sheet.py <source_dir>

<source_dir> is a directory containing input fastq files with typical
filenames as used in the legacy (non-automated) process. Output will
be a sample sheet in YAML, for example:

  sample1_id:
    R1:
      path/to/sample_R1.fq.gz
    R2:
      path/to/sample_R2.fq.gz
  sample2_id:
  ...
"""

import argparse
import pathlib
import re
import yaml
import pandas as pd


fq_pattern = re.compile("(.*?)(?:_S\d+_|_S\d+.|_|\.)(?:p)?R?(1|2)(?:_.*\.|\..*\.|\.)f(ast)?q(\.gz)?")
fasta_pattern = re.compile("(.*?).fasta?")

def add_species_to_sample(sample, samples_dic, species_file):
    try: 
        samples_dic[sample].update(species_file[sample])
    except:
        pass

def add_metadata(metadata_file, samples_dic):
    assert metadata_file.is_file(), "Provided metadata file does not exist"
    # Load species file
    species_file = pd.read_csv(metadata_file, index_col = 0, dtype={'Sample':str})
    species_file.index = species_file.index.map(str)
    species_file['genus-abbr'] = species_file['Genus'].apply(lambda x: x[0])
    species_file['species-mlst7'] = species_file['genus-abbr'] + species_file["Species"]
    species_file['species-mlst7'] = species_file['species-mlst7'].apply(lambda x: x.lower())
    species_file = species_file.transpose().to_dict()
    # Update dictionary with species
    for sample_name in samples_dic :
        add_species_to_sample(sample_name, samples_dic, species_file)


def main(args):
    assert args.dir.is_dir(), "Argument must be a directory."

    # If input comes from Juno, the directory with fastq and 
    # the one with fasta should be different
    if(args.dir.joinpath("clean_fastq").is_dir()):
        fastq_dir = args.dir.joinpath("clean_fastq")
    else:
        fastq_dir = args.dir
    
    if (args.dir.joinpath("de_novo_assembly_filtered").is_dir()):
        fasta_dir = args.dir.joinpath("de_novo_assembly_filtered")
    else:
        fasta_dir = args.dir


    samples = {}
    for file_ in fastq_dir.iterdir():
        if file_.is_dir():
            continue
        match = fq_pattern.fullmatch(file_.name)
        if match:
            sample = samples.setdefault(match.group(1), {})
            sample["R{}".format(match.group(2))] = str(file_)
            sample['species-mlst7'] = "NotProvided"
    for file_ in fasta_dir.iterdir():
        if file_.is_dir():
            continue
        match = fasta_pattern.fullmatch(file_.name)
        if match:
            sample = samples.setdefault(match.group(1), {})
            sample["assembly"] = str(file_)
            sample['species-mlst7'] = "NotProvided"

    if args.metadata is not None :
        add_metadata(args.metadata, samples)
        
    print(yaml.dump(samples, default_flow_style=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dir", type=pathlib.Path, 
                       help="Directory where input files are located")
    parser.add_argument("--metadata", type=pathlib.Path, 
                       help=".csv file containing at least 3 columns: 'Sample', 'Genus' and 'Species'")
    main(parser.parse_args())

