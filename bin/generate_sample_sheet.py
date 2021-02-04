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

fq_pattern = re.compile("(.*?)(?:_S\d+_|_S\d+.|_|\.)R?(1|2)(?:_.*\.|\..*\.|\.)f(ast)?q(\.gz)?")
fasta_pattern = re.compile("(.*?)(?:_.*\.|\..*\.|\.).fasta?")

def main(args):
    assert args.dir.is_dir(), "Argument must be a directory."

    samples = {}
    for file_ in args.dir.iterdir():
        if file_.is_dir():
            continue
        match = fq_pattern.fullmatch(file_.name)
        if match:
            sample = samples.setdefault(match.group(1), {})
            sample["R{}".format(match.group(2))] = str(file_)
    if len(samples) == 0 :
        for file_ in args.dir.iterdir():
            if file_.is_dir():
                continue
            match = fasta_pattern.fullmatch(file_.name)
            if match:
                sample = samples.setdefault(match.group(1), {})
                sample["assembly"] = str(file_)

    if args.metadata is not None :
        assert args.metadata.is_file(), "Provided metadata file does not exist"
        # Load species file
        species_file = pd.read_csv(args.metadata, index_col = 0, dtype={'Sample':str})
        species_file.index = species_file.index.map(str)
        species_file = species_file.transpose().to_dict()
        for sample_name in samples :
            try :
                samples[sample_name]["species"] = species_file[sample_name]["Genus"].strip().lower()[0] + species_file[sample_name]["Species"].strip()
            except:
                pass
                    
    print(yaml.dump(samples, default_flow_style=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dir", type=pathlib.Path, 
                       help="Directory where input files are located")
    parser.add_argument("--metadata", type=pathlib.Path, 
                       help=".csv file containing at least 3 columns: 'Sample', 'Genus' and 'Species'")
    main(parser.parse_args())
