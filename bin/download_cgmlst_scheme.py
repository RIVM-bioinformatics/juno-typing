import argparse
import base_juno_pipeline.helper_functions
import dask.bag as db
from datetime import datetime
import json
import os
import pathlib
import subprocess
import warnings
import yaml


class cgMLSTSchemes(base_juno_pipeline.helper_functions.JunoHelpers):
    """Class containing information and functions to download cgMLST schemes"""

    def __init__(self, 
                species_list,
                output_dir='output',
                threads=2,
                download_loci=True):

        self.output_dir = pathlib.Path(output_dir)
        self.threads = int(threads)
        self.species_list = list(map(str.title, species_list))
        self.download_loci=download_loci
        self.date_and_time = datetime.now().strftime('%d-%m-%Y %H:%M:%S')
        self.schemes = {'Salmonella': \
                            {'url': 'https://rest.pubmlst.org/db/pubmlst_salmonella_seqdef/schemes/4',
                            'scheme_summary_file': 'salmonella_scheme.json'},
                        'Escherichia': \
                            {'url': 'https://rest.pubmlst.org/db/pubmlst_escherichia_seqdef/schemes/6',
                            'scheme_summary_file': 'escherichia_scheme.json'},
                        'Test': \
                            {'url': 'https://rest.pubmlst.org/db/pubmlst_salmonella_seqdef/schemes/2',
                            'scheme_summary_file': 'test_scheme.json'}}
        species_to_download = [species for species in self.species_list if species in self.schemes]
        species_not_found = [species for species in self.species_list if species not in self.schemes]
        
        some_species_not_supported = len(species_not_found) > 0
        if some_species_not_supported:
            warning_message = f'The following species are currently not supported for cgMLST: {species_not_found}\n'\
                'The schemes for those species/genera will not be downloaded'
            warnings.warn(self.error_formatter(warning_message))
        no_species_supported = len(species_to_download) == 0
        if no_species_supported:
            error_message = 'None of the provided species is supported for cgMLST' \
                f'The supported genera/species are: {self.schemes.keys()}'
            raise ValueError(self.error_formatter(error_message))

        downloaded_schemes = {species:info for (species, info) in \
                self.download_pubmlst_scheme(species_list=species_to_download)}
        self.schemes = {species: {**self.schemes[species], **downloaded_schemes[species]} \
                for species in downloaded_schemes}
        print(yaml.dump(self.schemes, default_flow_style=False))

    def download_pubmlst_locus(self, locus_url, output_dir_species):
        output_file = locus_url.split('/')[-1] + '.fasta'
        locus_url = locus_url + '/alleles_fasta'
        if self.download_loci:
            subprocess.run(['wget', '--quiet', '--output-document', output_file, locus_url], 
                            check = True, timeout=200, 
                            cwd=output_dir_species)
        return True

    def download_pubmlst_scheme(self, species_list):
        for species in species_list:
            output_dir_per_species = self.output_dir.joinpath(species.lower())
            os.makedirs(output_dir_per_species, exist_ok=True)
            scheme_url = self.schemes[species]['url']
            output_file = output_dir_per_species.\
                joinpath(self.schemes[species]['scheme_summary_file'])
            subprocess.run(['wget', '--quiet', '--output-document', output_file, scheme_url], 
                            check = True, timeout=60)
            with open(output_file) as scheme_definition:
                scheme = json.load(scheme_definition)
            loci = db.from_sequence(scheme['loci'], 
                                    npartitions = self.threads)
            loci.map(self.download_pubmlst_locus, 
                    output_dir_species=output_dir_per_species).compute()
            scheme_info = {'scheme_description': scheme['description'],
                            'locus_count': scheme['locus_count']}
            yield species, scheme_info


if __name__ == '__main__':
    argument_parser = argparse.ArgumentParser(description='Download cgMLST schemes.')
    argument_parser.add_argument('-s','--species', metavar='N', 
                        type=str, nargs='+',
                        default=['salmonella', 'escherichia'],
                        help='Species whose cgMLST scheme should be downloaded')
    argument_parser.add_argument('-o', '--output_dir', type=pathlib.Path, 
                        default='output',
                        help='Output directory.')
    argument_parser.add_argument('-t', '--threads', type=int, 
                        default=4,
                        help='Number of threads to be used.')
    argument_parser.add_argument('--no_download', dest='download_loci', 
                                action='store_false')
    args = argument_parser.parse_args()
    cgMLSTSchemes(threads=args.threads, species_list=args.species, 
                    download_loci=args.download_loci,
                    output_dir=args.output_dir)