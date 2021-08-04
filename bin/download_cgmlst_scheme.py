import argparse
import base_juno_pipeline.helper_functions
import bs4
import dask.bag as db
from datetime import datetime
import json
import os
import pathlib
import requests
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
        self.species_list = list(map(str.lower, species_list))
        self.download_loci=download_loci
        self.date_and_time = datetime.now().strftime('%d-%m-%Y %H:%M:%S')
        self.schemes = {'salmonella': \
                            {'source': 'pubmlst',
                            'url': 'https://rest.pubmlst.org/db/pubmlst_salmonella_seqdef/schemes/4',
                            'scheme_summary_file': 'salmonella_scheme.json'},
                        'escherichia': \
                            {'source': 'pubmlst',
                            'url': 'https://rest.pubmlst.org/db/pubmlst_escherichia_seqdef/schemes/6',
                            'scheme_summary_file': 'escherichia_scheme.json'},
                        'yersinia': \
                            {'source': 'enterobase',
                            'url': 'http://enterobase.warwick.ac.uk/schemes/Yersinia.cgMLSTv1/',
                            'scheme_summary_file': 'yersinia_scheme.txt'},
                        # 'vibrio': \
                        #     {'source': 'enterobase',
                        #     'url': 'http://enterobase.warwick.ac.uk/schemes/Vibrio.cgMLSTv1/',
                        #     'scheme_summary_file': 'vibrio_scheme.txt'},
                        'test_pubmlst': \
                            {'source': 'pubmlst',
                            'url': 'https://rest.pubmlst.org/db/pubmlst_salmonella_seqdef/schemes/2',
                            'scheme_summary_file': 'test_scheme.json'},
                        'test_enterobase': \
                            {'source': 'enterobase',
                            'url': 'https://enterobase.warwick.ac.uk/schemes/Yersinia.Achtman7GeneMLST/',
                            'scheme_summary_file': 'test_scheme.json'}}
        species_to_download = [species for species in self.species_list if species in self.schemes]
        species_not_found = [species for species in self.species_list if species not in self.schemes]
        
        some_species_not_supported = len(species_not_found) > 0
        if some_species_not_supported:
            warning_message = f'The following species are currently not supported for cgMLST: {species_not_found}. '\
                'The schemes for those species/genera will not be downloaded'
            warnings.warn(self.error_formatter(warning_message))
        no_species_supported = len(species_to_download) == 0
        if no_species_supported:
            error_message = 'None of the provided species is supported for cgMLST. ' \
                f'The supported genera/species are: {self.schemes.keys()}'
            raise ValueError(self.error_formatter(error_message))

        downloaded_schemes = {species:info for (species, info) in \
                self.download_cgmlst_scheme(species_list=species_to_download)}
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

    def download_enterobase_locus(self, locus_url, output_dir_species):
        output_file = locus_url.split('/')[-1]
        if self.download_loci:
            subprocess.run(['wget', '--quiet', '--output-document', output_file, locus_url], 
                            check = True, timeout=200, 
                            cwd=output_dir_species)
            subprocess.run(['gunzip', output_file], 
                            check = True, timeout=600, 
                            cwd=output_dir_species)
        return True

    def download_pubmlst_scheme(self, scheme_url, output_dir_per_species, species):
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
        return scheme_info

    # Adapted from: https://github.com/kristyhoran/coreuscan/blob/master/coreuscan/coreuscan.py
    def download_enterobase_scheme(self, scheme_url, output_dir_per_species):
        website_data = requests.get(scheme_url)
        website_data.raise_for_status()
        parsed_website_data = bs4.BeautifulSoup(website_data.text, 'html.parser')
        loci_list = []
        for line in parsed_website_data.find_all('a'):
            if line.contents[0] != '../':
                sp = str(line.contents[0])
                sp = sp.strip('[]')
                if sp.endswith('fasta.gz'):
                    locus_url = scheme_url + sp
                    loci_list.append(locus_url)
        loci = db.from_sequence(loci_list, 
                                npartitions = self.threads)
        loci.map(self.download_enterobase_locus, 
                output_dir_species=output_dir_per_species).compute()
        scheme_info = {'scheme_description': None,
                        'locus_count': len(loci_list)}
        return scheme_info

    def download_cgmlst_scheme(self, species_list):
        for species in species_list:
            output_dir_per_species = self.output_dir.joinpath(species)
            os.makedirs(output_dir_per_species, exist_ok=True)
            scheme_url = self.schemes[species]['url']
            if self.schemes[species]['source'] == 'pubmlst':
                species_scheme_info = self.download_pubmlst_scheme(scheme_url, output_dir_per_species, species)
            elif self.schemes[species]['source'] == 'enterobase':
                species_scheme_info = self.download_enterobase_scheme(scheme_url, output_dir_per_species)
            yield species, species_scheme_info



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