import argparse
from socket import timeout
import base_juno_pipeline.helper_functions as juno_helpers
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

cgmlst_schemes = {
    'salmonella': {
        'source': 'enterobase',
        'url': 'http://enterobase.warwick.ac.uk/schemes/Salmonella.cgMLSTv2/'
    },
    'campylobacter': {
        'source': 'pubmlst',
        'url': 'https://rest.pubmlst.org/db/pubmlst_campylobacter_seqdef/schemes/4'
    },
    'shigella': {
        'source': 'enterobase',
        'url': 'http://enterobase.warwick.ac.uk/schemes/Escherichia.cgMLSTv1/'
    },
    'escherichia': {
        'source': 'seqsphere',
        'url': 'https://www.cgmlst.org/ncs/schema/8896773/locus/'
    },
    'listeria': {
        'source': 'seqsphere',
        'url': 'https://www.cgmlst.org/ncs/schema/690488/locus/'
    },
    'listeria_optional': {
        'source': 'bigsdb_pasteur',
        'url': 'https://bigsdb.pasteur.fr/api/db/pubmlst_listeria_seqdef/schemes/3'
    },
    # Link not existing despite being used in official website
    # 'vibrio': \
    #     {'source': 'enterobase',
    #     'url': 'http://enterobase.warwick.ac.uk/schemes/VIBwgMLST.cgMLSTv1/'},
    'yersinia': {
        'source': 'enterobase',
        'url': 'http://enterobase.warwick.ac.uk/schemes/Yersinia.cgMLSTv1/'
    },
    'test_pubmlst': {
        'source': 'pubmlst',
        'url': 'https://rest.pubmlst.org/db/pubmlst_salmonella_seqdef/schemes/2'
    },
    'test_enterobase': {
        'source': 'enterobase',
        'url': 'https://enterobase.warwick.ac.uk/schemes/Yersinia.Achtman7GeneMLST/'
    },
    'test_bigsdb_pasteur': {
        'source': 'bigsdb_pasteur',
        'url': 'https://bigsdb.pasteur.fr/api/db/pubmlst_listeria_seqdef/schemes/2'
    }
}



class cgMLSTSchemes(juno_helpers.JunoHelpers):
    """Class containing information and functions to download cgMLST schemes"""

    def __init__(self, 
                genus_list,
                output_dir='output',
                threads=2,
                download_loci=True):

        self.output_dir = pathlib.Path(output_dir)
        self.threads = int(threads)
        self.genus_list = [genus.lower() for genus in genus_list]
        self.download_loci=download_loci
        self.date_and_time = datetime.now().strftime('%d-%m-%Y %H:%M:%S')
        self.schemes =  cgmlst_schemes
        genus_to_download = [genus for genus in self.genus_list if genus in self.schemes]
        # Some genera require more than one scheme (normal + optional)
        # so both need to be downloaded here
        if 'escherichia' in genus_to_download:
            genus_to_download.append('shigella')
        if 'shigella' in genus_to_download:
            genus_to_download.append('escherichia')
        if 'listeria' in genus_to_download:
            genus_to_download.append('listeria_optional')
        genus_to_download = list(set(genus_to_download))
        genus_not_found = [genus for genus in self.genus_list if genus not in self.schemes]
        
        some_genus_not_supported = len(genus_not_found) > 0
        if some_genus_not_supported:
            warnings.warn(
                self.error_formatter(
                    f'The following genus are currently not supported for cgMLST: {",".join(genus_not_found)}. '\
                    'The schemes for those genus/genera will not be downloaded'
                )
            )

        all_nonsupported = len(genus_to_download) == 0
        if all_nonsupported:
            supported_genera = [genus for genus in cgmlst_schemes.keys() if not genus.startswith('test_')]
            raise ValueError(
                self.error_formatter(
                    'None of the provided genera is supported for cgMLST. ' \
                    f'The supported genera/species are: {",".join(supported_genera)}'
                )
            )

        downloaded_schemes = {
            genus:info for (genus, info) in \
            self.download_cgmlst_scheme(genus_list=genus_to_download)
        }
        self.schemes = {
            genus: {**self.schemes[genus], **downloaded_schemes[genus]} \
            for genus in downloaded_schemes
        }
        print(yaml.dump(self.schemes, default_flow_style=False))

    def __download_wget(self, url, output_file, working_dir='.', timeout=200):
        subprocess.run(
            f'wget --quiet --output-document {output_file} {url}',
            shell=True, check=True, timeout=timeout,
            cwd=working_dir
        )
    def download_pubmlst_locus(self, locus_url, output_dir_genus):
        output_file = locus_url.split('/')[-1] + '.fasta'
        locus_url = locus_url + '/alleles_fasta'
        if self.download_loci:
            self.__download_wget(
                locus_url, output_file, working_dir=output_dir_genus
            )
        return True

    def download_enterobase_locus(self, locus_url, output_dir_genus):
        output_file = locus_url.split('/')[-1]
        self.__download_wget(
                locus_url, output_file, working_dir=output_dir_genus
            )
        subprocess.run(
            ['gunzip', output_file], 
            check = True, timeout=100, 
            cwd=output_dir_genus
        )
        return True

    def download_seqsphere_locus(self, locus_url, output_dir_genus):
        output_file = locus_url.split('/')[-1]
        if self.download_loci:
            self.__download_wget(
                locus_url, output_file, working_dir=output_dir_genus
            )
        return True

    def download_pubmlst_scheme(self, scheme_url, output_dir_per_genus, genus):
        output_file = output_dir_per_genus.\
                    joinpath('scheme_summary_file')
        # Need to make it again because it is only made before if download_loci=True
        os.makedirs(output_dir_per_genus, exist_ok=True)
        self.__download_wget(scheme_url, output_file, timeout=60)
        with open(output_file) as scheme_definition:
            scheme = json.load(scheme_definition)
        if self.download_loci:
            print(
                self.message_formatter(
                    f'Downloading {scheme["locus_count"]} loci from PubMLST server...'
                )
            )
            loci = db.from_sequence(
                scheme['loci'], npartitions = self.threads
            )
            loci.map(
                self.download_pubmlst_locus, 
                output_dir_genus=output_dir_per_genus
            ).compute()
        else:
            # Remove created files/directories if loci are not supposed to be downloaded
            os.unlink(output_file)
            try:
                # Do not remove the directories if they were not empty before
                # os.rmdir will fail if the directory is not empty
                os.rmdir(output_dir_per_genus)
                os.rmdir(self.output_dir)
            except:
                pass
        scheme_info = {
            'scheme_description': scheme['description'],
            'locus_count': scheme['locus_count']
        }
        return scheme_info

    # Adapted from: https://github.com/kristyhoran/coreuscan/blob/master/coreuscan/coreuscan.py
    def download_enterobase_scheme(self, scheme_url, output_dir_per_genus):
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
        if self.download_loci:
            print(
                self.message_formatter(
                    f'Downloading {len(loci_list)} loci from Enterobase server...'
                )
            )
            loci = db.from_sequence(loci_list, 
                                npartitions = self.threads)
            loci.map(self.download_enterobase_locus, 
                    output_dir_genus=output_dir_per_genus).compute()
        scheme_info = {
            'scheme_description': None,
            'locus_count': len(loci_list)
        }
        return scheme_info

    def download_seqsphere_scheme(self, scheme_url, output_dir_per_genus):
        website_data = requests.get(scheme_url)
        website_data.raise_for_status()
        parsed_website_data = bs4.BeautifulSoup(website_data.text, 'html.parser')
        loci_list = []
        for line in parsed_website_data.find_all('a'):
            if line.contents[0] != '../':
                sp = str(line.contents[0])
                sp = sp.strip('[]')
                if sp.endswith('.fasta'):
                    locus_url = scheme_url + sp
                    loci_list.append(locus_url)
        if self.download_loci:
            print(self.message_formatter(f'Downloading {len(loci_list)} loci from SeqSphere+ server...'))
            loci = db.from_sequence(loci_list, 
                                    npartitions = self.threads)
            loci.map(self.download_seqsphere_locus, 
                    output_dir_genus=output_dir_per_genus).compute()
        scheme_info = {
            'scheme_description': None,
            'locus_count': len(loci_list)
        }
        return scheme_info

    def download_cgmlst_scheme(self, genus_list):
        for genus in genus_list:
            source = self.schemes[genus]['source']
            print(
                self.message_formatter(
                    f'Collecting cgMLST scheme for {genus.title()} from {source.title()}...'
                )
            )
            output_dir_per_genus = self.output_dir.joinpath(genus)
            if self.download_loci:
                os.makedirs(output_dir_per_genus, exist_ok=True)
            scheme_url = self.schemes[genus]['url']
            if source == 'pubmlst':
                genus_scheme_info = self.download_pubmlst_scheme(scheme_url, output_dir_per_genus, genus)
            elif source == 'bigsdb_pasteur':
                # The BIGSDB Pasteur database works exactly the same than PubMLST
                genus_scheme_info = self.download_pubmlst_scheme(scheme_url, output_dir_per_genus, genus)
            elif source == 'enterobase':
                genus_scheme_info = self.download_enterobase_scheme(scheme_url, output_dir_per_genus)
            elif source == 'seqsphere':
                genus_scheme_info = self.download_seqsphere_scheme(scheme_url, output_dir_per_genus)
            yield genus, genus_scheme_info


def main():
    argument_parser = argparse.ArgumentParser(
        description='Download cgMLST schemes.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    choices_scheme = [genus for genus in cgmlst_schemes.keys() if not genus.startswith('test_')]
    argument_parser.add_argument(
        '-g','--genus', metavar='N', 
        type=str, nargs='+',
        choices=choices_scheme,
        default=[genus for genus in choices_scheme if not genus.endswith('_optional')],
        help='Genus/genera whose cgMLST scheme should be downloaded. The supported genera are: %(choices)s.'
    )
    argument_parser.add_argument(
        '-o', '--output_dir', type=pathlib.Path, 
        default='output',
        help='Output directory.'
    )
    argument_parser.add_argument(
        '-t', '--threads', type=int, 
        default=4,
        help='Number of threads to be used.'
    )
    argument_parser.add_argument(
        '--no_download', dest='download_loci', 
        action='store_false'
    )
    args = argument_parser.parse_args()
    cgMLSTSchemes(
        threads=args.threads, genus_list=args.genus, 
        download_loci=args.download_loci,
        output_dir=args.output_dir
    )


if __name__ == '__main__':
    main()