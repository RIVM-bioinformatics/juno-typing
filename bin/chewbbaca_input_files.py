import argparse
import base_juno_pipeline.helper_functions
from os import system
import pathlib
import subprocess
from yaml import safe_load


class inputChewBBACA(base_juno_pipeline.helper_functions.JunoHelpers):
    """Class to produce and enlist (in a dictionary) the parameters necessary
    to run chewBBACA for different genera. The samples can be run also by calling
    the runChewBBACA method
    """

    def __init__(self, 
                best_hit_kmerfinder_files,
                sample_sheet='config/sample_sheet.yaml',
                output_dir='output/cgmlst/'):

        self.supported_genera = ['campylobacter', 'escherichia', 'listeria', 
                                'listeria_optional', 'salmonella','shigella', 
                                'yersinia']
        self.best_hit_kmerfinder_files = best_hit_kmerfinder_files
        self.sample_sheet = pathlib.Path(sample_sheet)
        assert self.sample_sheet.is_file(), f'The provided sample sheet {str(self.sample_sheet)} does not exist.'
        self.output_dir = pathlib.Path(output_dir)


    def __read_sample_sheet(self):
        print("Reading sample sheet...\n")
        with open(self.sample_sheet) as sample_sheet_file:
            self.sample_sheet_dict = safe_load(sample_sheet_file)

    def __read_best_hit_kmerfinder_file(self, file_path):
        sample_name = pathlib.Path(file_path).parent.name
        with open(file_path) as file_:
            genus = [file_.readlines()[0].split(' ')[0]]
        if genus[0] == 'listeria':
            genus.append('listeria_optional')
        elif genus[0] == 'escherichia':
            genus.append('shigella')
        elif genus[0] == 'shigella':
            genus.append('escherichia')
        return sample_name, genus

    def __get_genus_per_sample(self):
        for genus in self.best_hit_kmerfinder_files:
            print(f'Getting genus {genus}...\n')
            yield self.__read_best_hit_kmerfinder_file(genus)

    def __get_samples_dict(self):
        print(f'Getting list of all genera included in this sample set...\n')
        self.samples_dict = {sample:genus for sample, genus in self.__get_genus_per_sample()}

    def __get_genera_set(self):
        print(f'Getting list of all genera included in this sample set...\n')
        self.__get_samples_dict()
        genera_set = set([genus for genus_list in self.samples_dict.values() for genus in genus_list if genus in self.supported_genera])
        self.genera_set = genera_set

    def __enlist_samples_per_genus(self):
        print(f'Getting list of samples per genus found...\n')
        self.__get_genera_set()
        genera_input_dict = {}
        for genus in self.genera_set:
            samples_in_genus = [sample for sample in self.samples_dict if genus in self.samples_dict[sample]]
            genera_input_dict[genus] = {'samples': samples_in_genus}
        print(f'genera input: {genera_input_dict}')
        return genera_input_dict

    def make_file_with_samples_per_genus(self):
        self.__read_sample_sheet()
        genera_input_dict = self.__enlist_samples_per_genus()
        print(f'genera input: {genera_input_dict}')
        for genus in genera_input_dict:
            genus_file = self.output_dir.joinpath(genus + '_samples.txt')
            with open(genus_file, 'w') as file_:
                for sample in genera_input_dict[genus]['samples']:
                    assembly_file = self.sample_sheet_dict[sample]['assembly']
                    file_.write(assembly_file+'\n')
            genera_input_dict[genus]['genus_file'] = str(genus_file)
        print(f'Files with samples per genus will be written in {self.output_dir}!\n')
        self.genera_input_dict = genera_input_dict
        return genera_input_dict


if __name__ == '__main__':
    argument_parser = argparse.ArgumentParser(description='run chewBBACA using the appropriate cgMLST scheme according to genus of sample(s)')
    argument_parser.add_argument('-f','--best-hit-kmerfinder-files', 
                                metavar='N', type=str, nargs='+',
                                help='Files (one per sample) containing the best hit for kmerfinder. They must contain only two words: one for the genus and one for the species of the sample')
    argument_parser.add_argument('-s','--sample-sheet', type=pathlib.Path, 
                                default='config/sample_sheet.yaml',
                                help='Sample sheet (yaml file) containing a dictionary with the samples and at least one key:value pair assembly:file_path_to_assembly_file.\
                                    Example {sample1: {assembly: input_dir/sample1.fasta}}')
    argument_parser.add_argument('-o', '--output-dir', type=pathlib.Path, 
                                default='output/cgmlst',
                                help='Output directory for the chewBBACA results. A subfolder per genus will be created inside the output directory.')
    args = argument_parser.parse_args()
    chewbbaca_run = inputChewBBACA(best_hit_kmerfinder_files=args.best_hit_kmerfinder_files,
                                    sample_sheet=args.sample_sheet,
                                    output_dir=args.output_dir)
    chewbbaca_run.make_file_with_samples_per_genus()
