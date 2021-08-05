import argparse
import base_juno_pipeline.helper_functions
from os import system
import pathlib
import subprocess
from yaml import safe_load


class runChewBBACA(base_juno_pipeline.helper_functions.JunoHelpers):
    """Class to produce and enlist (in a dictionary) the parameters necessary
    to run chewBBACA for different genera. The samples can be run also by calling
    the runChewBBACA method
    """

    def __init__(self, 
                best_hit_kmerfinder_files,
                sample_sheet='config/sample_sheet.yaml',
                output_dir='output/cgmlst/',
                cgmlst_db_dir='/mnt/db/juno/typing_db/cgmlst/prepared_schemes',
                threads=8):

        self.best_hit_kmerfinder_files = best_hit_kmerfinder_files
        self.sample_sheet = pathlib.Path(sample_sheet)
        assert self.sample_sheet.is_file(), f'The provided sample sheet {str(self.sample_sheet)} does not exist.'
        self.output_dir = pathlib.Path(output_dir)
        assert self.output_dir.is_dir(), f'The provided output directory {str(self.output_dir)} does not exist.'
        self.cgmlst_db_dir = pathlib.Path(cgmlst_db_dir)
        assert self.cgmlst_db_dir.is_dir(), f'The provided directory containing the cgMLST database {str(self.cgmlst_db_dir)} does not exist.'
        self.threads = int(threads)

        with open(self.sample_sheet) as sample_sheet_file:
            self.sample_sheet_dict = safe_load(sample_sheet_file)
        
        self.genera_set = self.get_genera_set()
        self.chewbbaca_parameters = self.get_chewbbaca_parameters_per_genus()

    def read_best_hit_kmerfinder_file(self, file_path):
        sample_name = pathlib.Path(file_path).parent.name
        with open(file_path) as file_:
            genus = file_.readlines()[0].split(' ')[0]
        return sample_name, genus

    def get_genus_per_sample(self):
        for genus in self.best_hit_kmerfinder_files:
            yield self.read_best_hit_kmerfinder_file(genus)

    def get_genera_set(self):
        samples_dict = {sample:genus for sample, genus in self.get_genus_per_sample()}
        genera_set = set(samples_dict.values())
        return genera_set
        
    def enlist_samples_per_genus(self):
        samples_dict = {sample:genus for sample, genus in self.get_genus_per_sample()}
        samples_per_genus_dict = {}
        for genus in self.genera_set:
            samples_per_genus_dict[genus] = {'samples': [sample for sample in samples_dict if samples_dict[sample] == genus]}
        return samples_per_genus_dict

    def make_file_with_samples_per_genus(self):
        samples_per_genus_dict = self.enlist_samples_per_genus()
        for genus in samples_per_genus_dict:
            genus_file = self.output_dir.joinpath(genus)
            with open(genus_file, 'w') as file_:
                for sample in samples_per_genus_dict[genus]['samples']:
                    assembly_file = self.sample_sheet_dict[sample]['assembly']
                    file_.write(assembly_file)
            samples_per_genus_dict[genus]['genus_file'] = str(genus_file)
        return samples_per_genus_dict

    def get_chewbbaca_parameters_per_genus(self):
        chewbbaca_parameters = self.make_file_with_samples_per_genus()
        for genus in chewbbaca_parameters:
            genus_cgmlst_db = self.cgmlst_db_dir.joinpath(genus)
            if genus_cgmlst_db.is_dir():
                chewbbaca_parameters[genus]['cgmlst_scheme'] = str(genus_cgmlst_db)
            else:
                print(self.error_formatter(f'There is no recognizable cgMLST scheme for the genus {genus.title()} in the {self.cgmlst_db_dir} directory.'))
                chewbbaca_parameters[genus]['cgmlst_scheme'] = None
        return chewbbaca_parameters

    def run_chewbbaca(self):
        for genus in self.chewbbaca_parameters:
            if self.chewbbaca_parameters[genus]['cgmlst_scheme'] is not None:
                print(self.message_formatter(f'Running samples for genus: {genus}'))
                chewbbaca_command = ['chewBBACA.py', 'AlleleCall', 
                                    '--cpu',  str(self.threads), 
                                    '--input-files', 
                                    self.chewbbaca_parameters[genus]['genus_file'], 
                                    '--schema-directory', 
                                    self.chewbbaca_parameters[genus]['cgmlst_scheme'],
                                    '--output-directory', 
                                    str(self.output_dir), 
                                    '--force-reset']
                print(f'Command used to run chewbbaca: \n{" ".join(chewbbaca_command)}')
                genus_cgmlst_output_dir = self.output_dir.joinpath(genus + '_samples')
                if not genus_cgmlst_output_dir.joinpath('results_alleles.tsv').exists():
                    output_dir_content_before_chewbbaca = [subdir for subdir in self.output_dir.iterdir()]
                    subprocess.run(chewbbaca_command,
                                    check = True, timeout=3600)
                    output_dir_content_after_chewbbaca = [subdir for subdir in self.output_dir.iterdir() if subdir not in output_dir_content_before_chewbbaca]
                    result_chewbbaca_dir = [subdir for subdir in output_dir_content_after_chewbbaca if subdir not in output_dir_content_before_chewbbaca]
                    assert len(result_chewbbaca_dir) == 1, f'Something went wrong while detecting the output directory of chewbbaca for genus {genus}'
                    if genus_cgmlst_output_dir.exists():
                        system(f'rm -rf {str(genus_cgmlst_output_dir)}')
                    result_chewbbaca_dir[0].rename(genus_cgmlst_output_dir)
            pathlib.Path(self.chewbbaca_parameters[genus]['genus_file']).unlink()


if __name__ == '__main__':
    argument_parser = argparse.ArgumentParser(description='run chewBBACA using the appropriate cgMLST scheme according to genus of sample(s)')
    argument_parser.add_argument('-f','--best-hit-kmerfinder-files', 
                                metavar='N', type=str, nargs='+',
                                help='Files (one per sample) containing the best hit for kmerfinder. They must contain only two words: one for the genus and one for the species of the sample')
    argument_parser.add_argument('-s','--sample-sheet', type=pathlib.Path, 
                                default='output/cgmlst',
                                help='Sample sheet (yaml file) containing a dictionary with the samples and at least one key:value pair assembly:file_path_to_assembly_file.\
                                    Example {sample1: {assembly: input_dir/sample1.fasta}}')
    argument_parser.add_argument('-o', '--output-dir', type=pathlib.Path, 
                                default='output/cgmlst',
                                help='Output directory for the chewBBACA results. A subfolder per genus will be created inside the output directory.')
    argument_parser.add_argument('-d', '--cgmlst-db-dir', type=pathlib.Path, 
                                default='/mnt/db/juno/typing_db/cgmlst/prepared_schemes',
                                help='Output directory for the chewBBACA results. A subfolder per genus will be created inside the output directory.')
    argument_parser.add_argument('-t', '--threads', type=int, 
                                default=8,
                                help='Number of threads to be used.')
    args = argument_parser.parse_args()
    chewbbaca_run = runChewBBACA(best_hit_kmerfinder_files=args.best_hit_kmerfinder_files)
    chewbbaca_run.run_chewbbaca()
