import os
import pathlib
from sys import path
import unittest
from pandas import read_csv

main_script_path = str(pathlib.Path(pathlib.Path(__file__).parent.absolute()).parent.absolute())
path.insert(0, main_script_path)
from bin import general_juno_pipeline
from bin import download_dbs
from bin import serotypefinder_multireport
import juno_typing


class TestPipelineStartup(unittest.TestCase):

    def setUpClass(): 
        fake_dirs = ['fake_dir_empty', 
                    'fake_dir_wsamples', 
                    'fake_dir_incomplete',
                    'fake_dir_juno', 
                    'fake_dir_juno/clean_fastq', 
                    'fake_dir_juno/de_novo_assembly_filtered']

        fake_files = ['fake_dir_wsamples/sample1_R1.fastq',
                    'fake_dir_wsamples/sample1_R2.fastq.gz',
                    'fake_dir_wsamples/sample2_R1_filt.fq',
                    'fake_dir_wsamples/sample2_R2_filt.fq.gz', 
                    'fake_dir_wsamples/sample1.fasta',
                    'fake_dir_wsamples/sample2.fasta',
                    'fake_dir_incomplete/sample1_R1.fastq',
                    'fake_dir_incomplete/sample1_R2.fastq.gz',
                    'fake_dir_incomplete/sample2_R1_filt.fq',
                    'fake_dir_incomplete/sample2_R2_filt.fq.gz',
                    'fake_dir_incomplete/sample2.fasta',
                    'fake_dir_juno/clean_fastq/1234_R1.fastq.gz',
                    'fake_dir_juno/clean_fastq/1234_R2.fastq.gz', 
                    'fake_dir_juno/de_novo_assembly_filtered/1234.fasta']     
                    
        for folder in fake_dirs:
            pathlib.Path(folder).mkdir(exist_ok = True)
        for file_ in fake_files:
            pathlib.Path(file_).touch(exist_ok = True)

    def tearDownClass():
        fake_dirs = ['fake_dir_empty', 
                    'fake_dir_wsamples', 
                    'fake_dir_incomplete',
                    'fake_dir_juno', 
                    'fake_dir_juno/clean_fastq', 
                    'fake_dir_juno/de_novo_assembly_filtered']

        for folder in fake_dirs:
            os.system('rm -rf {}'.format(str(folder)))

    def test_emptydir(self):
        self.assertRaises(ValueError, general_juno_pipeline.PipelineStartup, pathlib.Path('fake_dir_empty'), 'both')

    def test_incompletedir(self):
        self.assertRaises(KeyError, general_juno_pipeline.PipelineStartup, pathlib.Path('fake_dir_incomplete'), 'both')

    def test_correctdir_wdifffastqextensions(self):
        expected_output = {'sample1': {'R1': 'fake_dir_wsamples/sample1_R1.fastq', 
                                        'R2': 'fake_dir_wsamples/sample1_R2.fastq.gz', 
                                        'assembly': 'fake_dir_wsamples/sample1.fasta'}, 
                            'sample2': {'R1': 'fake_dir_wsamples/sample2_R1_filt.fq', 
                                        'R2': 'fake_dir_wsamples/sample2_R2_filt.fq.gz', 
                                        'assembly': 'fake_dir_wsamples/sample2.fasta'}}
        pipeline = general_juno_pipeline.PipelineStartup(pathlib.Path('fake_dir_wsamples'), 'both')
        self.assertDictEqual(pipeline.sample_dict, expected_output)

    def test_junodir_wnumericsamplenames(self):
        expected_output = {'1234': {'R1': 'fake_dir_juno/clean_fastq/1234_R1.fastq.gz', 
                                        'R2': 'fake_dir_juno/clean_fastq/1234_R2.fastq.gz', 
                                        'assembly': 'fake_dir_juno/de_novo_assembly_filtered/1234.fasta'}}
                
        pipeline = general_juno_pipeline.PipelineStartup(pathlib.Path('fake_dir_juno'), 'both')
        self.assertDictEqual(pipeline.sample_dict, expected_output)

    def test_string_accepted_as_inputdir(self):
        expected_output = {'1234': {'R1': 'fake_dir_juno/clean_fastq/1234_R1.fastq.gz', 
                                        'R2': 'fake_dir_juno/clean_fastq/1234_R2.fastq.gz', 
                                        'assembly': 'fake_dir_juno/de_novo_assembly_filtered/1234.fasta'}}
                
        pipeline = general_juno_pipeline.PipelineStartup('fake_dir_juno', 'both')
        self.assertDictEqual(pipeline.sample_dict, expected_output)


# class TestDownloadDbs(unittest.TestCase):

#     def setUpClass():
#         pathlib.Path('fake_db').mkdir(exist_ok = True)

#     def tearDownClass():
#         os.system('rm -rf fake_db')

#     def test_download_kmerfinder_software(self):
#         path_to_db = pathlib.Path('fake_db', 'kmerfinder_soft')
#         download_dbs.download_software_kmerfinder(path_to_db)
#         self.assertTrue(path_to_db.joinpath('kmerfinder.py').exists())

#     def test_download_mlst7_software(self):
#         path_to_db = pathlib.Path('fake_db', 'mlst7_soft')
#         download_dbs.download_software_mlst7(path_to_db)
#         self.assertTrue(path_to_db.joinpath('mlst.py').exists())

#     def test_download_kmerfinder_db(self):
#         path_to_db = pathlib.Path('fake_db', 'kmerfinder_db')
#         download_dbs.download_db_kmerfinder(path_to_db)
#         self.assertTrue(path_to_db.joinpath('config').exists())

#     def test_download_mlst7_db(self):
#         path_to_db = pathlib.Path('fake_db', 'mlst7_db')
#         download_dbs.download_db_mlst7(path_to_db)
#         self.assertTrue(path_to_db.joinpath('senterica', 'senterica.length.b').exists())

#     def test_download_serotypefinder_db(self):
#         path_to_db = pathlib.Path('fake_db', 'serotypefinder_db')
#         download_dbs.download_db_serotypefinder(path_to_db)
#         self.assertTrue(path_to_db.joinpath('H_type.seq.b').exists())

#     def test_download_seroba_db(self):
#         path_to_db = pathlib.Path('fake_db', 'seroba_db')
#         download_dbs.download_db_seroba(path_to_db, kmersize = 71)
#         self.assertTrue(path_to_db.joinpath('database', 'cdhit_cluster').exists())

#     def test_download_all_dbs(self):
#         path_to_db = pathlib.Path('fake_db')
#         # bin_path = pathlib.Path(__file__).parent.absolute().joinpath('../bin')
#         download_dbs.get_downloads_juno_typing(path_to_db, path_to_db, False, 71)
#         self.assertTrue(path_to_db.joinpath('kmerfinder', 'kmerfinder.py').exists())
#         self.assertTrue(path_to_db.joinpath('cge-mlst', 'mlst.py').exists())
#         self.assertTrue(path_to_db.joinpath('kmerfinder_db', 'config').exists())
#         self.assertTrue(path_to_db.joinpath('mlst7_db', 'senterica', 'senterica.length.b').exists())
#         self.assertTrue(path_to_db.joinpath('serotypefinder_db', 'H_type.seq.b').exists())
#         self.assertTrue(path_to_db.joinpath('seroba_db', 'database', 'cdhit_cluster').exists())


class TestJunoTypingPipeline(unittest.TestCase):
    def setUpClass():
        fake_dirs = ['fake_dir_wsamples', 
                    'fake_dir_juno', 
                    'fake_dir_juno/clean_fastq', 
                    'fake_dir_juno/de_novo_assembly_filtered',
                    'fake_db',
                    'fake_db/kmerfinder_db',
                    'fake_db/kmerfinder_db/bacteria',
                    'fake_db/mlst7_db/',
                    'fake_db/mlst7_db/senterica/']

        fake_files = ['fake_dir_wsamples/sample1_R1.fastq',
                    'fake_dir_wsamples/sample1_R2.fastq.gz',
                    'fake_dir_wsamples/sample2_R1_filt.fq',
                    'fake_dir_wsamples/sample2_R2_filt.fq.gz', 
                    'fake_dir_wsamples/sample1.fasta',
                    'fake_dir_wsamples/sample2.fasta',
                    'fake_dir_juno/clean_fastq/1234_R1.fastq.gz',
                    'fake_dir_juno/clean_fastq/1234_R2.fastq.gz', 
                    'fake_dir_juno/de_novo_assembly_filtered/1234.fasta',
                    'fake_db/kmerfinder_db/bacteria/bacteria.ATG.length.b',
                    'fake_db/mlst7_db/senterica/senterica.length.b']
                            
        for folder in fake_dirs:
            pathlib.Path(folder).mkdir(exist_ok = True)
        for file_ in fake_files:
            pathlib.Path(file_).touch(exist_ok = True)

    def tearDownClass():
        fake_dirs = ['fake_dir_wsamples', 
                    'fake_dir_juno', 
                    'fake_db']
        for folder in fake_dirs:
            os.system('rm -rf {}'.format(str(folder)))
    
    def test_junotyping_dryrun(self):
        raised = False
        try:
            juno_typing.JunoTypingRun(input_dir = 'fake_dir_wsamples', 
                                    output_dir = pathlib.Path('output'), 
                                    db_dir = pathlib.Path('fake_db'),
                                    dryrun = True)
        except InputError:
            raised = True
        self.assertFalse(raised, 'Exception raised when running a dryrun')

class TestSerotypeFinderMultireport(unittest.TestCase):
    def setUpClass():
        assert pathlib.Path('tests/example_input/fake_result_serotype1.csv').exists, "Missing input file for TestSerotypeFinderMultireport test"
    
    def test_findallele(self):
        result_csv = read_csv('tests/example_input/fake_result_serotype1.csv', index_col = 0)
        allele = serotypefinder_multireport.find_allele(result_csv, 'wzx')
        self.assertEqual(allele, 'O50/O2')

    def test_getsampleserotype(self):
        alleles = serotypefinder_multireport.get_sample_serotype('tests/example_input/fake_result_serotype1.csv')
        for allele in ['O50/O2', 'O2','H6']:
            self.assertTrue(allele in alleles, 'The result_serotype.csv obtained by SerotypeFinder is not being properly processed')

    def test_mergemultiplereports(self):
        list_files = ['tests/example_input/fake_result_serotype1.csv',
                    'tests/example_input/fake_result_serotype2.csv',
                    'tests/example_input/fake_result_serotype3.csv']
        sample_names = ['sample1', 'sample2', 'sample3']
        multireport = serotypefinder_multireport.merge_multiple_reports(list_files, sample_names)
        self.assertEqual(multireport.loc['sample1','O type'].iat[0], 'O2/O50')
        self.assertEqual(multireport.loc['sample2','O type'].iat[0], 'O128')
        self.assertEqual(multireport.loc['sample3','O type'].iat[0], 'Error! No O type found')

if __name__ == '__main__':
    unittest.main(exit=False)
    os.system('rm -rf fake_d*')