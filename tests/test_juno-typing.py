import csv
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
    """Testing the pipeline startup (generating dict with samples) from general Juno pipelines"""
    
    def setUpClass(): 
        """Making fake directories and files to test different case scenarios for starting pipeline"""

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
        """Removing fake directories/files"""

        fake_dirs = ['fake_dir_empty', 
                    'fake_dir_wsamples', 
                    'fake_dir_incomplete',
                    'fake_dir_juno', 
                    'fake_dir_juno/clean_fastq', 
                    'fake_dir_juno/de_novo_assembly_filtered']

        for folder in fake_dirs:
            os.system('rm -rf {}'.format(str(folder)))

    def test_emptydir(self):
        """Testing the pipeline startup fails if the input directory does not have expected files"""
        self.assertRaises(ValueError, general_juno_pipeline.PipelineStartup, pathlib.Path('fake_dir_empty'), 'both')

    def test_incompletedir(self):
        """Testing the pipeline startup fails if the input directory is missing some of the fasta files for the fastq files"""
        self.assertRaises(KeyError, general_juno_pipeline.PipelineStartup, pathlib.Path('fake_dir_incomplete'), 'both')

    def test_correctdir_wdifffastqextensions(self):
        """Testing the pipeline startup accepts fastq and fastq.gz files"""

        expected_output = {'sample1': {'R1': 'fake_dir_wsamples/sample1_R1.fastq', 
                                        'R2': 'fake_dir_wsamples/sample1_R2.fastq.gz', 
                                        'assembly': 'fake_dir_wsamples/sample1.fasta'}, 
                            'sample2': {'R1': 'fake_dir_wsamples/sample2_R1_filt.fq', 
                                        'R2': 'fake_dir_wsamples/sample2_R2_filt.fq.gz', 
                                        'assembly': 'fake_dir_wsamples/sample2.fasta'}}
        pipeline = general_juno_pipeline.PipelineStartup(pathlib.Path('fake_dir_wsamples'), 'both')
        self.assertDictEqual(pipeline.sample_dict, expected_output)

    def test_junodir_wnumericsamplenames(self):
        """Testing the pipeline startup converts numeric file names to string"""

        expected_output = {'1234': {'R1': 'fake_dir_juno/clean_fastq/1234_R1.fastq.gz', 
                                        'R2': 'fake_dir_juno/clean_fastq/1234_R2.fastq.gz', 
                                        'assembly': 'fake_dir_juno/de_novo_assembly_filtered/1234.fasta'}}
                
        pipeline = general_juno_pipeline.PipelineStartup(pathlib.Path('fake_dir_juno'), 'both')
        self.assertDictEqual(pipeline.sample_dict, expected_output)

    def test_string_accepted_as_inputdir(self):
        """Testing the pipeline startup accepts string (not only pathlib.Path) as input"""

        expected_output = {'1234': {'R1': 'fake_dir_juno/clean_fastq/1234_R1.fastq.gz', 
                                        'R2': 'fake_dir_juno/clean_fastq/1234_R2.fastq.gz', 
                                        'assembly': 'fake_dir_juno/de_novo_assembly_filtered/1234.fasta'}}
                
        pipeline = general_juno_pipeline.PipelineStartup('fake_dir_juno', 'both')
        self.assertDictEqual(pipeline.sample_dict, expected_output)



# @unittest.skipIf(not pathlib.Path('/mnt/scratch_dir/hernanda').exists(),
#                      "Skipped in non-RIVM environments (for sake of time)")
# class TestDownloadDbs(unittest.TestCase):
#     """Testing the downloading of databases and software used by the pipeline"""

#     def setUpClass():
#         pathlib.Path('fake_db').mkdir(exist_ok = True)

#     def tearDownClass():
#         os.system('rm -rf fake_db')

#     def test_download_kmerfinder_software(self):
#         """Testing kmerfinder is properly downloaded and set up"""

#         path_to_db = pathlib.Path('fake_db', 'kmerfinder_soft')
#         download_dbs.download_software_kmerfinder(path_to_db)
#         self.assertTrue(path_to_db.joinpath('kmerfinder.py').exists())

#     def test_download_mlst7_software(self):
#         """Testing cge-mlst is properly downloaded and set up"""

#         path_to_db = pathlib.Path('fake_db', 'mlst7_soft')
#         download_dbs.download_software_mlst7(path_to_db)
#         self.assertTrue(path_to_db.joinpath('mlst.py').exists())

#     def test_download_kmerfinder_db(self):
#         """Testing kmerfinder database is properly downloaded and set up"""

#         path_to_db = pathlib.Path('fake_db', 'kmerfinder_db')
#         download_dbs.download_db_kmerfinder(path_to_db)
#         self.assertTrue(path_to_db.joinpath('config').exists())

#     def test_download_mlst7_db(self):
#         """Testing cge-mlst database is properly downloaded and set up"""

#         path_to_db = pathlib.Path('fake_db', 'mlst7_db')
#         download_dbs.download_db_mlst7(path_to_db)
#         self.assertTrue(path_to_db.joinpath('senterica', 'senterica.length.b').exists())

#     def test_download_serotypefinder_db(self):
#         """Testing SerotypeFinder database is properly downloaded and set up"""

#         path_to_db = pathlib.Path('fake_db', 'serotypefinder_db')
#         download_dbs.download_db_serotypefinder(path_to_db)
#         self.assertTrue(path_to_db.joinpath('H_type.seq.b').exists())

#     def test_download_seroba_db(self):
#         """Testing seroba database is properly downloaded and set up"""

#         path_to_db = pathlib.Path('fake_db', 'seroba_db')
#         download_dbs.download_db_seroba(path_to_db, kmersize = 71)
#         self.assertTrue(path_to_db.joinpath('database', 'cdhit_cluster').exists())

#     def test_download_all_dbs(self):
#         """Testing all databases/software are downloaded and set up"""

#         path_to_db = pathlib.Path('fake_db')
#         download_dbs.get_downloads_juno_typing(path_to_db, path_to_db, False, 71)
#         self.assertTrue(path_to_db.joinpath('kmerfinder', 'kmerfinder.py').exists())
#         self.assertTrue(path_to_db.joinpath('cge-mlst', 'mlst.py').exists())
#         self.assertTrue(path_to_db.joinpath('kmerfinder_db', 'config').exists())
#         self.assertTrue(path_to_db.joinpath('mlst7_db', 'senterica', 'senterica.length.b').exists())
#         self.assertTrue(path_to_db.joinpath('serotypefinder_db', 'H_type.seq.b').exists())
#         self.assertTrue(path_to_db.joinpath('seroba_db', 'database', 'cdhit_cluster').exists())



class TestJunoTypingDryRun(unittest.TestCase):
    """Testing the JunoTyping class (code specific for this pipeline)"""

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

        with open('fake_dir_wsamples/fake_metadata.csv', mode='w') as metadata_file:
            metadata_writer = csv.writer(metadata_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

            metadata_writer.writerow(['Sample', 'Genus', 'Species'])
            metadata_writer.writerow(['sample1', 'Salmonella', 'enterica'])
            metadata_writer.writerow(['sample2', 'Escherichia', 'coli'])


    def tearDownClass():
        fake_dirs = ['fake_dir_wsamples', 
                    'fake_dir_juno', 
                    'fake_db']
        for folder in fake_dirs:
            os.system('rm -rf {}'.format(str(folder)))
    
    def test_junotyping_dryrun(self):
        """Testing the pipeline runs properly as a dry run"""
        raised = False
        try:
            juno_typing.JunoTypingRun(input_dir = 'fake_dir_wsamples', 
                                    metadata= None,
                                    output_dir = pathlib.Path('output'), 
                                    db_dir = pathlib.Path('fake_db'),
                                    dryrun = True)
        except:
            raised = True
            raise
        self.assertFalse(raised, 'Exception raised when running a dryrun')

    def test_junotyping_dryrun_wMetadata(self):
        """Testing the pipeline runs properly as a dry run when providing a metadata file"""
        raised = False
        try:
            juno_typing.JunoTypingRun(input_dir = 'fake_dir_wsamples', 
                                    metadata= 'fake_dir_wsamples/fake_metadata.csv',
                                    output_dir = pathlib.Path('output'), 
                                    db_dir = pathlib.Path('fake_db'),
                                    dryrun = True)
        except:
            raised = True
            raise
        self.assertFalse(raised, 'Exception raised when running a dryrun and providing a metadata file')



@unittest.skipIf(not pathlib.Path('/data/BioGrid/hernanda/test_data_per_pipeline/Enteric/Juno-typing/').exists(),
                    "Skipped in non-RIVM environments (for sake of time)")
class TestJunoTypingPipeline(unittest.TestCase):
    """Testing the JunoTyping class (code specific for this pipeline)"""

    def setUpClass():
        os.system('rm -rf test_output')

    def tearDownClass():
        os.system('rm -rf test_output')
    
    def test_junotyping_run(self):
        """Testing the pipeline runs properly with real samples"""
        raised = False
        try:
            juno_typing.JunoTypingRun(input_dir = '/data/BioGrid/hernanda/test_data_per_pipeline/Enteric/Juno-typing/', 
                                    db_dir = '/mnt/db/juno/typing_db',
                                    output_dir = pathlib.Path('test_output'))
        except:
            raised = True
            raise
        self.assertFalse(raised, 'Exception raised when running a dryrun')

    def test_junotyping_run_wMetadata(self):
        """Testing the pipeline runs properly with real samples when providing a metadata file"""
        raised = False
        try:
            juno_typing.JunoTypingRun(input_dir = '/data/BioGrid/hernanda/test_data_per_pipeline/Enteric/Juno-typing/', 
                                    db_dir = '/mnt/db/juno/typing_db',
                                    metadata= '/data/BioGrid/hernanda/test_data_per_pipeline/Enteric/Juno-typing/metadata.csv',
                                    output_dir = pathlib.Path('test_output'))
        except:
            raised = True
            raise
        self.assertFalse(raised, 'Exception raised when running a dryrun')



class TestSerotypeFinderMultireport(unittest.TestCase):
    def setUpClass():
        assert pathlib.Path('tests/example_input/fake_result_serotype1.csv').exists, "Missing input file for TestSerotypeFinderMultireport test"
    
    def test_findallele(self):
        """The code should look for the allele result (wzx in this case) in the results of SerotypeFinder (E. coli serotyper)"""

        result_csv = read_csv('tests/example_input/fake_result_serotype1.csv', index_col = 0)
        allele = serotypefinder_multireport.find_allele(result_csv, 'wzx')
        self.assertEqual(allele, 'O50/O2')

    def test_getsampleserotype(self):
        """The code should return all the possible serotypes for every locus from the results of SerotypeFinder (E. coli serotyper)"""

        alleles = serotypefinder_multireport.get_sample_serotype('tests/example_input/fake_result_serotype1.csv')
        for allele in ['O50/O2', 'O2','H6']:
            self.assertTrue(allele in alleles, 'The result_serotype.csv obtained by SerotypeFinder is not being properly processed')

    def test_mergemultiplereports(self):
        """Testing properly reporting of O-type. Especial attention to multiple alleles arranged from smaller to larger from SerotypeFinder (E. coli serotyper)"""

        list_files = ['tests/example_input/fake_result_serotype1.csv',
                    'tests/example_input/fake_result_serotype2.csv',
                    'tests/example_input/fake_result_serotype3.csv']
        sample_names = ['sample1', 'sample2', 'sample3']
        multireport = serotypefinder_multireport.merge_multiple_reports(list_files, sample_names)
        self.assertEqual(multireport.loc['sample1','O type'].iat[0], 'O2/O50')
        self.assertEqual(multireport.loc['sample2','O type'].iat[0], 'O128')
        self.assertEqual(multireport.loc['sample3','O type'].iat[0], 'Error! No O type found')

if __name__ == '__main__':
	unittest.main()