import os
from pandas import read_csv
import pathlib
from sys import path
import unittest

main_script_path = str(pathlib.Path(pathlib.Path(__file__).parent.absolute()).parent.absolute())
path.insert(0, main_script_path)
from bin import serotypefinder_multireport
from bin import download_cgmlst_scheme
from bin import chewbbaca_per_genus


class TestSerotypeFinderMultireport(unittest.TestCase):
    def setUpClass():
        assert pathlib.Path('tests/example_output/expected_result_serotype1.csv').exists, \
            "Missing input file for TestSerotypeFinderMultireport test"
    
    def test_findallele(self):
        """The code should look for the allele result (wzx in this case) in \
        the results of SerotypeFinder (E. coli serotyper)
        """

        result_csv = read_csv('tests/example_output/expected_result_serotype1.csv', index_col = 0)
        allele = serotypefinder_multireport.find_allele(result_csv, 'wzx')
        self.assertEqual(allele, 'O50/O2')

    def test_getsampleserotype(self):
        """The code should return all the possible serotypes for every locus
        from the results of SerotypeFinder (E. coli serotyper)
        """

        alleles = serotypefinder_multireport.get_sample_serotype('tests/example_output/expected_result_serotype1.csv')
        for allele in ['O50/O2', 'O2','H6']:
            self.assertTrue(allele in alleles, 'The result_serotype.csv obtained by SerotypeFinder is not being properly processed')

    def test_mergemultiplereports(self):
        """Testing properly reporting of O-type. Especial attention to \
        multiple alleles arranged from smaller to larger from SerotypeFinder (E. coli serotyper)"""

        list_files = ['tests/example_output/expected_result_serotype1.csv',
                    'tests/example_output/expected_result_serotype2.csv',
                    'tests/example_output/expected_result_serotype3.csv']
        sample_names = ['sample1', 'sample2', 'sample3']
        multireport = serotypefinder_multireport.merge_multiple_reports(list_files, sample_names)
        self.assertEqual(multireport.loc['sample1','O type'].iat[0], 'O2/O50')
        self.assertEqual(multireport.loc['sample2','O type'].iat[0], 'O128')
        self.assertEqual(multireport.loc['sample3','O type'].iat[0], 'Error! No O type found')


class TestDownloadcgMLSTSchemes(unittest.TestCase):

    def setUpClass():
        os.system('mkdir -p test_output')

    def tearDownClass():
        os.system('rm -rf test_output')

    # def test_output_dir_is_created(self):
    #     """Testing there is no error if output dir does not exist and that
    #     it is created then
    #     """

    #     os.system('rm -rf test_output')
    #     cgMLSTschemes_result = download_cgmlst_scheme.cgMLSTSchemes(['salmonella'],
    #                                                     output_dir='test_output',
    #                                                     threads=2,
    #                                                     download_loci=False)
    #     self.assertTrue(pathlib.Path('test_output').is_dir())
        
    # def test_salmonella_and_escherichia_cgMLSTschemes(self):
    #     """The Salmonella and Escherichia schemes should be properly
    #     found from pubmlst
    #     """

    #     cgMLSTschemes_result = download_cgmlst_scheme.cgMLSTSchemes(['salmonella', 'escherichia'],
    #                                                     output_dir='test_output',
    #                                                     threads=2,
    #                                                     download_loci=False)
    #     expected_result = {'salmonella': \
    #                         {'source': 'pubmlst',
    #                         'url': 'https://rest.pubmlst.org/db/pubmlst_salmonella_seqdef/schemes/4',
    #                         'scheme_summary_file': 'salmonella_scheme.json',
    #                         'locus_count': 3002,
    #                         'scheme_description': 'cgMLST v2 (Enterobase)'},
    #                     'escherichia': \
    #                         {'source': 'pubmlst',
    #                         'url': 'https://rest.pubmlst.org/db/pubmlst_escherichia_seqdef/schemes/6',
    #                         'scheme_summary_file': 'escherichia_scheme.json',
    #                         'locus_count': 2513,
    #                         'scheme_description': 'cgMLST'}}
    #     self.assertTrue(len(cgMLSTschemes_result.schemes) == 2)
    #     self.assertTrue('salmonella' in cgMLSTschemes_result.schemes)
    #     self.assertTrue('escherichia' in cgMLSTschemes_result.schemes)
    #     self.assertTrue(len(cgMLSTschemes_result.schemes['salmonella']) == 5)
    #     self.assertTrue(len(cgMLSTschemes_result.schemes['escherichia']) == 5)
    #     self.assertEqual(cgMLSTschemes_result.schemes, expected_result)

    # def test_warning_when_some_species_not_supported_for_cgMLST(self):
    #     """If at least one of the species/genera provided for downloading the
    #     cgMLST scheme is not supported (not on the list of downloadable
    #     schemes) then a warning should be thrown
    #     """

    #     self.assertWarnsRegex(UserWarning, 'not supported', 
    #                             download_cgmlst_scheme.cgMLSTSchemes, 
    #                             ['salmonella', 'no_real_species'], 
    #                             output_dir='output', threads=2,
    #                             download_loci=False)
    
    # def test_error_when_none_of_listed_species_are_supported(self):
    #     """If none of the species/genera provided for downloading the 
    #     cgMLST scheme is supported (none of them is on the list of 
    #     downloadable schemes) then an error should be thrown
    #     """

    #     self.assertRaisesRegex(ValueError, 
    #                             'None of the provided species is supported for cgMLST', 
    #                             download_cgmlst_scheme.cgMLSTSchemes, 
    #                             ['no_real_species'], 
    #                             output_dir='output', threads=2,
    #                             download_loci=False)

    # def test_genus_given_in_capital_letters(self):
    #     """The Salmonella scheme should be found even if the genus name
    #     is given in small letters (no first capital)
    #     """

    #     cgMLSTschemes_result = download_cgmlst_scheme.cgMLSTSchemes(['SALMONELLA'],
    #                                                     output_dir='test_output',
    #                                                     threads=2,
    #                                                     download_loci=False)
    #     expected_result = {'salmonella': \
    #                         {'source': 'pubmlst',
    #                         'url': 'https://rest.pubmlst.org/db/pubmlst_salmonella_seqdef/schemes/4',
    #                         'scheme_summary_file': 'salmonella_scheme.json',
    #                         'locus_count': 3002,
    #                         'scheme_description': 'cgMLST v2 (Enterobase)'}}
    #     self.assertTrue(len(cgMLSTschemes_result.schemes) == 1)
    #     self.assertTrue('salmonella' in cgMLSTschemes_result.schemes)
    #     self.assertTrue(len(cgMLSTschemes_result.schemes['salmonella']) == 5)
    #     self.assertEqual(cgMLSTschemes_result.schemes, expected_result)

    # def test_genus_given_in_firstcapital_letter(self):
    #     """The Salmonella scheme should be found even if the genus name
    #     is given in all capital letters (no first capital and other small)
    #     """

    #     cgMLSTschemes_result = download_cgmlst_scheme.cgMLSTSchemes(['Salmonella'],
    #                                                     output_dir='test_output',
    #                                                     threads=2,
    #                                                     download_loci=False)
    #     expected_result = {'salmonella': \
    #                         {'source': 'pubmlst',
    #                         'url': 'https://rest.pubmlst.org/db/pubmlst_salmonella_seqdef/schemes/4',
    #                         'scheme_summary_file': 'salmonella_scheme.json',
    #                         'locus_count': 3002,
    #                         'scheme_description': 'cgMLST v2 (Enterobase)'}}
    #     self.assertTrue(len(cgMLSTschemes_result.schemes) == 1)
    #     self.assertTrue('salmonella' in cgMLSTschemes_result.schemes)
    #     self.assertTrue(len(cgMLSTschemes_result.schemes['salmonella']) == 5)
    #     self.assertEqual(cgMLSTschemes_result.schemes, expected_result)

    @unittest.skipIf(not pathlib.Path('/mnt/scratch_dir/hernanda').exists(),
                    "Skipped in non-RIVM environments (for sake of time)")
    def test_pubmlst_scheme_is_properly_downloaded(self):
        """A test subfolder should be created and a fasta files per locus
        should be downloaded
        """

        cgMLSTschemes_result = download_cgmlst_scheme.cgMLSTSchemes(['test_pubmlst'],
                                                        output_dir='test_output',
                                                        threads=1,
                                                        download_loci=True)
        expected_result = {'test_pubmlst': \
                            {'source': 'pubmlst',
                            'url': 'https://rest.pubmlst.org/db/pubmlst_salmonella_seqdef/schemes/2',
                            'scheme_summary_file': 'test_scheme.json',
                            'locus_count': 7,
                            'scheme_description': 'MLST'}}
        self.assertTrue(len(cgMLSTschemes_result.schemes) == 1)
        self.assertTrue('test_pubmlst' in cgMLSTschemes_result.schemes)
        self.assertTrue(len(cgMLSTschemes_result.schemes['test_pubmlst']) == 5)
        self.assertEqual(cgMLSTschemes_result.schemes, expected_result)
        dir_with_downloaded_scheme = pathlib.Path('test_output', 'test_pubmlst')
        self.assertTrue(dir_with_downloaded_scheme.is_dir())
        files_in_downloaded_scheme = os.listdir(dir_with_downloaded_scheme)
        fasta_files_in_downloaded_scheme = [file_ for file_ in files_in_downloaded_scheme if str(file_).endswith('.fasta')]
        len(fasta_files_in_downloaded_scheme)
        self.assertEqual(len(fasta_files_in_downloaded_scheme), expected_result['test_pubmlst']['locus_count'])

    @unittest.skipIf(not pathlib.Path('/mnt/scratch_dir/hernanda').exists(),
                    "Skipped in non-RIVM environments (for sake of time)")
    def test_enterobase_scheme_is_properly_downloaded(self):
        """A test subfolder should be created and a fasta files per locus
        should be downloaded
        """

        cgMLSTschemes_result = download_cgmlst_scheme.cgMLSTSchemes(['test_enterobase'],
                                                        output_dir='test_output',
                                                        threads=1,
                                                        download_loci=True)
        expected_result = {'test_enterobase': \
                            {'source': 'enterobase',
                            'url': 'https://enterobase.warwick.ac.uk/schemes/Yersinia.Achtman7GeneMLST/',
                            'scheme_summary_file': 'test_scheme.json',
                            'locus_count': 7,
                            'scheme_description': None}}
        self.assertTrue(len(cgMLSTschemes_result.schemes) == 1)
        self.assertTrue('test_enterobase' in cgMLSTschemes_result.schemes)
        self.assertTrue(len(cgMLSTschemes_result.schemes['test_enterobase']) == 5)
        self.assertEqual(cgMLSTschemes_result.schemes, expected_result)
        dir_with_downloaded_scheme = pathlib.Path('test_output', 'test_enterobase')
        self.assertTrue(dir_with_downloaded_scheme.is_dir())
        files_in_downloaded_scheme = os.listdir(dir_with_downloaded_scheme)
        fasta_files_in_downloaded_scheme = [file_ for file_ in files_in_downloaded_scheme if str(file_).endswith('.fasta')]
        len(fasta_files_in_downloaded_scheme)
        self.assertEqual(len(fasta_files_in_downloaded_scheme), expected_result['test_enterobase']['locus_count'])


class TestChewbbacaPerGenus(unittest.TestCase):
    
    def setUpClass():
        os.system('mkdir -p test_chewbbaca_per_genus')
        os.system('mkdir -p test_chewbbaca_per_genus/sample1')
        os.system('mkdir -p test_chewbbaca_per_genus/sample2')
        os.system('mkdir -p test_chewbbaca_per_genus/sample3')
        os.system('mkdir -p test_chewbbaca_per_genus/output')
        os.system('mkdir -p test_chewbbaca_per_genus/db_cgmlst/salmonella')
        os.system('mkdir -p test_chewbbaca_per_genus/db_cgmlst/escherichia')
        with open('test_chewbbaca_per_genus/sample_sheet.yaml', 'w') as file_:
            file_.write('sample1:\n  assembly: sample1.fasta\nsample2:\n  assembly: sample2.fasta\nsample3:\n  assembly: sample3.fasta')
        with open('test_chewbbaca_per_genus/sample1/best_species_hit.txt', 'w') as file_:
            file_.write('salmonella')
        with open('test_chewbbaca_per_genus/sample2/best_species_hit.txt', 'w') as file_:
            file_.write('escherichia coli')
        with open('test_chewbbaca_per_genus/sample3/best_species_hit.txt', 'w') as file_:
            file_.write('streptococcus pneumoniae')
            
    def tearDownClass():
        os.system('rm -rf test_chewbbaca_per_genus')

    def test_chewbbaca_parameters(self):
        """Testing whether the parameters used to run chewbbaca are properly 
        created for given best_species_hit. Ran in the best 'perfect' scenario
        """

        best_hit_files = ['test_chewbbaca_per_genus/sample2/best_species_hit.txt']
        chewbbaca_run = chewbbaca_per_genus.runChewBBACA(best_hit_kmerfinder_files=best_hit_files,
                                                        sample_sheet='test_chewbbaca_per_genus/sample_sheet.yaml',
                                                        output_dir='test_chewbbaca_per_genus/output',
                                                        cgmlst_db_dir='test_chewbbaca_per_genus/db_cgmlst/',
                                                        threads=1)
        expected_output = {'escherichia': {'samples': ['sample2'], 
                                            'genus_file': 'test_chewbbaca_per_genus/output/escherichia', 
                                            'cgmlst_scheme': 'test_chewbbaca_per_genus/db_cgmlst/escherichia'}}
        actual_output = chewbbaca_run.chewbbaca_parameters
        self.assertEqual(expected_output, actual_output)
        
    def test_chewbbaca_parameters_when_non_supported_species(self):
        """Testing whether the parameters used to run chewbbaca are properly 
        created for given best_species_hit. The sample tested here is not supported
        """

        best_hit_files = ['test_chewbbaca_per_genus/sample3/best_species_hit.txt']
        chewbbaca_run = chewbbaca_per_genus.runChewBBACA(best_hit_kmerfinder_files=best_hit_files,
                                                        sample_sheet='test_chewbbaca_per_genus/sample_sheet.yaml',
                                                        output_dir='test_chewbbaca_per_genus/output',
                                                        cgmlst_db_dir='test_chewbbaca_per_genus/db_cgmlst/',
                                                        threads=1)
        expected_output = {'streptococcus': {'samples': ['sample3'], 
                                            'genus_file': 'test_chewbbaca_per_genus/output/streptococcus', 
                                            'cgmlst_scheme': None}}
        actual_output = chewbbaca_run.chewbbaca_parameters
        self.assertEqual(expected_output, actual_output)

    def test_chewbbaca_parameters_work_when_only_genus_no_species_in_besthit_files(self):
        """Testing whether the parameters used to run chewbbaca are properly 
        created even when the besthit file just contains the genus and not
        the species
        """

        best_hit_files = ['test_chewbbaca_per_genus/sample1/best_species_hit.txt']
        chewbbaca_run = chewbbaca_per_genus.runChewBBACA(best_hit_kmerfinder_files=best_hit_files,
                                                        sample_sheet='test_chewbbaca_per_genus/sample_sheet.yaml',
                                                        output_dir='test_chewbbaca_per_genus/output',
                                                        cgmlst_db_dir='test_chewbbaca_per_genus/db_cgmlst/',
                                                        threads=1)
        expected_output = {'salmonella': {'samples': ['sample1'], 
                                            'genus_file': 'test_chewbbaca_per_genus/output/salmonella', 
                                            'cgmlst_scheme': 'test_chewbbaca_per_genus/db_cgmlst/salmonella'}}
        actual_output = chewbbaca_run.chewbbaca_parameters
        self.assertEqual(expected_output, actual_output)
        
    def test_chewbbaca_parameters_multifiles_multigenera(self):
        """Testing whether the parameters used to run chewbbaca are properly 
        created for different best_hit files for multiple genera in one run
        """

        best_hit_files = ['test_chewbbaca_per_genus/sample1/best_species_hit.txt', 
                            'test_chewbbaca_per_genus/sample2/best_species_hit.txt',
                            'test_chewbbaca_per_genus/sample3/best_species_hit.txt']
        chewbbaca_run = chewbbaca_per_genus.runChewBBACA(best_hit_kmerfinder_files=best_hit_files,
                                                        sample_sheet='test_chewbbaca_per_genus/sample_sheet.yaml',
                                                        output_dir='test_chewbbaca_per_genus/output',
                                                        cgmlst_db_dir='test_chewbbaca_per_genus/db_cgmlst/',
                                                        threads=1)
        expected_output = {'escherichia': {'samples': ['sample2'], 
                                            'genus_file': 'test_chewbbaca_per_genus/output/escherichia', 
                                            'cgmlst_scheme': 'test_chewbbaca_per_genus/db_cgmlst/escherichia'}, 
                            'streptococcus': {'samples': ['sample3'], 
                                            'genus_file': 'test_chewbbaca_per_genus/output/streptococcus', 
                                            'cgmlst_scheme': None}, 
                            'salmonella': {'samples': ['sample1'], 
                                            'genus_file': 'test_chewbbaca_per_genus/output/salmonella', 
                                            'cgmlst_scheme': 'test_chewbbaca_per_genus/db_cgmlst/salmonella'}}
        actual_output = chewbbaca_run.chewbbaca_parameters
        self.assertEqual(expected_output, actual_output)


if __name__ == '__main__':
	unittest.main()