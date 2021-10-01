from Bio import SeqIO
from numpy import nan
import os
import pandas as pd
import pathlib
from sys import path
import unittest

main_script_path = str(pathlib.Path(pathlib.Path(__file__).parent.absolute()).parent.absolute())
path.insert(0, main_script_path)
from bin import download_cgmlst_scheme
from bin import chewbbaca_input_files
from bin import serotyper_multireport
from bin import get_allele_hashes


class TestSerotypeFinderMultireport(unittest.TestCase):
    def setUpClass():
        assert pathlib.Path('tests/example_output/expected_result_serotype1.csv').exists, \
            "Missing input file for TestSerotypeFinderMultireport test"
        pathlib.Path('test_output').mkdir(exist_ok=True)
    
    def tearDownClass():
        pathlib.Path('test_output').rmdir()

    def test_findallele(self):
        """The code should look for the allele result (wzx in this case) in \
        the results of SerotypeFinder (E. coli serotyper)
        """

        serotyper = serotyper_multireport.SerotypeFinderMultireport(input_files=['tests/example_output/expected_result_serotype1.csv'], 
                                                                    output_file='test_output/serotypefinder_multireport.csv')
        result_csv = pd.read_csv('tests/example_output/expected_result_serotype1.csv', index_col = 0)
        allele = serotyper.find_allele(result_csv, 'wzx')
        self.assertEqual(allele, 'O50/O2')

    def test_getsampleserotype(self):
        """The code should return all the possible serotypes for every locus
        from the results of SerotypeFinder (E. coli serotyper)
        """

        serotyper = serotyper_multireport.SerotypeFinderMultireport(input_files=['tests/example_output/expected_result_serotype1.csv'], 
                                                                    output_file='test_output/serotypefinder_multireport.csv')
        alleles = serotyper.get_sample_serotype('tests/example_output/expected_result_serotype1.csv')
        for allele in ['O50/O2', 'O2','H6']:
            self.assertTrue(allele in alleles, 'The result_serotype.csv obtained by SerotypeFinder is not being properly processed')

    def test_mergemultiplereports(self):
        """Testing properly reporting of O-type. Especial attention to \
        multiple alleles arranged from smaller to larger from SerotypeFinder (E. coli serotyper)"""

        list_files = ['tests/example_output/expected_result_serotype1.csv',
                    'tests/example_output/expected_result_serotype2.csv',
                    'tests/example_output/expected_result_serotype3.csv']
        sample_names = ['sample1', 'sample2', 'sample3']
        serotyper = serotyper_multireport.SerotypeFinderMultireport(input_files=['tests/example_output/expected_result_serotype1.csv',
                                                                                'tests/example_output/expected_result_serotype2.csv',
                                                                                'tests/example_output/expected_result_serotype3.csv'], 
                                                                    output_file='test_output/serotypefinder_multireport.csv')
        serotyper.make_multireport()
        self.assertEqual(serotyper.multireport.loc[0,'O type'], 'O2/O50', serotyper.multireport)
        self.assertEqual(serotyper.multireport.loc[1,'O type'], 'O128')
        self.assertEqual(serotyper.multireport.loc[2,'O type'], 'Error! No O type found')


class TestDownloadcgMLSTSchemes(unittest.TestCase):

    def setUpClass():
        os.system('mkdir -p test_output')

    def tearDownClass():
        os.system('rm -rf test_output')

    def test_output_dir_is_removed_if_download_is_false(self):
        """Testing there is no error if output dir does not exist and that
        it is created then
        """

        os.system('rm -rf test_output')
        cgMLSTschemes_result = download_cgmlst_scheme.cgMLSTSchemes(['salmonella'],
                                                        output_dir='test_output',
                                                        threads=2,
                                                        download_loci=False)
        dir_with_downloaded_scheme = pathlib.Path('test_output', 'salmonella')
        self.assertFalse(dir_with_downloaded_scheme.is_dir())
        self.assertFalse(pathlib.Path('test_output').is_dir())
        

    def test_output_dir_is_not_removed_if_not_empty(self):
        """Testing there is no error if output dir does not exist and that
        it is created then
        """

        os.system('touch test_output/file.txt')
        os.mkdir('test_output/salmonella')
        os.system('touch test_output/salmonella/file.txt')
        cgMLSTschemes_result = download_cgmlst_scheme.cgMLSTSchemes(['salmonella'],
                                                        output_dir='test_output',
                                                        threads=2,
                                                        download_loci=False)
        self.assertTrue(pathlib.Path('test_output').is_dir())
        dir_with_downloaded_scheme = pathlib.Path('test_output', 'salmonella')
        self.assertTrue(dir_with_downloaded_scheme.is_dir())

    def test_salmonella_and_shigella_cgMLSTschemes(self):
        """The Salmonella and Shigella schemes should be properly
        found from Enterobase
        """

        cgMLSTschemes_result = download_cgmlst_scheme.cgMLSTSchemes(['salmonella', 'shigella'],
                                                        output_dir='test_output',
                                                        threads=2,
                                                        download_loci=False)
        expected_result = {'salmonella': \
                            {'source': 'enterobase',
                            'url': 'http://enterobase.warwick.ac.uk/schemes/Salmonella.cgMLSTv2/',
                            'locus_count': 3002,
                            'scheme_description': None},
                        'shigella': \
                            {'source': 'enterobase',
                            'url': 'http://enterobase.warwick.ac.uk/schemes/Escherichia.cgMLSTv1/',
                            'locus_count': 2513,
                            'scheme_description': None}}
        self.assertEqual(len(cgMLSTschemes_result.schemes), 2)
        self.assertTrue('salmonella' in cgMLSTschemes_result.schemes)
        self.assertTrue('shigella' in cgMLSTschemes_result.schemes)
        self.assertEqual(cgMLSTschemes_result.schemes, expected_result)
        dir_with_downloaded_scheme = pathlib.Path('test_output', 'salmonella')
        self.assertFalse(dir_with_downloaded_scheme.is_dir())
        dir_with_downloaded_scheme = pathlib.Path('test_output', 'shigella')
        self.assertFalse(dir_with_downloaded_scheme.is_dir())

    def test_escherichia_cgMLSTschemes(self):
        """Escherichia should download 2 schemes: the default seqsphere and 
        the optional enterobase one (shared with shigella)
        """

        cgMLSTschemes_result = download_cgmlst_scheme.cgMLSTSchemes(['escherichia'],
                                                        output_dir='test_output',
                                                        threads=2,
                                                        download_loci=False)
        expected_result = {'escherichia': \
                            {'source': 'seqsphere',
                            'url': 'https://www.cgmlst.org/ncs/schema/5064703/locus/',
                            'locus_count': 2513,
                            'scheme_description': None},
                        'shigella': \
                            {'source': 'enterobase',
                            'url': 'http://enterobase.warwick.ac.uk/schemes/Escherichia.cgMLSTv1/',
                            'locus_count': 2513,
                            'scheme_description': None}}
        self.assertEqual(len(cgMLSTschemes_result.schemes), 2)
        self.assertTrue('escherichia' in cgMLSTschemes_result.schemes)
        self.assertTrue('shigella' in cgMLSTschemes_result.schemes)
        self.assertEqual(cgMLSTschemes_result.schemes, expected_result)
        dir_with_downloaded_scheme = pathlib.Path('test_output', 'escherichia')
        self.assertFalse(dir_with_downloaded_scheme.is_dir())
        dir_with_downloaded_scheme = pathlib.Path('test_output', 'shigella')
        self.assertFalse(dir_with_downloaded_scheme.is_dir())

    def test_listeria_cgMLSTschemes(self):
        """Listeria should download 2 schemes: the default seqsphere and 
        the optional from BigsDB Pasteur.
        """

        cgMLSTschemes_result = download_cgmlst_scheme.cgMLSTSchemes(['listeria'],
                                                        output_dir='test_output',
                                                        threads=2,
                                                        download_loci=False)
        expected_result = {'listeria': \
                            {'source': 'seqsphere',
                            'url': 'https://www.cgmlst.org/ncs/schema/690488/locus/',
                            'locus_count': 1701,
                            'scheme_description': None},
                        'listeria_optional': \
                            {'source': 'bigsdb_pasteur',
                            'url': 'https://bigsdb.pasteur.fr/api/db/pubmlst_listeria_seqdef/schemes/3',
                            'locus_count': 1748,
                            'scheme_description': 'cgMLST1748'}}
        self.assertEqual(len(cgMLSTschemes_result.schemes), 2)
        self.assertTrue('listeria' in cgMLSTschemes_result.schemes)
        self.assertTrue('listeria_optional' in cgMLSTschemes_result.schemes)
        self.assertEqual(cgMLSTschemes_result.schemes, expected_result)
        dir_with_downloaded_scheme = pathlib.Path('test_output', 'listeria')
        self.assertFalse(dir_with_downloaded_scheme.is_dir())
        dir_with_downloaded_scheme = pathlib.Path('test_output', 'listeria_optional')
        self.assertFalse(dir_with_downloaded_scheme.is_dir())

    def test_warning_when_some_species_not_supported_for_cgMLST(self):
        """If at least one of the species/genera provided for downloading the
        cgMLST scheme is not supported (not on the list of downloadable
        schemes) then a warning should be thrown
        """

        self.assertWarnsRegex(UserWarning, 'not supported', 
                                download_cgmlst_scheme.cgMLSTSchemes, 
                                ['salmonella', 'no_real_species'], 
                                output_dir='output', threads=2,
                                download_loci=False)
    
    def test_error_when_none_of_listed_species_are_supported(self):
        """If none of the species/genera provided for downloading the 
        cgMLST scheme is supported (none of them is on the list of 
        downloadable schemes) then an error should be thrown
        """

        self.assertRaisesRegex(ValueError, 
                                'None of the provided genera is supported for cgMLST', 
                                download_cgmlst_scheme.cgMLSTSchemes, 
                                ['no_real_species'], 
                                output_dir='output', threads=2,
                                download_loci=False)

    def test_genus_given_in_capital_letters(self):
        """The Salmonella scheme should be found even if the genus name
        is given in small letters (no first capital)
        """

        cgMLSTschemes_result = download_cgmlst_scheme.cgMLSTSchemes(['SALMONELLA'],
                                                        output_dir='test_output',
                                                        threads=2,
                                                        download_loci=False)
        expected_result = {'salmonella': \
                            {'source': 'enterobase',
                            'url': 'http://enterobase.warwick.ac.uk/schemes/Salmonella.cgMLSTv2/',
                            'locus_count': 3002,
                            'scheme_description': None}}
        self.assertEqual(len(cgMLSTschemes_result.schemes), 1)
        self.assertTrue('salmonella' in cgMLSTschemes_result.schemes)
        self.assertEqual(cgMLSTschemes_result.schemes, expected_result)

    def test_genus_given_in_firstcapital_letter(self):
        """The Salmonella scheme should be found even if the genus name
        is given in all capital letters (no first capital and other small)
        """

        cgMLSTschemes_result = download_cgmlst_scheme.cgMLSTSchemes(['Salmonella'],
                                                        output_dir='test_output',
                                                        threads=2,
                                                        download_loci=False)
        expected_result = {'salmonella': \
                            {'source': 'enterobase',
                            'url': 'http://enterobase.warwick.ac.uk/schemes/Salmonella.cgMLSTv2/',
                            'locus_count': 3002,
                            'scheme_description': None}}
        self.assertTrue(len(cgMLSTschemes_result.schemes) == 1)
        self.assertTrue('salmonella' in cgMLSTschemes_result.schemes)
        self.assertEqual(cgMLSTschemes_result.schemes, expected_result)

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
                            'locus_count': 7,
                            'scheme_description': 'MLST'}}
        self.assertTrue(len(cgMLSTschemes_result.schemes) == 1)
        self.assertTrue('test_pubmlst' in cgMLSTschemes_result.schemes)
        self.assertEqual(cgMLSTschemes_result.schemes, expected_result)
        dir_with_downloaded_scheme = pathlib.Path('test_output', 'test_pubmlst')
        self.assertTrue(dir_with_downloaded_scheme.is_dir())
        files_in_downloaded_scheme = os.listdir(dir_with_downloaded_scheme)
        fasta_files_in_downloaded_scheme = [file_ for file_ in files_in_downloaded_scheme if str(file_).endswith('.fasta')]
        len(fasta_files_in_downloaded_scheme)
        self.assertEqual(len(fasta_files_in_downloaded_scheme), expected_result['test_pubmlst']['locus_count'])

    def test_bigsdb_pasteur_scheme_is_properly_downloaded(self):
        """A test subfolder should be created and a fasta files per locus
        should be downloaded
        """

        cgMLSTschemes_result = download_cgmlst_scheme.cgMLSTSchemes(['test_bigsdb_pasteur'],
                                                        output_dir='test_output',
                                                        threads=1,
                                                        download_loci=True)
        expected_result = {'test_bigsdb_pasteur': \
                            {'source': 'bigsdb_pasteur',
                            'url': 'https://bigsdb.pasteur.fr/api/db/pubmlst_listeria_seqdef/schemes/2',
                            'locus_count': 7,
                            'scheme_description': 'MLST'}}
        self.assertTrue(len(cgMLSTschemes_result.schemes) == 1)
        self.assertTrue('test_bigsdb_pasteur' in cgMLSTschemes_result.schemes)
        self.assertEqual(cgMLSTschemes_result.schemes, expected_result)
        dir_with_downloaded_scheme = pathlib.Path('test_output', 'test_bigsdb_pasteur')
        self.assertTrue(dir_with_downloaded_scheme.is_dir())
        files_in_downloaded_scheme = os.listdir(dir_with_downloaded_scheme)
        fasta_files_in_downloaded_scheme = [file_ for file_ in files_in_downloaded_scheme if str(file_).endswith('.fasta')]
        len(fasta_files_in_downloaded_scheme)
        self.assertEqual(len(fasta_files_in_downloaded_scheme), expected_result['test_bigsdb_pasteur']['locus_count'])

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
                            'locus_count': 7,
                            'scheme_description': None}}
        self.assertTrue(len(cgMLSTschemes_result.schemes) == 1)
        self.assertTrue('test_enterobase' in cgMLSTschemes_result.schemes)
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
        os.system('mkdir -p test_chewbbaca_per_genus/sample4')
        os.system('mkdir -p test_chewbbaca_per_genus/sample5')
        os.system('mkdir -p test_chewbbaca_per_genus/sample6')
        os.system('mkdir -p test_chewbbaca_per_genus/output')
        with open('test_chewbbaca_per_genus/sample_sheet.yaml', 'w') as file_:
            file_.write('sample1:\n  assembly: sample1.fasta\nsample2:\n  assembly: sample2.fasta\nsample3:\n  assembly: sample3.fasta\nsample4:\n  assembly: sample4.fasta\nsample5:\n  assembly: sample5.fasta\nsample6:\n  assembly: sample6.fasta')
        with open('test_chewbbaca_per_genus/sample1/best_species_hit.txt', 'w') as file_:
            file_.write('salmonella enterica')
        with open('test_chewbbaca_per_genus/sample2/best_species_hit.txt', 'w') as file_:
            file_.write('salmonella')
        with open('test_chewbbaca_per_genus/sample3/best_species_hit.txt', 'w') as file_:
            file_.write('listeria monocytogenes')
        with open('test_chewbbaca_per_genus/sample4/best_species_hit.txt', 'w') as file_:
            file_.write('escherichia coli')
        with open('test_chewbbaca_per_genus/sample5/best_species_hit.txt', 'w') as file_:
            file_.write('shigella flexnerii')
        with open('test_chewbbaca_per_genus/sample6/best_species_hit.txt', 'w') as file_:
            file_.write('fakegenus fakespecies')
        
    def tearDownClass():
        os.system('rm -rf test_chewbbaca_per_genus')
        
    def test_dict_with_samples_per_genus(self):
        """Testing whether the dictionary with samples belonging to a
        genus is created properly for a simple sample
        """

        best_hit_files = ['test_chewbbaca_per_genus/sample1/best_species_hit.txt']
        input_chewbbaca = chewbbaca_input_files.inputChewBBACA(best_hit_kmerfinder_files=best_hit_files,
                                                        sample_sheet='test_chewbbaca_per_genus/sample_sheet.yaml',
                                                        output_dir='test_chewbbaca_per_genus/output')
        input_chewbbaca.make_file_with_samples_per_genus()
        expected_output = {'salmonella': {'samples': ['sample1'], 
                                            'genus_file': 'test_chewbbaca_per_genus/output/salmonella_samples.txt'}}
        actual_output = input_chewbbaca.genera_input_dict
        self.assertEqual(expected_output, actual_output)

    def test_dict_with_samples_per_genus_if_besthitfile_contains_only_genus(self):
        """Testing whether the dictionary with samples belonging to a
        genus is created properly for a sample where only the genus (Salmonella)
        and not the species (enterica) is provided
        """

        best_hit_files = ['test_chewbbaca_per_genus/sample2/best_species_hit.txt']
        input_chewbbaca = chewbbaca_input_files.inputChewBBACA(best_hit_kmerfinder_files=best_hit_files,
                                                        sample_sheet='test_chewbbaca_per_genus/sample_sheet.yaml',
                                                        output_dir='test_chewbbaca_per_genus/output')
        input_chewbbaca.make_file_with_samples_per_genus()
        expected_output = {'salmonella': {'samples': ['sample2'], 
                                            'genus_file': 'test_chewbbaca_per_genus/output/salmonella_samples.txt'}}
        actual_output = input_chewbbaca.genera_input_dict
        self.assertEqual(expected_output, actual_output)

    def test_dict_with_samples_per_genus_if_two_cgmlst_schemes_needed(self):
        """Testing whether the dictionary with samples belonging to a
        genus is created properly for a Listeria, where two cgMLST schemes are needed
        """

        best_hit_files = ['test_chewbbaca_per_genus/sample3/best_species_hit.txt']
        input_chewbbaca = chewbbaca_input_files.inputChewBBACA(best_hit_kmerfinder_files=best_hit_files,
                                                        sample_sheet='test_chewbbaca_per_genus/sample_sheet.yaml',
                                                        output_dir='test_chewbbaca_per_genus/output')
        input_chewbbaca.make_file_with_samples_per_genus()
        expected_output = {'listeria': {'samples': ['sample3'], 
                                        'genus_file': 'test_chewbbaca_per_genus/output/listeria_samples.txt'},
                            'listeria_optional': {'samples': ['sample3'], 
                                                'genus_file': 'test_chewbbaca_per_genus/output/listeria_optional_samples.txt'}}
        actual_output = input_chewbbaca.genera_input_dict
        self.assertEqual(expected_output, actual_output)

    def test_dict_with_samples_for_escherichia(self):
        """Testing whether the dictionary with samples belonging to a
        genus is created properly for a Escherichia, where two cgMLST schemes 
        are needed one of them being Escherichia and the other one Shigella
        """

        best_hit_files = ['test_chewbbaca_per_genus/sample4/best_species_hit.txt']
        input_chewbbaca = chewbbaca_input_files.inputChewBBACA(best_hit_kmerfinder_files=best_hit_files,
                                                        sample_sheet='test_chewbbaca_per_genus/sample_sheet.yaml',
                                                        output_dir='test_chewbbaca_per_genus/output')
        input_chewbbaca.make_file_with_samples_per_genus()
        expected_output = {'escherichia': {'samples': ['sample4'], 
                                        'genus_file': 'test_chewbbaca_per_genus/output/escherichia_samples.txt'},
                            'shigella': {'samples': ['sample4'], 
                                        'genus_file': 'test_chewbbaca_per_genus/output/shigella_samples.txt'}}
        actual_output = input_chewbbaca.genera_input_dict
        self.assertEqual(expected_output, actual_output)

    def test_dict_with_samples_for_escherichia_if_also_shigella_samples_present(self):
        """Testing whether the dictionary with samples belonging to a
        genus is created properly for a Escherichia and Shigella samples
        combined. They should not overwrite each other
        """
        # TODO: This might change if we improve identification and then Shigella only needs to run one schema

        best_hit_files = ['test_chewbbaca_per_genus/sample4/best_species_hit.txt',
                            'test_chewbbaca_per_genus/sample5/best_species_hit.txt']
        input_chewbbaca = chewbbaca_input_files.inputChewBBACA(best_hit_kmerfinder_files=best_hit_files,
                                                        sample_sheet='test_chewbbaca_per_genus/sample_sheet.yaml',
                                                        output_dir='test_chewbbaca_per_genus/output')
        input_chewbbaca.make_file_with_samples_per_genus()
        expected_output = {'escherichia': {'samples': ['sample4', 'sample5'], 
                                        'genus_file': 'test_chewbbaca_per_genus/output/escherichia_samples.txt'},
                            'shigella': {'samples': ['sample4', 'sample5'], 
                                        'genus_file': 'test_chewbbaca_per_genus/output/shigella_samples.txt'}}
        actual_output = input_chewbbaca.genera_input_dict
        self.assertEqual(expected_output, actual_output)

    def test_genus_dict_for_multifiles_multigenera(self):
        """Testing whether the dictionary with sample files is created properly
        when multiple samples from multiple genera are included
        """
        # TODO: This might change if we improve identification and then Shigella only needs to run one schema

        best_hit_files = ['test_chewbbaca_per_genus/sample1/best_species_hit.txt', 
                            'test_chewbbaca_per_genus/sample2/best_species_hit.txt',
                            'test_chewbbaca_per_genus/sample3/best_species_hit.txt',
                            'test_chewbbaca_per_genus/sample4/best_species_hit.txt',
                            'test_chewbbaca_per_genus/sample5/best_species_hit.txt']
        input_chewbbaca = chewbbaca_input_files.inputChewBBACA(best_hit_kmerfinder_files=best_hit_files,
                                                        sample_sheet='test_chewbbaca_per_genus/sample_sheet.yaml',
                                                        output_dir='test_chewbbaca_per_genus/output')
        input_chewbbaca.make_file_with_samples_per_genus()
        expected_output = {'escherichia': {'samples': ['sample4', 'sample5'], 
                                            'genus_file': 'test_chewbbaca_per_genus/output/escherichia_samples.txt'},
                            'listeria': {'samples': ['sample3'], 
                                        'genus_file': 'test_chewbbaca_per_genus/output/listeria_samples.txt'},
                            'listeria_optional': {'samples': ['sample3'], 
                                                'genus_file': 'test_chewbbaca_per_genus/output/listeria_optional_samples.txt'},
                            'salmonella': {'samples': ['sample1', 'sample2'], 
                                            'genus_file': 'test_chewbbaca_per_genus/output/salmonella_samples.txt'},
                            'shigella': {'samples': ['sample4', 'sample5'], 
                                        'genus_file': 'test_chewbbaca_per_genus/output/shigella_samples.txt'}}
        actual_output = input_chewbbaca.genera_input_dict
        self.assertEqual(expected_output, actual_output)

    def test_genus_dict_multifiles_multigenera_ignoring_nonsupported(self):
        """Testing whether the dictionary with sample files is created properly
        when multiple samples from multiple genera are included. However, 
        unsupported genera (in this case fakegenus) should be ignored.
        """

        # TODO: This might change if we improve identification and then Shigella only needs to run one schema

        best_hit_files = ['test_chewbbaca_per_genus/sample1/best_species_hit.txt', 
                            'test_chewbbaca_per_genus/sample2/best_species_hit.txt',
                            'test_chewbbaca_per_genus/sample3/best_species_hit.txt',
                            'test_chewbbaca_per_genus/sample4/best_species_hit.txt',
                            'test_chewbbaca_per_genus/sample5/best_species_hit.txt',
                            'test_chewbbaca_per_genus/sample6/best_species_hit.txt']
        input_chewbbaca = chewbbaca_input_files.inputChewBBACA(best_hit_kmerfinder_files=best_hit_files,
                                                        sample_sheet='test_chewbbaca_per_genus/sample_sheet.yaml',
                                                        output_dir='test_chewbbaca_per_genus/output')
        input_chewbbaca.make_file_with_samples_per_genus()
        expected_output = {'escherichia': {'samples': ['sample4', 'sample5'], 
                                            'genus_file': 'test_chewbbaca_per_genus/output/escherichia_samples.txt'},
                            'listeria': {'samples': ['sample3'], 
                                        'genus_file': 'test_chewbbaca_per_genus/output/listeria_samples.txt'},
                            'listeria_optional': {'samples': ['sample3'], 
                                                'genus_file': 'test_chewbbaca_per_genus/output/listeria_optional_samples.txt'},
                            'salmonella': {'samples': ['sample1', 'sample2'], 
                                            'genus_file': 'test_chewbbaca_per_genus/output/salmonella_samples.txt'},
                            'shigella': {'samples': ['sample4', 'sample5'], 
                                        'genus_file': 'test_chewbbaca_per_genus/output/shigella_samples.txt'}}
        actual_output = input_chewbbaca.genera_input_dict
        self.assertEqual(expected_output, actual_output)

class TestCgmlstHasher(unittest.TestCase):
    
    def test_allele_hashing(self):
        allele = get_allele_hashes.Allele(seq = 'ACTGACTGACTG')
        allele.update_hash()
        expected_output = 'fdfc2bd20f914826b284ce70cfb09e8343d8bba2'
        self.assertEqual(allele.hash, expected_output)

    def test_hash_for_allele_number_is_found(self):
        with open('tests/example_input/locusx.fasta') as file_:
            parsed_fasta = SeqIO.parse(file_, 'fasta')
            fasta_record = next(parsed_fasta)
        allele = get_allele_hashes.Allele(fasta_record=fasta_record)
        allele.extract_attr_from_fasta_record()
        expected_output = '73d83c30aac82425db66e6f074b3278c5ec01404'
        self.assertEqual(allele.hash, expected_output)

    def test_locus_get_subset_alleles(self):
        locus = get_allele_hashes.Locus(fasta_file='tests/example_input/locusx.fasta')
        locus.get_locus_name_from_fasta_file()
        locus.get_subset_alleles_from_fasta(id_numbers=[1,3])
        hash_tbl = locus.make_table_with_hashes(subset=True)
        expected_output = pd.DataFrame.from_dict({'id_number': [1,3], 
                                                    'hash': ['73d83c30aac82425db66e6f074b3278c5ec01404',
                                                                '0273d1c3e65880eacb77e0e9385e7e51c0bdfe9c']})
        self.assertTrue(isinstance(locus.alleles['locusx_1'], get_allele_hashes.Allele))
        self.assertTrue(len(locus.alleles), 2)
        self.assertTrue('locusx_3' in locus.alleles.keys())
        self.assertTrue(all(hash_tbl == expected_output))

    def test_locus_get_all_alleles(self):
        locus = get_allele_hashes.Locus(fasta_file='tests/example_input/locusx.fasta')
        locus.get_locus_name_from_fasta_file()
        locus.get_all_alleles_from_fasta()
        hash_tbl = locus.make_table_with_hashes()
        expected_output = pd.DataFrame.from_dict({'id_number': [1,2,3,4,5], 
                                                    'hash': ['73d83c30aac82425db66e6f074b3278c5ec01404',
                                                                '8e049af5c551058534b88b9dbddbd56d521a6372',
                                                                '0273d1c3e65880eacb77e0e9385e7e51c0bdfe9c',
                                                                '627ede0f675b9236075a32e0958493d7519181d0',
                                                                '0535a35b12c4f532fba151d686b1c2ba41b47024']})
        self.assertTrue(isinstance(locus.alleles['locusx_1'], get_allele_hashes.Allele))
        self.assertTrue(len(locus.alleles), 5)
        self.assertTrue('locusx_2' in locus.alleles.keys())
        self.assertTrue('locusx_3' in locus.alleles.keys())
        self.assertTrue('locusx_4' in locus.alleles.keys())
        self.assertTrue('locusx_5' in locus.alleles.keys())
        self.assertTrue(all(hash_tbl == expected_output))

    @unittest.skipIf(not pathlib.Path('/mnt/db/juno/typing_db/cgmlst/prepared_schemes').exists(),
                    "Skipped in non-RIVM environments because database is not present")
    def test_translate_chewbbaca_report_to_hashes(self):
        expected_output = pd.read_csv('tests/example_output/expected_result_chewbbaca_to_hashes.tsv',
                                    sep="\t", index_col=0)
        sample2_hashes = get_allele_hashes.CgmlstResult(chewbbaca_result_file='tests/example_output/example_result_chewbbaca.tsv', 
                                                        scheme_name='shigella')
        hashed_report = sample2_hashes.get_hashed_cgmlst_report()
        self.assertTrue(all(hashed_report == expected_output))



if __name__ == '__main__':
	unittest.main()