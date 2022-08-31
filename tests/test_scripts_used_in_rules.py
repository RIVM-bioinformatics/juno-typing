from Bio import SeqIO
from numpy import nan
import os
import pandas as pd
import pathlib
from sys import path
import unittest

main_script_path = str(
    pathlib.Path(pathlib.Path(__file__).parent.absolute()).parent.absolute()
)
path.insert(0, main_script_path)
from bin import serotyper_multireport


class TestSerotypeFinderMultireport(unittest.TestCase):
    def setUpClass():
        assert pathlib.Path(
            "tests/example_output/expected_result_serotype1.csv"
        ).exists, "Missing input file for TestSerotypeFinderMultireport test"
        pathlib.Path("test_output").mkdir(exist_ok=True)

    def tearDownClass():
        pathlib.Path("test_output").rmdir()

    def test_findallele(self):
        """The code should look for the allele result (wzx in this case) in \
        the results of SerotypeFinder (E. coli serotyper)
        """

        serotyper = serotyper_multireport.SerotypeFinderMultireport(
            input_files=["tests/example_output/expected_result_serotype1.csv"],
            output_file="test_output/serotypefinder_multireport.csv",
        )
        result_csv = pd.read_csv(
            "tests/example_output/expected_result_serotype1.csv", index_col=0
        )
        allele = serotyper.find_allele(result_csv, "wzx")
        self.assertEqual(allele, "O50/O2")

    def test_getsampleserotype(self):
        """The code should return all the possible serotypes for every locus
        from the results of SerotypeFinder (E. coli serotyper)
        """

        serotyper = serotyper_multireport.SerotypeFinderMultireport(
            input_files=["tests/example_output/expected_result_serotype1.csv"],
            output_file="test_output/serotypefinder_multireport.csv",
        )
        alleles = serotyper.get_sample_serotype(
            "tests/example_output/expected_result_serotype1.csv"
        )
        for allele in ["O50/O2", "O2", "H6"]:
            self.assertTrue(
                allele in alleles,
                "The result_serotype.csv obtained by SerotypeFinder is not being properly processed",
            )

    def test_mergemultiplereports(self):
        """Testing properly reporting of O-type. Especial attention to \
        multiple alleles arranged from smaller to larger from SerotypeFinder (E. coli serotyper)"""

        list_files = [
            "tests/example_output/expected_result_serotype1.csv",
            "tests/example_output/expected_result_serotype2.csv",
            "tests/example_output/expected_result_serotype3.csv",
        ]
        sample_names = ["sample1", "sample2", "sample3"]
        serotyper = serotyper_multireport.SerotypeFinderMultireport(
            input_files=[
                "tests/example_output/expected_result_serotype1.csv",
                "tests/example_output/expected_result_serotype2.csv",
                "tests/example_output/expected_result_serotype3.csv",
            ],
            output_file="test_output/serotypefinder_multireport.csv",
        )
        serotyper.make_multireport()
        self.assertEqual(
            serotyper.multireport.loc[0, "O type"], "O2/O50", serotyper.multireport
        )
        self.assertEqual(serotyper.multireport.loc[1, "O type"], "O128")
        self.assertEqual(
            serotyper.multireport.loc[2, "O type"], "Error! No O type found"
        )


if __name__ == "__main__":
    unittest.main()
