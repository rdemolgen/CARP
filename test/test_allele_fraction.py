import unittest
from unittest.mock import MagicMock
import argparse
import sys, os, shutil

sys.path.insert(0, f'{os.getcwd()}/src')
from allele_fraction import Allele_Fraction

class TestAlleleFraction(unittest.TestCase):
    testOutDir = f"{os.getcwd()}/test/output/"
    
    @classmethod
    def setUpClass(cls):
        # Mock argparse Namespace object to simulate command-line arguments
        cls.args = MagicMock(spec=argparse.Namespace)
        
        # Set default values for the mock arguments (you can change these as needed in each test)
        cls.args.vcfFile = "test/data/test_allele_fraction.vcf.gz"
        cls.args.location = "1:1-2000"
        cls.args.samples = None
        cls.args.genotypes = None
        cls.args.no_filtering = False
        cls.args.min_qual = 30
        cls.args.min_dp = 10
        cls.args.min_gq = 20
        cls.args.outDir = 'test/output/'

        # output for test plots
        os.mkdir(f'{cls.testOutDir}')

    def instantiate_AF(self):
        # Instantiate Allele_Fraction with the mock arguments
        return Allele_Fraction(
            vcfFile=self.args.vcfFile,
            location=self.args.location,
            samples=self.args.samples,
            genotypes=self.args.genotypes,
            no_filtering=self.args.no_filtering,
            min_qual=self.args.min_qual,
            min_dp=self.args.min_dp,
            min_gq=self.args.min_gq,
            outDir=self.args.outDir
        )

    def test_default_samples(self):
        AF = self.instantiate_AF()

        self.assertEqual(len(AF.samples), 3)

    def test_tweaked_samples(self):
        self.args.samples = "sample1 sample2"
        AF = self.instantiate_AF()

        self.assertEqual(len(AF.samples), 2)

    @classmethod
    def tearDownClass(cls):
        print(f"Deleting... {cls.testOutDir}")
        shutil.rmtree(cls.testOutDir)

# To run the tests
if __name__ == '__main__':
    unittest.main()
