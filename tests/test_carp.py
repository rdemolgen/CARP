import argparse, datetime, os, unittest, sys

from pathlib import Path
from src.carp import Carp
from unittest.mock import MagicMock

class Test(unittest.TestCase):

    def setUp(self):
        now = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        self.testDir = Path(f"{os.getcwd()}/tests")
        self.cwd = Path(os.getcwd())

        # define arguments
        self.args = MagicMock(spec=argparse.Namespace)
        self.args.mode = None # required
        self.args.location = None # required
        self.args.proband_id = "WGS_EX2601234" # required
        self.args.prefix = "unittest" # required
        self.args.inDir = self.testDir / "inputs"
        self.args.outDir = self.testDir / "outputs"
        self.args.samples = None
        self.args.genotypes = None
        self.args.vcfFile = None
        self.args.cyto = Path(f"{os.getcwd()}/resources/web_ClinGen_region_curation_list_GRCh38_20250425.tsv")
        self.args.isca = Path(f"{os.getcwd()}/resources/hg38_cytoBand.txt")
        self.args.no_filtering = False
        self.args.min_qual = 30
        self.args.min_dp = 10
        self.args.min_gq = 20
        self.args.min_mq = 40
        self.args.min_qd = 2
    
    def test_(self):
        """Test for"""
        self.args.proband_id = "TwEx_EX2601234"
        carp = Carp(self.args)

        self.assertEqual(1, 2)

    def test_mode_dosage(self):
        """Test that the dosage mode has not be developed yet."""
        self.args.mode = "dosage"
        
        carp = Carp(self.args)

        with self.assertRaises(NotImplementedError):
            carp.run()

    def test_mode_fail(self):
        """Test that the dosage mode has not be developed yet."""
        self.args.mode = "XXX"
        
        carp = Carp(self.args)

        with self.assertRaises(RuntimeError):
            carp.run()

