import datetime, logging, os, pysam, unittest

from pathlib import Path
from src.baf import Baf
from src.utility import Utility

class Test(unittest.TestCase):

    def setUp(self):
        self.now = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        self.testDir = Path(f"{os.getcwd()}/tests")
        self.cwd = Path(os.getcwd())
        self.utility = Utility(self.testDir / "output", "test_", "TwEx_EX2601234", self.now)
        self.filters = {
            "no_filter": False,
            "min_qual": 30,
            "min_dp": 10,
            "min_gq": 20,
            "min_mq": 40,
            "min_qd": 2,
        }
        self.vcfPath = Path(self.testDir / "data" / "TwEx2_EX2601743-TwEx2_EX2601744-TwEx2_EX2601745.vcf.gz")
        self.vcfFile = pysam.VariantFile(str(self.vcfPath))    
        self.baf = Baf(self.utility.logger, self.filters, self.vcfFile)

    def tearDown(self):
        # Reconstruct the log name used in setup_logging
        logName = f"test_TwEx_EX2601234_{self.now}"
        logger = logging.getLogger(logName)

        # Properly close and remove handlers
        for handler in logger.handlers[:]:
            handler.close()
            logger.removeHandler(handler)
    
    def cleanUp(self, in_dir: Path):
        """Remove directory and contents"""
        for f in in_dir.glob("*"):
            f.unlink(missing_ok=True)
        in_dir.rmdir()

    def test_load_vcf(self):
        """
            Test loading a vcf
        """
        vcfFile = self.baf.load_vcf(self.vcfPath)
        self.assertEqual(type(vcfFile), pysam.VariantFile)

    def test_get_samples(self):
        """
            Get samples from user input or vcf
        """
        sample_str = "TwEx2_EX2601743 TwEx2_EX2601744 TwEx2_EX2601745"

        userInput = self.baf.get_samples(sample_str, self.vcfFile)
        fromVcf = self.baf.get_samples(sample_str, self.vcfFile)
        expOutput = sample_str.split(" ")
        self.assertEqual(userInput, expOutput)
        self.assertEqual(fromVcf, expOutput)

    def test_get_location(self):
        """
            Test different location inputs
        """
        chr, start, end = self.baf.get_location("all")
        self.assertEqual((chr, start, end), (None, None, None))
        with self.assertRaises(ValueError):
            chr, start, end = self.baf.get_location("")

    def test_chr_len(self):
        """
            Return chromosome length from vcf header
        """
        chr1 = self.baf.get_chr_len(1, self.vcfFile)
        chr21 = self.baf.get_chr_len(21, self.vcfFile)
        chrX = self.baf.get_chr_len("X", self.vcfFile)
        self.assertEqual(chr1, 248956422)
        self.assertEqual(chr21, 46709983)
        self.assertEqual(chrX, 156040895)
        with self.assertRaises(ValueError):
            self.baf.get_chr_len("A", self.vcfFile)

    def test_get_variants(self):
        """
            Return variant position based on given location criteria, genotype and quality thresholds
        """
        results = self.baf.get_variants(
            "TwEx2_EX2601743",
            21,
            10413729,
            14226797,
            "0/1"
        )
        self.assertEqual(results["id"], "TwEx2_EX2601743")
        self.assertEqual(len(results["positions"]), 67)
        self.assertEqual(results["genotype"], "0/1")
        # Test each filter
        self.filters["no_filter"] = True
        results = self.baf.get_variants(
            "TwEx2_EX2601743",
            21,
            10413729,
            14226797,
            "0/1"
        )
        self.assertEqual(len(results["positions"]), 68)
        self.filters["no_filter"] = False
        self.filters["min_qual"] = 60
        results = self.baf.get_variants(
            "TwEx2_EX2601743",
            21,
            10413729,
            14226797,
            "1/1"
        )
        self.assertEqual(len(results["positions"]), 14)
        self.filters["min_qual"] = 30
        self.filters["min_dp"] = 400
        results = self.baf.get_variants(
            "TwEx2_EX2601743",
            21,
            10413729,
            14226797,
            "0/1"
        )
        self.assertEqual(len(results["positions"]), 6)
        self.filters["min_dp"] = 10
        self.filters["min_gq"] = 99
        results = self.baf.get_variants(
            "TwEx2_EX2601743",
            21,
            10413729,
            14226797,
            "1/1"
        )
        self.assertEqual(len(results["positions"]), 12)
        self.filters["min_gq"] = 20
        self.filters["min_mq"] = 60
        results = self.baf.get_variants(
            "TwEx2_EX2601743",
            21,
            10413729,
            14226797,
            "0/1"
        )
        self.assertEqual(len(results["positions"]), 9)
        self.filters["min_mq"] = 40
        self.filters["min_qd"] = 20
        results = self.baf.get_variants(
            "TwEx2_EX2601743",
            21,
            10413729,
            14226797,
            "0/1"
        )
        self.assertEqual(len(results["positions"]), 4)

    def test_genotypes(self):
        """
            Return list of genotyes if given
        """
        self.assertEqual(self.baf.genotypes(None), None)
        self.assertEqual(self.baf.genotypes("0/1 1/1 0/0"), ["0/1", "1/1", "0/0"])

    def test_get_genotype(self):
        """
            Returns genotype from string as either tuple or verbose string
        """
        ref = self.baf.get_genotype("0/0", False)
        refL = self.baf.get_genotype("0/0", True)
        het = self.baf.get_genotype("0/1", False)
        hetL = self.baf.get_genotype("0/1", True)
        hom = self.baf.get_genotype("1/1", False)
        homL = self.baf.get_genotype("1/1", True)
        unknownL = self.baf.get_genotype("1/2", True)
        self.assertEqual(ref, (0, 0))
        self.assertEqual(refL, "reference")
        self.assertEqual(het, (0, 1))
        self.assertEqual(hetL, "heterozygous")
        self.assertEqual(hom, (1, 1))
        self.assertEqual(homL, "homozygous")
        with self.assertRaises(ValueError):
            self.baf.get_genotype("1/2", False)
        self.assertEqual(unknownL, "all")

    def test_genotype_combinations(self):
        """
            Return genotype combinations to automatically generate standard plots
        """
        twoSamples = [
                ['0/1', '0/0'],
                ['0/1', '0/1'],
                ['0/1', '1/1'],
                ['1/1', '0/0'],
                ['1/1', '0/1'],
                ['1/1', '1/1']
            ]  
        threeSamples = [
                ['0/1', '0/0', '0/1'],
                ['0/1', '0/0', '1/1'],
                ['1/1', '0/0', '0/1'],
                ['1/1', '0/0', '1/1'],
                ['0/1', '0/1', '0/0'],
                ['0/1', '1/1', '0/0'],
                ['1/1', '0/1', '0/0'],
                ['1/1', '1/1', '0/0']
            ]
        self.assertEqual(self.baf.genotype_combinations(2), twoSamples)
        self.assertEqual(self.baf.genotype_combinations(3), threeSamples)
        with self.assertRaises(ValueError):
            self.baf.genotype_combinations("4")
            self.baf.genotype_combinations(None)
