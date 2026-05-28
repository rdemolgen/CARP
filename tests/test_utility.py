import datetime, logging, os, unittest, re

from pathlib import Path
from src.utility import Utility

class Test(unittest.TestCase):

    def setUp(self):
        self.now = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        self.date = datetime.datetime.now().strftime("%Y-%m-%d")
        self.cwd = Path(os.getcwd())
        self.testDir = self.cwd / "tests"
        self.testData = Path(f"/mnt/data1/resources/test_data/carp")
        self.utility = Utility("test_utility_", "TwEx_EX2601234", self.now, self.testDir / "output")
    
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

    def test_get_genotype(self):
        """
            Returns genotype from string as either tuple or verbose string
        """
        ref = self.utility.get_genotype("0/0", False)
        refL = self.utility.get_genotype("0/0", True)
        het = self.utility.get_genotype("0/1", False)
        hetL = self.utility.get_genotype("0/1", True)
        hom = self.utility.get_genotype("1/1", False)
        homL = self.utility.get_genotype("1/1", True)
        unknown = self.utility.get_genotype("1/2", False)
        unknownL = self.utility.get_genotype("1/2", True)
        self.assertEqual(ref, (0, 0))
        self.assertEqual(refL, "reference")
        self.assertEqual(het, (0, 1))
        self.assertEqual(hetL, "heterozygous")
        self.assertEqual(hom, (1, 1))
        self.assertEqual(homL, "homozygous")  
        self.assertEqual(unknown, None)        
        self.assertEqual(unknownL, "all")

    def test_capture(self):
        """Test for setting capture"""
        self.assertEqual(self.utility.set_capture("TwEx_EX2601234"), "exome")
        self.assertEqual(self.utility.set_capture("TwEx2_EX2601234"), "exome")
        self.assertEqual(self.utility.set_capture("WGS_EX2601234"), "genome")
        with self.assertRaises(ValueError):
            self.utility.set_capture("X_EX2601234")

    def test_verify_dir(self):
        """Test for safely setting IO folders"""
        try:
            newDir = f"{str(self.cwd)}/tests/output/new1/new2"
            noDir = None
        
            self.assertEqual(type(newDir), str)
            self.assertEqual(self.utility.verify_dir(newDir), Path(newDir))
            self.assertEqual(self.utility.verify_dir(noDir), self.cwd)
        finally:
            self.cleanUp(Path(newDir))

    def test_find_files_vcfs(self):
        """Test for setting capture"""
        try:
            proband_id = "TwEx2_EX2601746"
            in_dir = Path(self.cwd) / "tests" / "input" / "vcfs"
            in_dir.mkdir(mode=0o777, parents=True, exist_ok=True)

            sample_id = r"(?:WGS_EX\d{7}|TwEx\d*_EX\d{7})"
            vcf_pattern = (
                rf"^(?!.*gnomad_filtered)"
                rf"(?=.*(?:^|-){re.escape(proband_id)}(?:-|\.vcf\.gz$))"
                rf"{sample_id}"
                rf"(?:-{sample_id})*"
                rf"\.vcf\.gz$"
            )

            # no matches
            with self.assertRaises(FileNotFoundError):
                self.utility.find_file(in_dir, vcf_pattern)
            # file with no proband_id in filename
            vcf_file0 = in_dir / "WGS_EX2601234.vcf.gz"
            vcf_file0.write_text("test vcf")           
            with self.assertRaises(FileNotFoundError):
                self.utility.find_file(in_dir, vcf_pattern)           
            # correct match
            vcf_file = in_dir / "TwEx2_EX2601746-TwEx2_EX2601747-TwEx2_EX2601748.vcf.gz"
            vcf_file.write_text("test vcf")
            self.assertEqual(
                self.utility.find_file(in_dir, vcf_pattern).resolve(),
                vcf_file.resolve()
                )
            # matches > 1
            vcf_file2 = in_dir / "TwEx2_EX2601746-TwEx2_EX2601749-TwEx2_EX2601750.vcf.gz"
            vcf_file2.write_text("test vcf")
            with self.assertRaises(RuntimeError):
                self.utility.find_file(in_dir, vcf_pattern)

            data_pattern = rf'^{re.escape(proband_id)}.*\.20000.data$'
            data_file = in_dir / "TwEx2_EX2601746_23GWV3LT3.coverageBinner6.20000.data"
            data_file.write_text("test data")
            self.assertEqual(
                self.utility.find_file(in_dir, data_pattern).resolve(),
                data_file.resolve()
                )
            cnvs_pattern = rf'^cnvs_{re.escape(proband_id)}.*\.20000$'
            cnvs_file = in_dir / "cnvs_TwEx2_EX2601746_22WF3LLT4.20000"
            cnvs_file.write_text("test cnvs")
            self.assertEqual(
                self.utility.find_file(in_dir, cnvs_pattern).resolve(),
                cnvs_file.resolve()
                )
        finally:
            self.cleanUp(in_dir)

    def test_logging(self):
        """Test that logging works"""
        logDir = self.testDir / "logging"
        try:
            logger = self.utility.setup_logging(
                "test_",
                "TwEx_EX2601234",
                self.now,
                logDir
            )
            self.assertEqual(type(logger), logging.Logger)
            logger.info("Add something to log file")
            logFile = logger.handlers[0].baseFilename
            with open(logFile, "r") as f:
                lines = f.readlines()
            self.assertEqual(len(lines), 2)
            self.assertTrue(lines[1].endswith("Add something to log file\n"))
        finally:
            self.cleanUp(logDir)
