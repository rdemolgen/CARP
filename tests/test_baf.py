import datetime, logging, os, pysam, unittest

from pandas.core.groupby.generic import DataFrameGroupBy
from pathlib import Path
from typing import Optional
from src.baf import Baf
from src.utility import Utility

class Test(unittest.TestCase):

    def setUp(self):
        self.now = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        self.date = datetime.datetime.now().strftime("%Y-%m-%d")
        self.cwd = Path(os.getcwd())
        self.testDir = self.cwd / "tests"
        self.testData = Path(f"/mnt/data1/resources/test_data/carp")
        self.chrs = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
        self.utility = Utility("test_baf_", "TwEx_EX2601234", self.now, self.testDir / "output")
        self.outDir = self.utility.verify_dir(self.testDir / f"output_{self.date}")
        self.filters = {
            "no_filtering": False,
            "min_qual": 30,
            "min_dp": 10,
            "min_gq": 20,
            "min_mq": 40,
            "min_qd": 2,
        }
        self.vcfPath = Path(self.testData / "TwEx2_EX2601743-TwEx2_EX2601744-TwEx2_EX2601745.vcf.gz")   
        self.baf = Baf(self.utility.logger, self.filters, self.vcfPath, self.chrs)
        self.vcf = self.baf.vcf

    def tearDown(self):
        # Reconstruct the log name used in setup_logging
        logName = f"test_TwEx_EX2601234_{self.now}"
        logger = logging.getLogger(logName)

        # Properly close and remove handlers
        for handler in logger.handlers[:]:
            handler.close()
            logger.removeHandler(handler)
    
    def cleanUp(self, in_dir: Path, content_only: Optional[bool]=False):
        """Remove directory and contents"""
        for f in in_dir.glob("*"):
            f.unlink(missing_ok=True)
        if not content_only:
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

        userInput = self.baf.get_samples(self.vcf, sample_str)
        fromVcf = self.baf.get_samples(self.vcf, None)
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
        chr, start, end = self.baf.get_location("1:1234-5678")
        self.assertEqual((chr, start, end), ("1", "1234", "5678"))

    def test_chr_len(self):
        """
            Return chromosome length from vcf header
        """
        chr1 = self.baf.get_chr_len(self.vcf, 1)
        chr21 = self.baf.get_chr_len(self.vcf, 21)
        chrX = self.baf.get_chr_len(self.vcf, "X")
        self.assertEqual(chr1, 248956422)
        self.assertEqual(chr21, 46709983)
        self.assertEqual(chrX, 156040895)
        with self.assertRaises(ValueError):
            self.baf.get_chr_len(self.vcf, "A")

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
        self.assertEqual(len(results["positions"]), 67)
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
        # test that ./. genotypes are ignored
        for gt in ["0/0", "0/1", "1/1", "all"]:
            results = self.baf.get_variants(
                "TwEx2_EX2601743",
                21,
                45990497,
                45990497,
                gt
            )  
            self.assertEqual(len(results["positions"]), 0)

    def test_genotypes(self):
        """
            Return list of genotyes if given
        """
        self.assertEqual(self.baf.genotypes(None), None)
        self.assertEqual(self.baf.genotypes("0/1 1/1 0/0"), ["0/1", "1/1", "0/0"])

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

    def test_intersect_sample_variants(self):
        """
            Find matching variant positions between samples
        """
        samples = ["TwEx2_EX2601743", "TwEx2_EX2601744"]
        genotypes = ["1/1", "0/1"]
        chrom, start, end = "21", 10413729, 14226797,
        samples_var_pos = []
        for sample, genotype in zip(samples, genotypes):
            samples_var_pos.append(self.baf.get_variants(sample, chrom, start, end, genotype))
        shared_positions = self.baf.intersect_sample_variants(samples_var_pos)
        self.assertEqual(type(shared_positions), list)
        self.assertEqual(len(shared_positions), 9)

    def test_calc_baf(self):
        """
            Returns baf and variant position
        """
        positions = [10413783, 10602110, 14108913, 14109044, 14185882, 14186025, 14186143, 14210899, 14211011]
        results = self.baf.calc_baf("TwEx2_EX2601744", "21", positions)
        baf = results[0]
        pos = results[1]
        self.assertEqual(type(results), tuple)
        self.assertEqual(len(baf), 10)
        self.assertEqual(len(pos), 10)
        self.assertEqual(baf, [0.17229729729729729, 0.8170289855072463, 0.7935943060498221, 0.5125, 0.4462809917355372, 0.42857142857142855, 0.40707964601769914, 0.27906976744186046, 0.46206896551724136, 0.6])
        # indel and snp that are covered by the same position, 10413783 
        self.assertEqual(pos, [10413783, 10413783, 10602110, 14108913, 14109044, 14185882, 14186025, 14186143, 14210899, 14211011])

    def test_get_plot_data(self):
        """
            Returns list of variat positions and their b-allele frequency 
        """
        samples = self.baf.get_samples(self.vcf, None)
        sample = samples[0]
        chrom = "21"
        start = 10413729
        end = 14226797
        genotypes = None
        af, pos = self.baf.get_plot_data(sample, chrom, start, end, genotypes)
        self.assertEqual(len(af), 86)
        self.assertEqual(len(pos), 86)
        self.assertEqual((af[0], pos[0]), (0.2125984251968504,10413733))

    def test_run_single_plots(self):
        """
           Automatically generate different genotype plots for each sample  
        """
        samples = self.baf.get_samples(self.vcf, None)
        chrom = "21"
        start = 10413729
        end = 14226797
        outDir = self.utility.verify_dir(self.outDir / "run_single_plot")
        genotypes = None
        self.baf.run_single_plots(samples, chrom, start, end, outDir, genotypes)
        plots = [f for f in outDir.glob("*")]
        plots.sort()
        self.assertEqual(len(plots), 9)
        self.assertAlmostEqual(plots[0].name, "TwEx2_EX2601743_all_chr21.10413729-14226797_BAF.png")
        self.cleanUp(outDir, True)
        samples = samples[:2]
        genotypes = ["0/1", "1/1"]
        self.baf.run_single_plots(samples, chrom, start, end, outDir, genotypes)
        plots = [f for f in outDir.glob("*")]
        plots.sort()
        self.assertEqual(len(plots), 4)
        self.assertAlmostEqual(plots[0].name, "TwEx2_EX2601743_het_chr21.10413729-14226797_BAF.png")
        self.cleanUp(outDir)

    def test_run_joint_call_plots(self):
        """
            Generate joint baf plots. This can be a user defined by (sample, genotypes) input or all possible genotype combinations
        """
        samples = self.baf.get_samples(self.vcf, None)
        indexSample = samples[0]
        chrom = "21"
        start = 10413729
        end = 14226797
        outDir = self.utility.verify_dir(self.outDir / "joint_call_plots")
        genotypes = None
        self.baf.run_joint_call_plots(indexSample, samples, chrom, start, end, outDir, genotypes)
        plots = [f for f in outDir.glob("*")]
        plots.sort()
        self.assertEqual(len(plots), 8)
        self.assertAlmostEqual(plots[0].name, "TwEx2_EX2601743_het_TwEx2_EX2601744_het_TwEx2_EX2601745_ref_chr21.10413729-14226797_BAF.png")
        self.cleanUp(outDir, True)
        genotypes = ["0/1", "0/0", "1/1"]
        self.baf.run_joint_call_plots(indexSample, samples, chrom, start, end, outDir, genotypes)
        plots = [f for f in outDir.glob("*")]
        self.assertEqual(len(plots), 1)
        self.assertAlmostEqual(plots[0].name, "TwEx2_EX2601743_het_TwEx2_EX2601744_ref_TwEx2_EX2601745_hom_chr21.10413729-14226797_BAF.png")
        self.cleanUp(outDir)
    
    def test_get_genome_wide_baf(self):
        """
            Return baf data for whole genome grouped by chromosome
        """
        samples = self.baf.get_samples(self.vcf, None)[0]
        results = self.baf.get_genome_wide_baf(samples)
        self.assertEqual(type(results), DataFrameGroupBy)
        self.assertEqual(len(results), 24)
        self.assertEqual(list(results.groups.keys()), self.chrs)

