import argparse
import datetime
import os
import unittest

from pathlib import Path
from unittest.mock import MagicMock, patch

from src.carp import Carp, main


class Test(unittest.TestCase):

    def setUp(self):
        now = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        self.testDir = Path(f"{os.getcwd()}/tests")
        self.cwd = Path(os.getcwd())

        self.args = MagicMock(spec=argparse.Namespace)
        self.args.mode = "baf"
        self.args.location = "1"
        self.args.proband_id = "WGS_EX2601234"
        self.args.prefix = "unittest"
        self.args.inDir = self.testDir / "input"
        self.args.outDir = self.testDir / "output"
        self.args.samples = None
        self.args.genotypes = None
        self.args.vcfFile = None
        self.args.cyto = self.cwd / "resources" / "hg38_cytoBand.txt"
        self.args.isca = self.cwd / "resources" / "web_ClinGen_region_curation_list_GRCh38_20250425.tsv"
        self.args.no_filtering = False
        self.args.min_qual = 30
        self.args.min_dp = 10
        self.args.min_gq = 20
        self.args.min_mq = 40
        self.args.min_qd = 2

    def make_utility_mock(self):
        utility_instance = MagicMock()
        utility_instance.verify_dir.side_effect = lambda x: x
        utility_instance.set_capture.return_value = "genome"
        utility_instance.find_file.return_value = Path("/fake/path/test.vcf.gz")
        utility_instance.logger = MagicMock()
        return utility_instance

    def make_baf_mock(self):
        baf_instance = MagicMock()
        baf_instance.vcf = MagicMock()
        baf_instance.get_genome_wide_baf.return_value = "genome_baf_data"
        baf_instance.get_plot_data.return_value = "chr_baf_data"
        return baf_instance

    @patch("src.carp.Baf")
    @patch("src.carp.Utility")
    def test_init(self, mock_utility_cls, mock_baf_cls):
        mock_utility = self.make_utility_mock()
        mock_utility_cls.return_value = mock_utility

        mock_baf = self.make_baf_mock()
        mock_baf_cls.return_value = mock_baf
        mock_baf_cls.get_samples.return_value = ["WGS_EX2601234"]
        mock_baf_cls.genotypes.return_value = None
        mock_baf_cls.get_location.return_value = ("1", 1, None)

        carp = Carp(self.args)

        self.assertEqual(carp.mode, "baf")
        self.assertEqual(carp.proband_id, "WGS_EX2601234")
        self.assertEqual(carp.prefix, "unittest")
        self.assertEqual(carp.samples, ["WGS_EX2601234"])
        self.assertEqual(carp.genotypes, None)
        self.assertEqual(carp.chrom, "1")
        self.assertEqual(carp.start, 1)
        self.assertEqual(carp.end, None)

        mock_utility_cls.assert_called_once()
        mock_baf_cls.assert_called_once()
        mock_baf_cls.get_samples.assert_called_once()
        mock_baf_cls.genotypes.assert_called_once_with(self.args.genotypes)
        mock_baf_cls.get_location.assert_called_once_with(self.args.location)

    @patch("src.carp.Baf")
    @patch("src.carp.Utility")
    def test_mode_baf_all_location_raises_systemexit(self, mock_utility_cls, mock_baf_cls):
        self.args.mode = "baf"
        self.args.location = "all"

        mock_utility = self.make_utility_mock()
        mock_utility_cls.return_value = mock_utility

        mock_baf = self.make_baf_mock()
        mock_baf_cls.return_value = mock_baf
        mock_baf_cls.get_samples.return_value = ["WGS_EX2601234"]
        mock_baf_cls.genotypes.return_value = None
        mock_baf_cls.get_location.return_value = ("all", None, None)

        carp = Carp(self.args)

        with self.assertRaises(SystemExit):
            carp.run()

    @patch("src.carp.Baf")
    @patch("src.carp.Utility")
    def test_mode_baf_single_sample_no_genotype(self, mock_utility_cls, mock_baf_cls):
        self.args.mode = "baf"

        mock_utility = self.make_utility_mock()
        mock_utility_cls.return_value = mock_utility

        mock_baf = self.make_baf_mock()
        mock_baf_cls.return_value = mock_baf
        mock_baf_cls.get_samples.return_value = ["WGS_EX2601234"]
        mock_baf_cls.genotypes.return_value = None
        mock_baf_cls.get_location.return_value = ("1", 1, None)

        carp = Carp(self.args)
        carp.run()

        mock_baf.run_single_plots.assert_called_once_with(
            ["WGS_EX2601234"], "1", 1, None, self.args.outDir
        )
        mock_baf.run_joint_call_plots.assert_not_called()

    @patch("src.carp.Baf")
    @patch("src.carp.Utility")
    def test_mode_baf_multiple_samples_no_genotype(self, mock_utility_cls, mock_baf_cls):
        self.args.mode = "baf"

        mock_utility = self.make_utility_mock()
        mock_utility_cls.return_value = mock_utility

        mock_baf = self.make_baf_mock()
        mock_baf_cls.return_value = mock_baf
        mock_baf_cls.get_samples.return_value = ["WGS_EX2601234", "WGS_EX2601235"]
        mock_baf_cls.genotypes.return_value = None
        mock_baf_cls.get_location.return_value = ("1", 1, None)

        carp = Carp(self.args)
        carp.run()

        mock_baf.run_single_plots.assert_called_once_with(
            ["WGS_EX2601234", "WGS_EX2601235"], "1", 1, None, self.args.outDir
        )
        mock_baf.run_joint_call_plots.assert_called_once_with(
            ["WGS_EX2601234", "WGS_EX2601235"], "1", 1, None, self.args.outDir
        )

    @patch("src.carp.Baf")
    @patch("src.carp.Utility")
    def test_mode_baf_single_sample_with_single_genotype(self, mock_utility_cls, mock_baf_cls):
        self.args.mode = "baf"
        self.args.genotypes = "0/1"

        mock_utility = self.make_utility_mock()
        mock_utility_cls.return_value = mock_utility

        mock_baf = self.make_baf_mock()
        mock_baf_cls.return_value = mock_baf
        mock_baf_cls.get_samples.return_value = ["WGS_EX2601234"]
        mock_baf_cls.genotypes.return_value = ["het"]
        mock_baf_cls.get_location.return_value = ("1", 1, None)

        carp = Carp(self.args)
        carp.run()

        mock_baf.run_single_plots.assert_called_once_with(
            ["WGS_EX2601234"], "1", 1, None, self.args.outDir, ["het"]
        )
        mock_baf.run_joint_call_plots.assert_not_called()

    @patch("src.carp.Baf")
    @patch("src.carp.Utility")
    def test_mode_baf_matching_samples_and_genotypes_runs_joint(self, mock_utility_cls, mock_baf_cls):
        self.args.mode = "baf"
        self.args.genotypes = "0/1 1/1"

        mock_utility = self.make_utility_mock()
        mock_utility_cls.return_value = mock_utility

        mock_baf = self.make_baf_mock()
        mock_baf_cls.return_value = mock_baf
        mock_baf_cls.get_samples.return_value = ["WGS_EX2601234", "WGS_EX2601235"]
        mock_baf_cls.genotypes.return_value = ["het", "hom"]
        mock_baf_cls.get_location.return_value = ("1", 1, None)

        carp = Carp(self.args)
        carp.run()

        mock_baf.run_joint_call_plots.assert_called_once_with(
            ["WGS_EX2601234", "WGS_EX2601235"], "1", 1, None, self.args.outDir, ["het", "hom"]
        )
        mock_baf.run_single_plots.assert_not_called()

    @patch("src.carp.Baf")
    @patch("src.carp.Utility")
    def test_mode_baf_mismatched_samples_and_genotypes_raises_systemexit(self, mock_utility_cls, mock_baf_cls):
        self.args.mode = "baf"
        self.args.genotypes = "0/1"

        mock_utility = self.make_utility_mock()
        mock_utility_cls.return_value = mock_utility

        mock_baf = self.make_baf_mock()
        mock_baf_cls.return_value = mock_baf
        mock_baf_cls.get_samples.return_value = ["WGS_EX2601234", "WGS_EX2601235"]
        mock_baf_cls.genotypes.return_value = ["het"]
        mock_baf_cls.get_location.return_value = ("1", 1, None)

        carp = Carp(self.args)

        with self.assertRaises(SystemExit):
            carp.run()

    @patch("src.carp.Plots")
    @patch("src.carp.Dosage")
    @patch("src.carp.Baf")
    @patch("src.carp.Utility")
    def test_mode_ideogram_all_calls_pdf_report(self, mock_utility_cls, mock_baf_cls, mock_dosage_cls, mock_plots_cls):
        self.args.mode = "ideogram"
        self.args.location = "all"

        mock_utility = self.make_utility_mock()
        mock_utility.find_file.side_effect = [
            Path("/fake/path/proband.vcf.gz"),
            Path("/fake/path/sample1.20000.data"),
            Path("/fake/path/cnvs_sample1.20000"),
        ]
        mock_utility_cls.return_value = mock_utility

        mock_baf = self.make_baf_mock()
        mock_baf_cls.return_value = mock_baf
        mock_baf_cls.get_samples.return_value = ["WGS_EX2601234"]
        mock_baf_cls.genotypes.return_value = None
        mock_baf_cls.get_location.return_value = ("all", None, None)

        mock_dosage = MagicMock()
        mock_dosage_cls.return_value = mock_dosage

        carp = Carp(self.args)
        carp.run()

        mock_baf.get_genome_wide_baf.assert_called_once_with("WGS_EX2601234")
        mock_dosage_cls.assert_called_once()
        mock_plots_cls.pdf_report.assert_called_once()

    @patch("src.carp.Plots")
    @patch("src.carp.Dosage")
    @patch("src.carp.Baf")
    @patch("src.carp.Utility")
    def test_mode_ideogram_single_chromosome_saves_png(self, mock_utility_cls, mock_baf_cls, mock_dosage_cls, mock_plots_cls):
        self.args.mode = "ideogram"
        self.args.location = "1"

        mock_utility = self.make_utility_mock()
        mock_utility.find_file.side_effect = [
            Path("/fake/path/proband.vcf.gz"),
            Path("/fake/path/sample1.20000.data"),
            Path("/fake/path/cnvs_sample1.20000"),
        ]
        mock_utility_cls.return_value = mock_utility

        mock_baf = self.make_baf_mock()
        mock_baf_cls.return_value = mock_baf
        mock_baf_cls.get_samples.return_value = ["WGS_EX2601234"]
        mock_baf_cls.genotypes.return_value = None
        mock_baf_cls.get_location.return_value = ("1", 1, None)

        mock_dosage = MagicMock()
        mock_dosage_cls.return_value = mock_dosage

        mock_plot_obj = MagicMock()
        mock_plots_cls.make_chromosome_ideograms.return_value = mock_plot_obj

        carp = Carp(self.args)
        carp.run()

        mock_baf.get_plot_data.assert_called_once_with("WGS_EX2601234", "1", 1, None, "all")
        mock_plots_cls.make_chromosome_ideograms.assert_called_once()
        mock_plot_obj.savefig.assert_called_once()
        mock_plot_obj.close.assert_called_once_with("all")

    @patch("src.carp.Baf")
    @patch("src.carp.Utility")
    def test_mode_dosage_raises_not_implemented(self, mock_utility_cls, mock_baf_cls):
        self.args.mode = "dosage"

        mock_utility = self.make_utility_mock()
        mock_utility_cls.return_value = mock_utility

        mock_baf = self.make_baf_mock()
        mock_baf_cls.return_value = mock_baf
        mock_baf_cls.get_samples.return_value = ["WGS_EX2601234"]
        mock_baf_cls.genotypes.return_value = None
        mock_baf_cls.get_location.return_value = ("1", 1, None)

        carp = Carp(self.args)

        with self.assertRaises(NotImplementedError):
            carp.run()

    @patch("src.carp.Baf")
    @patch("src.carp.Utility")
    def test_mode_invalid_raises_runtime_error(self, mock_utility_cls, mock_baf_cls):
        self.args.mode = "xxx"

        mock_utility = self.make_utility_mock()
        mock_utility_cls.return_value = mock_utility

        mock_baf = self.make_baf_mock()
        mock_baf_cls.return_value = mock_baf
        mock_baf_cls.get_samples.return_value = ["WGS_EX2601234"]
        mock_baf_cls.genotypes.return_value = None
        mock_baf_cls.get_location.return_value = ("1", 1, None)

        carp = Carp(self.args)

        with self.assertRaises(RuntimeError):
            carp.run()

    @patch("src.carp.Carp")
    @patch("argparse.ArgumentParser.parse_args")
    def test_main(self, mock_parse_args, mock_carp_cls):
        mock_parse_args.return_value = self.args
        mock_carp = MagicMock()
        mock_carp_cls.return_value = mock_carp

        main()

        mock_carp_cls.assert_called_once_with(self.args)
        mock_carp.run.assert_called_once()