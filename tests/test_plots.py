import datetime, os, logging, unittest
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.axes import Axes
from pathlib import Path
from types import SimpleNamespace

from src.plots import Plots
from src.utility import Utility

class Test(unittest.TestCase):

    def setUp(self):
        self.now = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        self.date = datetime.datetime.now().strftime("%Y-%m-%d")
        self.cwd = Path(os.getcwd())
        self.testDir = self.cwd / "tests"
        self.outDir = self.testDir / f"output_{self.date}"
        self.outDir.mkdir(parents=True, exist_ok=True)

        self.utility = Utility("test_baf_", "TwEx_EX2601234", self.now, self.testDir / "output")
        self.plots = Plots()

        # Synthetic grouped BAF-like data
        self.baf_chr1 = pd.DataFrame({
            "chrom": ["1", "1"],
            "position": [100, 200],
            "allele_fraction": [0.5, 0.67],
            "genome_coordinate": [100, 200]
        })

        self.baf_chr2 = pd.DataFrame({
            "chrom": ["2", "2"],
            "position": [100, 200],
            "allele_fraction": [0.33, 0.5],
            "genome_coordinate": [300, 400]
        })

        self.grouped_baf = [
            ("1", self.baf_chr1),
            ("2", self.baf_chr2),
        ]

        self.dosage_df = pd.DataFrame({
            "chrom": ["1", "1", "2", "2"],
            "bin_start": [0, 100, 0, 100],
            "bin_end": [100, 200, 100, 200],
            "dosage": [1.0, 1.1, 0.9, 1.2],
            "stdev": [0.1, 0.2, 0.1, 0.2],
            "stdev_pos": [1.1, 1.2, 1.1, 1.2],
            "stdev_neg": [0.9, 0.8, 0.9, 0.8],
            "genome_coordinate": [100, 200, 300, 400]
        })

        self.grouped_dosage = self.dosage_df.groupby("chrom", sort=False)

        self.cnvs = pd.DataFrame({
            "chrom": ["1", "2"],
            "start": [50, 120],
            "end": [150, 220],
            "type": ["Deletion", "Duplication"],
            "size": [100, 100]
        })

        self.regions = pd.DataFrame({
            "chrom": ["1", "2"],
            "start": [60, 130],
            "end": [160, 230]
        })

        # Mock dosage-like object
        self.mock_dosage = SimpleNamespace(
            grouped_read_depth=self.grouped_dosage,
            prefix="test_sample",
            noiseCutoff=0.3,
            cnvs=self.cnvs,
            regions=self.regions
        )

    def tearDown(self):
        plt.close("all")
        logName = f"test_TwEx_EX2601234_{self.now}"
        logger = logging.getLogger(logName)
        for handler in logger.handlers[:]:
            handler.close()
            logger.removeHandler(handler)

    def cleanUp(self, in_dir: Path, content_only: bool = False):
        for f in in_dir.glob("*"):
            f.unlink(missing_ok=True)
        if not content_only:
            in_dir.rmdir()

    def test_plot_baf_ideogram(self):
        fig, ax = plt.subplots()
        result = Plots.plot_baf_ideogram(ax, self.grouped_baf, "genome")

        self.assertIs(result, ax)
        self.assertIsInstance(result, Axes)
        # scatter plots create collections
        self.assertGreater(len(ax.collections), 0)

    def test_make_genome_ideogram(self):
        result = Plots.make_genome_ideogram(self.grouped_baf, self.mock_dosage, "genome")

        self.assertIs(result, plt)
        fig = plt.gcf()
        self.assertEqual(len(fig.axes), 2)
        self.assertEqual(fig.axes[0].get_ylim(), (-0.1, 2.2))

    def test_plot_ideogram_ax(self):
        fig, ax = plt.subplots()
        result = Plots.plot_ideogram_ax(ax, "genome", self.mock_dosage)

        self.assertIs(result, ax)
        self.assertEqual(ax.get_title(), "test_sample, noise cut-off = 0.3")
        self.assertEqual([tick.get_text() for tick in ax.get_xticklabels()], ["1", "2"])

    def test_plot_chr_ideogram_ax(self):
        fig, ax = plt.subplots()
        chr1_dosage = self.grouped_dosage.get_group("1")

        result = Plots.plot_chr_ideogram_ax(
            "sample1",
            ax,
            chr1_dosage,
            "genome",
            "1",
            0,
            self.cnvs,
            self.regions
        )

        self.assertIs(result, ax)
        self.assertEqual(ax.get_ylabel(), "Dosage")
        self.assertEqual(ax.get_xlim()[0], 0)
        # text labels added for sample and tracks
        self.assertGreater(len(ax.texts), 0)
        # patches include track box + cnv/isca rectangles
        self.assertGreater(len(ax.patches), 0)

    def test_plot_chr_ideogram_ax_no_isca_when_not_first(self):
        fig, ax = plt.subplots()
        chr1_dosage = self.grouped_dosage.get_group("1")

        Plots.plot_chr_ideogram_ax(
            "sample1",
            ax,
            chr1_dosage,
            "genome",
            "1",
            1,  # outerind != 0
            self.cnvs,
            self.regions
        )

        texts = [t.get_text() for t in ax.texts]
        self.assertIn("CNVs (1)", texts)
        self.assertNotIn("ISCA regions (1)", texts)

    def test_plot_baf_with_axis(self):
        fig, ax = plt.subplots()
        plot_data = ([0.5, 0.67, 0.33], [100, 200, 300])

        result = Plots.plot_baf(
            plot_data=plot_data,
            sample="sample1",
            chrom="1",
            capture="genome",
            ax=ax,
            genotype="all"
        )

        self.assertIs(result, ax)
        self.assertEqual(ax.get_ylabel(), "Allele Fraction")
        self.assertEqual(ax.get_ylim(), (-0.1, 1.1))
        self.assertGreater(len(ax.collections), 0)
        self.assertEqual(len(ax.lines), 3)  # 3 dashed reference lines

    def test_plot_baf_with_dataframe_input_and_axis(self):
        fig, ax = plt.subplots()

        result = Plots.plot_baf(
            plot_data=self.baf_chr1,
            sample="sample1",
            chrom="1",
            capture="genome",
            ax=ax,
            genotype="all"
        )

        self.assertIs(result, ax)
        self.assertGreater(len(ax.collections), 0)

    def test_plot_baf_saves_png(self):
        plot_data = ([0.5, 0.67, 0.33], [100, 200, 300])

        Plots.plot_baf(
            plot_data=plot_data,
            sample="sample1",
            chrom="1",
            capture="targeted",
            genotype="all",
            outDir=self.outDir
        )

        pngs = list(self.outDir.glob("*_chr1_BAF.png"))
        self.assertEqual(len(pngs), 1)

    def test_get_ylims_default(self):
        ymin, ymax = Plots.get_ylims([0.8, 1.0, 1.5, 2.1])

        self.assertEqual(ymin, -0.1)
        self.assertEqual(ymax, 2.2)

    def test_get_ylims_high_value(self):
        ymin, ymax = Plots.get_ylims([0.8, 1.0, 2.35])

        self.assertEqual(ymin, -0.1)
        self.assertEqual(ymax, 2.45)

    def test_cnvs_track(self):
        fig, ax = plt.subplots()
        ax.set_ylim(0, 2.2)
        ax.set_xlim(0, 300)

        result = Plots.cnvs_track(ax, "1", self.cnvs)

        self.assertIs(result, ax)
        texts = [t.get_text() for t in ax.texts]
        self.assertIn("CNVs (1)", texts)
        self.assertGreater(len(ax.patches), 1)

    def test_isca_track(self):
        fig, ax = plt.subplots()
        ax.set_ylim(0, 2.2)
        ax.set_xlim(0, 300)

        result = Plots.isca_track(ax, "1", self.regions)

        self.assertIs(result, ax)
        texts = [t.get_text() for t in ax.texts]
        self.assertIn("ISCA regions (1)", texts)
        self.assertGreater(len(ax.patches), 1)

    def test_transform_y_point(self):
        fig, ax = plt.subplots()
        ax.set_ylim(0, 2.2)

        y0, h = Plots.transform_y_point(ax, 1.025, 0.1)

        self.assertIsInstance(y0, float)
        self.assertIsInstance(h, float)
        self.assertGreater(h, 0)

    def test_make_chromosome_ideograms_single_sample(self):
        samples = ["sample1"]
        sample_genome_baf = {
            "sample1": self.baf_chr1.groupby("chrom", sort=False)
        }
        dosage_dict = {
            "sample1": self.mock_dosage
        }

        result = Plots.make_chromosome_ideograms(
            samples=samples,
            sample_genome_baf=sample_genome_baf,
            dosage_dict=dosage_dict,
            chr="1",
            capture="genome"
        )

        self.assertIs(result, plt)
        fig = plt.gcf()
        self.assertGreaterEqual(len(fig.axes), 2)

    def test_make_chromosome_ideograms_multiple_samples(self):
        samples = ["sample1", "sample2"]

        baf_df = pd.DataFrame({
            "chrom": ["1", "1", "2", "2"],
            "position": [100, 200, 150, 250],
            "allele_fraction": [0.5, 0.67, 0.33, 0.5],
            "genome_coordinate": [100, 200, 300, 400]
        })

        sample_genome_baf = {
            "sample1": baf_df.groupby("chrom", sort=False),
            "sample2": baf_df.groupby("chrom", sort=False),
        }

        dosage_dict = {
            "sample1": self.mock_dosage,
            "sample2": self.mock_dosage,
        }

        result = Plots.make_chromosome_ideograms(
            samples=samples,
            sample_genome_baf=sample_genome_baf,
            dosage_dict=dosage_dict,
            chr="1",
            capture="genome"
        )

        self.assertIs(result, plt)

    def test_make_chromosome_ideograms_missing_chr_returns_empty_baf_df(self):
        samples = ["sample1", "sample2"]
        grouped_missing = self.baf_chr2.groupby("chrom", sort=False)  # only chr2 present

        sample_genome_baf = {
            "sample1": grouped_missing,
            "sample2": grouped_missing
        }
        dosage_dict = {
            "sample1": self.mock_dosage,
            "sample2": self.mock_dosage
        }

        result = Plots.make_chromosome_ideograms(
            samples=samples,
            sample_genome_baf=sample_genome_baf,
            dosage_dict=dosage_dict,
            chr="1",
            capture="genome"
        )