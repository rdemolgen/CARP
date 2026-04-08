import datetime, logging, os, unittest
import numpy as np
import pandas as pd

from collections import OrderedDict
from pandas.core.groupby.generic import DataFrameGroupBy
from pathlib import Path
from typing import Optional

from src.dosage import Dosage
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

        # initialise Dosage class
        self.sample = "WGS_EX2601287"
        self.prefix = "test_"
        self.cytobandsFile = self.cwd / "resources" / "hg38_cytoBand.txt"
        self.iscaFile = self.cwd / "resources" / "web_ClinGen_region_curation_list_GRCh38_20250425.tsv"
        self.readDepthFile = self.testData / "WGS_EX2601287_23GWV3LT3_chr21.coverageBinner6.20000.data"
        self.cnvsFile = self.testData / "cnvs_WGS_EX2601287_23GWV3LT3_chr21"
        self.noiseCutoff = 0.3
        self.dosage = Dosage(
            self.sample,
            self.prefix,
            self.readDepthFile,
            self.cnvsFile,
            self.cytobandsFile,
            self.iscaFile,
            self.chrs,
            self.noiseCutoff
        )

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

    def test_dosage_init(self):
        """

        """
        # read depth
        self.assertEqual(type(self.dosage.readDepth), pd.DataFrame)
        self.assertEqual(self.dosage.readDepth.shape, (1809, 10))
        self.assertEqual(list(self.dosage.readDepth.columns), ['chrom', 'bin_start', 'bin_end', 'dosage', 'stdev', 'unnorm_dosage', 'del_phred', 'dup_phred', 'stdev_pos', 'stdev_neg'])
        self.assertEqual(list(self.dosage.readDepth.iloc[2].tolist()), [21.0, 5240000.0, 5260000.0, 0.9520698308713224, 0.1817044925862543, 0.9526814113266492, -36.22369758765931, -16.34070527666398, 1.1817044925862543, 0.8182955074137457])
        # cnvs
        self.assertEqual(type(self.dosage.cnvs), pd.DataFrame)
        self.assertEqual(self.dosage.cnvs.shape, (1, 11))
        self.assertEqual(list(self.dosage.cnvs.columns), ['chrom', 'start', 'end', 'type', 'bins', 'total_width_bins', 'phred', 'phred_by_bin', 'dosage', 'sample', 'size'])
        self.assertEqual(list(self.dosage.cnvs.iloc[0].tolist()), ['21', np.int64(18680000), np.int64(18700000), 'Duplication', np.int64(1), np.int64(1), np.float64(89.32487792970007), np.float64(89.32487792970007), np.float64(1.3741937661126602), '/home/dnanexus/inputs/input2229880935862243092/WGS_EX2601287_23GWV3LT3.coverageBinner6', np.int64(20000)])
        # isca regions
        self.assertEqual(type(self.dosage.regions), pd.DataFrame)
        self.assertEqual(self.dosage.regions.shape, (56, 5))
        self.assertEqual(list(self.dosage.regions.columns), ['chrom', 'start', 'end', 'HI_evidence', 'TS_evidence'])
        self.assertEqual(list(self.dosage.regions.iloc[10].tolist()), ['15', np.int64(75339446), np.int64(75680568), 'SufficientEvidence', 'NoEvidence'])
        # cytobands
        self.assertEqual(type(self.dosage.cytobands), pd.DataFrame)
        self.assertEqual(self.dosage.cytobands.shape, (1223, 5))
        self.assertEqual(list(self.dosage.cytobands.columns), ['chrom', 'start', 'end', 'band', 'value', ])
        self.assertEqual(list(self.dosage.cytobands.iloc[2].tolist()), ['1', np.int64(5300000), np.int64(7100000), 'p36.31', '#F8F8F8'])
        # acens
        self.assertEqual(type(self.dosage.acens), pd.DataFrame)
        self.assertEqual(self.dosage.acens.shape, (48, 5))
        self.assertEqual(list(self.dosage.acens.columns), ['chrom', 'start', 'end', 'band', 'value', ])
        self.assertEqual(list(self.dosage.acens.iloc[2].tolist()), ['10', np.int64(38000000), np.int64(39800000), 'p11.1', '#E8E8E8'])
        # grouped read data
        grp_rd = self.dosage.grouped_read_depth
        self.assertEqual(type(grp_rd), DataFrameGroupBy)
        self.assertEqual(list(grp_rd.groups.keys()), [21])
        self.assertEqual(list(grp_rd.get_group(21).iloc[2].tolist())[:3], [21.0, 5240000.0, 5260000.0])

    def test_previous_value(self):
        data = OrderedDict([
            ("1", 100),
            ("2", 300),
            ("3", 600),
        ])

        self.assertEqual(Dosage.previous_value(data, "2"), 100)
        self.assertEqual(Dosage.previous_value(data, "3"), 300)

    def test_get_genome_coordinates(self):
        df = pd.DataFrame({
            "chrom": ["1", "1", "2", "2", "3"],
            "bin_end": [100, 200, 150, 300, 400]
        })

        result = Dosage.get_genome_coordinates(df, "bin_end")

        expected = OrderedDict([
            ("1", 200),
            ("2", 500),   # 300 + previous chr max (200)
            ("3", 900),   # 400 + previous cumulative value (500)
        ])

        self.assertEqual(result, expected)

    def test_apply_genome_coordinates_chr1(self):
        row = pd.Series({"chrom": "1", "bin_end": 150})
        genome_dict = OrderedDict([("1", 200), ("2", 500)])

        result = Dosage.apply_genome_coordinates(row, genome_dict, "bin_end")

        self.assertEqual(result, 150)

    def test_apply_genome_coordinates_non_chr1(self):
        row = pd.Series({"chrom": "2", "bin_end": 50})
        genome_dict = OrderedDict([("1", 200), ("2", 500)])

        result = Dosage.apply_genome_coordinates(row, genome_dict, "bin_end")

        self.assertEqual(result, 250)

    def test_remove_centromeres(self):
        read_depth = pd.DataFrame({
            "chrom": ["1", "1", "1", "2"],
            "bin_start": [100, 150, 250, 100],
            "bin_end":   [120, 170, 270, 120],
            "stdev":     [0.1, 0.1, 0.1, 0.1]
        })

        centromeres = pd.DataFrame({
            "chrom": ["1", "1"],
            "start": [140, 160],
            "end":   [180, 220]
        })

        result = self.dosage.remove_centromeres(read_depth, centromeres, ["1"])

        # bin_start 150 falls inside centromere and should be removed
        self.assertEqual(len(result), 3)
        self.assertListEqual(result["bin_start"].tolist(), [100, 250, 100])

    def test_limit_noise(self):
        read_depth = pd.DataFrame({
            "chrom": ["1", "1", "1", "2"],
            "stdev": [0.1, 0.29, 0.3, 0.5]
        })

        high_noise, clean = self.dosage.limit_noise(read_depth, 0.3)

        self.assertEqual(len(high_noise), 2)
        self.assertEqual(len(clean), 2)
        self.assertListEqual(high_noise["stdev"].tolist(), [0.3, 0.5])
        self.assertListEqual(clean["stdev"].tolist(), [0.1, 0.29])

    def test_group_and_get_cumulative(self):
        read_depth = pd.DataFrame({
            "chrom": ["1", "1", "2", "2"],
            "bin_end": [100, 200, 50, 100],
            "bin_start": [0, 100, 0, 50],
            "stdev": [0.1, 0.1, 0.1, 0.1]
        })

        genome_coord_dict = OrderedDict([
            ("1", 200),
            ("2", 300)
        ])

        grouped = Dosage.group_and_get_cumulative(read_depth, genome_coord_dict, "bin_end")

        self.assertIsInstance(grouped, DataFrameGroupBy)
        self.assertListEqual(list(grouped.groups.keys()), ["1", "2"])

        chr1 = grouped.get_group("1")
        chr2 = grouped.get_group("2")

        self.assertListEqual(chr1["genome_coordinate"].tolist(), [100, 200])
        self.assertListEqual(chr2["genome_coordinate"].tolist(), [250, 300])

    def test_process_read_depth_returns_groupby(self):
        result = self.dosage.process_read_depth()

        self.assertIsInstance(result, DataFrameGroupBy)

    def test_load_cnvs_adds_size_column(self):
        cnvs = self.dosage.load_cnvs(self.cnvsFile)

        self.assertIn("size", cnvs.columns)
        self.assertTrue((cnvs["size"] == (cnvs["end"] - cnvs["start"])).all())

    def test_format_cytobands_no_acen_stalk_gvar_in_main_df(self):
        cytobands, acens = self.dosage.format_cytobands(self.cytobandsFile)

        self.assertFalse(cytobands["value"].astype(str).str.contains("acen|stalk|gvar", regex=True).any())
        self.assertTrue((acens["value"] == "#E8E8E8").all())