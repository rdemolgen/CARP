import pysam, re, sys
import pandas as pd
from pathlib import Path
from typing import Optional, Tuple, Union

# script classes
from .dosage import Dosage
from .plots import Plots
from .utility import Utility

class Baf:
    """
        This class include functions to identify variants and generate allele fractions
    """

    def __init__(self, logger, filters, vcfFile, chrs):
        self.logger = logger
        self.filters = filters
        self.vcf = self.load_vcf(vcfFile)
        self.chrs = chrs

    @staticmethod
    def get_samples(vcfFile: pysam.VariantFile, samples: str) -> list:
        """
            Return list of samples, assuming proband is the first sample id
        """
        if samples == None:
            print(f"Samples from VCF file: {list(vcfFile.header.samples)}")
            return list(vcfFile.header.samples)
        else:
            print(f"Samples from user input: {samples.split(' ')}")
            return samples.split(' ')

    @staticmethod
    def load_vcf(vcfPath: Path) -> pysam.VariantFile:
        """
            Create a pysam object from vcf file
        """
        print(f"Loading {str(vcfPath)}")
        return pysam.VariantFile(vcfPath)

    @staticmethod
    def get_location(location: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
        """
            Return chrom, start, end from location string <chr:start-end>
        """
        if location == 'all':
            print(f"Genomic Location set to '{location}', processing whole genome")
            return None, None, None
    
        match = re.match(r"^(1?[0-9]|2[0-2]|X|Y|MT)(?::(\d+)-(\d+))?$", location)

        if not match:
            msg = f"Genomic location incorrectly formatted, '{location}'"
            raise ValueError(msg)
        
        return match.groups()

    @staticmethod
    def get_chr_len(vcfFile: pysam.VariantFile, chrom: Union[int, str]) -> int:
        """
            Return chromosome length from vcf header
        """
        chrom_dict = {str(i): i - 1 for i in range(1, 23)}
        chrom_dict.update({'X': 22, 'Y': 23, 'MT': 24})

        try:
            return vcfFile.header.contigs[chrom_dict[str(chrom)]].length
        except KeyError:
            msg = "Chromosome not recognised"
            raise ValueError(msg)

    def get_variants(self, sample: str, chrom: Union[int, str], start: Optional[int]=None, end: Optional[int]=None, genotype: Optional[str]=None) -> dict:
        """
            Return variant position based on given location criteria, genotype and quality thresholds
        """
        gt = Utility.get_genotype(genotype, False)
        positions = []
        if start is None: start = 1
        if end is None: end = Baf.get_chr_len(self.vcf, chrom)

        print(f"Searching {sample} for {gt} variants in chr{str(chrom)}:{int(start)}-{int(end)}")
        for rec in self.vcf.fetch(str(chrom), int(start), int(end)):
            sample_data = rec.samples[sample]

            if sample_data["GT"] != gt and gt != None:
                continue
            # Accept only "PASS" or "." variants
            if not self.filters["no_filtering"]:
                if "PASS" not in rec.filter.keys() and "." not in rec.filter.keys():
                    continue
            try:
                if rec.info['MQ'] < self.filters["min_mq"]:
                    # print("insufficient mapping quality")
                    continue
            except KeyError:
                continue
            try:
                if rec.info['QD'] < self.filters["min_qd"]:
                    # print("insufficient mapping quality")
                    continue
            except KeyError:
                continue
            if rec.qual is not None and rec.qual < self.filters["min_qual"]:
                # print("insufficient quality")
                continue
            try:
                if "DP" in sample_data and sample_data["DP"] < self.filters["min_dp"]:
                    # print("insufficient depth")
                    continue
            except:
                continue
            try:
                if "GQ" in sample_data and sample_data["GQ"] < self.filters["min_gq"]:
                    # print("insufficient GQ")
                    continue
            except Exception as e:
                continue

            positions.append(rec.pos)

        print(f"Number of matching variants: {len(positions)}")
        return {'id': sample, 'positions': positions, 'genotype': genotype}

    @staticmethod
    def genotypes(genotypes: Union[str, None]) -> list:
        """
            Return list of genotyes if given
        """
        if genotypes == None:
            print("No genotypes provided by user")
            return None
        else:
            print(f"Genotypes provided by user: {genotypes.split(' ')}")
            return genotypes.split(' ')

    def genotype_combinations(self, no_samples: int) -> list:
        """
            Return genotype combinations to automatically generate standard plots
        """
        if no_samples == 2:
            return [
                ['0/1', '0/0'],
                ['0/1', '0/1'],
                ['0/1', '1/1'],
                ['1/1', '0/0'],
                ['1/1', '0/1'],
                ['1/1', '1/1']
            ]
        elif no_samples == 3:
            return [
                ['0/1', '0/0', '0/1'],
                ['0/1', '0/0', '1/1'],
                ['1/1', '0/0', '0/1'],
                ['1/1', '0/0', '1/1'],
                ['0/1', '0/1', '0/0'],
                ['0/1', '1/1', '0/0'],
                ['1/1', '0/1', '0/0'],
                ['1/1', '1/1', '0/0']
            ]
        else:
            msg = f"Unexpected number of samples: {str(no_samples)}"
            raise ValueError(msg)

    def intersect_sample_variants(self, sample_pos: list) -> list:
        """
            Find matching variant positions between samples
        """
        intersect_pos = set(sample_pos[0]['positions'])

        for sample in sample_pos[1:]:
            intersect_pos &= set(sample['positions'])
        
        intersect_pos = list(intersect_pos)
        intersect_pos.sort()

        return intersect_pos

    def calc_baf(self, sample: str, chrom: str, positions: list) -> tuple:
        """
            Returns baf and variant position
        """
        allele_fractions = [] # Y axis
        variant_positions = [] # X axis

        print("Calculating variant BAFs")
        total = len(positions)
        bar_width = 40
        for i, pos in enumerate(positions):
            # Update progress bar
            progress = i / total
            filled = int(bar_width * progress)
            bar = "#" * filled + "-" * (bar_width - filled)
            sys.stdout.write(f"\r[{bar}] {i}/{total} ({progress:.0%})")
            sys.stdout.flush()
            for rec in self.vcf.fetch(str(chrom), pos -1, pos):
                ad = rec.samples[sample]['AD']
                if ad and sum(ad) > 0:  # Avoid division by zero
                    baf = ad[1] / sum(ad)  # Alt / (Ref + Alt) 
                    allele_fractions.append(baf)
                    variant_positions.append(pos)
        print("\n")
        return (allele_fractions, variant_positions)

    def get_plot_data(self, sample: str, chrom: str, start: Optional[int]=1, end: Optional[int]=None, genotype: Optional[str]=None) -> tuple:
        """
            Returns list of variat positions and their b-allele frequency
        """
        var_pos = self.get_variants(sample, str(chrom), start, end, genotype)
        plot_data = self.calc_baf(sample, str(chrom), var_pos['positions'])
        return plot_data

    def run_single_plots(self, samples: list, chrom: str, start: int, end: int, outDir: Path, genotypes: Optional[list]=None):
        """
            Automatically generate different genotype plots for each sample
        """
        if genotypes is None:
            genotypes = ['0/1', '1/1', 'all']

        for sample in samples:
            for genotype in genotypes:
                plot_data = self.get_plot_data(sample, chrom, start, end, genotype)
                Plots.plot_baf(plot_data, sample, str(chrom), capture='genome', start=start, end=end, genotype=genotype, outDir=outDir)

    def run_joint_call_plots(self, samples: list, chrom: str, start: int, end: int, outDir: Path, genotypes: Optional[list]=None):
        """
            Generate joint baf plots. This can be a user defined by (sample, genotypes) input or all possible genotype combinations
        """
        if genotypes is None:
            genotypes = self.genotype_combinations(len(samples))
        else:
            genotypes = [genotypes]
        
        for genotype in genotypes:
            samples_var_pos = []
            for sample, geno in zip(samples, genotype):
                samples_var_pos.append(self.get_variants(sample, str(chrom), start, end, geno))
            shared_positions = self.intersect_sample_variants(samples_var_pos)
            plot_data = self.calc_baf(samples[0], str(chrom), shared_positions)
            Plots.plot_baf(plot_data, samples, str(chrom), capture='genome', start=start, end=end, genotype=genotype, outDir=outDir)

    def get_genome_wide_baf(self, sample: str) -> pd.DataFrame: 
        """
            Return baf data for whole genome grouped by chromosome
        """
        chr_baf_dfs = []
        for c in self.chrs:
            baf = self.get_plot_data(sample, c)
            chrom_baf = pd.DataFrame({'chrom': c, 'position': baf[1], 'allele_fraction': baf[0]})
            chr_baf_dfs.append(chrom_baf)
        genome_baf = pd.concat(chr_baf_dfs)
        baf_cumulative_position = Dosage.get_genome_coordinates(genome_baf, 'position')
        grouped_baf = Dosage.group_and_get_cumulative(genome_baf, baf_cumulative_position, 'position')
        return grouped_baf