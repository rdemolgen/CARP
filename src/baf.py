import pysam, re
from pathlib import Path
from typing import Optional, Tuple, Union

class Baf:
    """
        This class include functions to identify variants and generate allele fractions
    """

    def __init__(self, logger, filters, vcfFile):
        self.logger = logger
        self.filters = filters
        self.vcfFile = vcfFile

    def get_samples(self, samples: str) -> list:
        """
            Return list of samples, assuming proband is the first sample id
        """
        if samples == None:
            self.logger.info(list(self.vcfFile.header.samples))
            return list(f"Samples from VCF file: {self.vcfFile.header.samples}")
        else:
            self.logger.info(f"Samples from user input: {samples.split(' ')}")
            return samples.split(' ')

    def load_vcf(self, vcfPath: Path) -> pysam.VariantFile:
        """
            Create a pysam object from vcf file
        """
        self.logger.info(f"Loading {str(vcfPath)}")
        return pysam.VariantFile(vcfPath)

    def get_location(self, location: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
        """
            Return chrom, start, end from location string <chr:start-end>
        """
        if location == 'all':
            self.logger.info(f"Genomic Location set to '{location}', processing whole genome")
            return None, None, None
    
        match = re.match(r"^(1?[0-9]|2[0-2]|X|Y|MT)(?::(\d+)-(\d+))?$", location)

        if not match:
            msg = f"Genomic location incorrectly formatted, '{location}'"
            self.logger.error(msg)
            raise ValueError(msg)
        
        return match.groups()

    def get_chr_len(self, chrom: Union[int, str]) -> int:
        """
            Return chromosome length from vcf header
        """
        chrom_dict = {str(i): i - 1 for i in range(1, 23)}
        chrom_dict.update({'X': 22, 'Y': 23, 'MT': 24})

        try:
            return self.vcfFile.header.contigs[chrom_dict[str(chrom)]].length
        except KeyError:
            msg = "Chromosome not recognised"
            self.logger.error(msg)
            raise ValueError(msg)

    def get_variants(self, sample: str, chrom: Union[int, str], start: Optional[int]=None, end: Optional[int]=None, genotype: Optional[str]=None) -> dict:
        """
            Return variant position based on given location criteria, genotype and quality thresholds
        """
        gt = self.get_genotype(genotype, False)
        positions = []
        if start is None: start = 1
        if end is None: end = self.get_chr_len(chrom)

        self.logger.info(f"Searching {sample} for {gt} variants in chr{str(chrom)}:{int(start)}-{int(end)} using the following filters:")
        for filter, value in self.filters.items():
            self.logger.info(f"\t{filter}: {value}")
        for rec in self.vcfFile.fetch(str(chrom), int(start), int(end)):
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

        self.logger.info(f"Number of matching variants: {len(positions)}")
        return {'id': sample, 'positions': positions, 'genotype': genotype}

    def genotypes(self, genotypes: Union[str, None]) -> list:
        """
            Return list of genotyes if given
        """
        if genotypes == None:
            self.logger.info("No genotypes provided by user")
            return None
        else:
            self.logger.info(f"Genotypes provided by user: {genotypes.split(' ')}")
            return genotypes.split(' ')

    def get_genotype(self, gt: str, label: bool) -> str:
        """
            Returns genotype from string as either tuple or verbose string
        """
        if not label:
            if gt == '0/0':
                return (0, 0)
            elif gt == '0/1':
                return (0, 1)
            elif gt == '1/1':
                return (1, 1)
            else:
                msg = "Unknown genotype"
                self.logger.error(msg)
                raise ValueError(msg)
        else:
            if gt == '0/0':
                return 'reference'
            elif gt == '0/1':
                return 'heterozygous'
            elif gt == '1/1':
                return 'homozygous'
            else:
                return 'all'

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
            self.logger.error(msg)
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

        for pos in positions:
            for rec in self.vcfFile.fetch(str(chrom), pos -1, pos):
                ad = rec.samples[sample]['AD']
                if ad and sum(ad) > 0:  # Avoid division by zero
                    baf = ad[1] / sum(ad)  # Alt / (Ref + Alt) 
                    allele_fractions.append(baf)
                    variant_positions.append(pos)

        return (allele_fractions, variant_positions)