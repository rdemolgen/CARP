import argparse
import matplotlib.pyplot as plt
import pysam
import re

class Allele_Fraction():
    '''This class include functions to identify variants and generate allele fractions'''

    def __init__(self, vcfFile, location, no_filtering, min_qual, min_dp, min_gq,
        samples=None, genotypes=None):
        self.vcfFile = vcfFile
        self.vcf = self.load_vcf()
        self.chrom, self.start, self.end = self.get_location(location)
        self.samples = self.get_samples(samples)
        self.genotypes = self.genotypes(genotypes)
        self.no_filtering = no_filtering
        self.min_qual = min_qual # Phred score 99.9%
        self.min_dp = min_dp
        self.min_gq = min_gq # Phred score 99.0%

        # Run automatic plot if no genotypes are given
        if self.genotypes is None:
            if len(self.samples) == 1:
                print(f"Automatically generating single sample BAF plots for {self.samples}.")
                self.run_single_plots()
            else:
                print(f"Automatically generating single sample and join BAF plots for {self.samples}.")
                self.run_single_plots()
                self.run_joint_call_plots()
        else:
            if len(self.samples) == 1 and len(self.genotypes) == 1:
                print(f"Generating single BAF plot for {self.samples} with {self.genotypes} genotype.")
                self.single_sample_baf(samples, self.chrom, self.start, self.end, self.genotypes[0])
            elif len(self.samples) == len(self.genotypes):
                print(f"Generating single and joint call BAF plots for {self.samples} with {self.genotypes} genotypes.")
                self.run_single_plots()
                self.joint_call_bcf(self.samples, self.chrom, self.start, self.end, self.genotypes)
            else:
                print("Number of samples and genotypes do not match")

    def get_samples(self, samples):
        '''Return list of samples, assuming proband is the first sample id'''
        if samples == None:
            return list(self.vcf.header.samples)
        else:
            return samples.split(' ')
    
    def genotypes(self, genotypes):
        '''Return list of genotyes if given'''
        if genotypes == None:
            return None
        else:
            return genotypes.split(' ')
    
    def load_vcf(self):
        '''Create a pysam object from vcf file'''
        return pysam.VariantFile(self.vcfFile)

    def get_location(self, location):
        '''Return chrom, start, end from location string <chr:start-end>'''
        try:
            match = re.match(r"^(1?[0-9]|2[0-2]|X|Y|MT)(?::(\d+)-(\d+))?$", location)
            return match.groups()
        except AttributeError as e:
            print(f"Genomic location incorrectly formatted, '{location}'")
            return None, None, None
    
    def get_genotype(self, gt, label):
        '''Returns genotype from string as either tuple or verbose string'''
        if not label:
            if gt == '0/0':
                return (0, 0)
            elif gt == '0/1':
                return (0, 1)
            elif gt == '1/1':
                return (1, 1)
            else:
                return None
        else:
            if gt == '0/0':
                return 'reference'
            elif gt == '0/1':
                return 'heterozygous'
            elif gt == '1/1':
                return 'homozygous'
            else:
                return 'all'
        
    def get_chr_len(self, chrom):
        '''Return chromosome length from vcf header'''
        chrom_dict = {str(i): i - 1 for i in range(1, 23)}
        chrom_dict.update({'X': 22, 'Y': 23, 'MT': 24})

        try:
            return self.vcf.header.contigs[chrom_dict[str(chrom)]].length
        except KeyError:
            print('Chromosome not recognised.')
            return None

    def get_variants(self, sample, chrom, start=None, end=None, genotype=None):
        '''Return variant position based on given location criteria, genotype and quality thresholds'''
        gt = self.get_genotype(genotype, False)
        positions = []
        if start is None: start = 1
        if end is None: end = self.get_chr_len(chrom)

        for rec in self.vcf.fetch(str(chrom), int(start), int(end)):
            sample_data = rec.samples[sample]

            if sample_data["GT"] != gt and gt != None:
                continue
            # Accept only "PASS" or "." variants
            if not self.no_filtering:
                if "PASS" not in rec.filter.keys() and "." not in rec.filter.keys():
                    continue
            if rec.qual is not None and rec.qual < self.min_qual:
                # print("insufficient quality")
                continue
            if "DP" in sample_data and sample_data["DP"] < self.min_dp:
                # print("insufficient depth")
                continue
            try:
                if "GQ" in sample_data and sample_data["GQ"] < self.min_gq:
                    # print("insufficient GQ")
                    continue
            except Exception as e:
                continue

            positions.append(rec.pos)

        return {'id': sample, 'positions': positions, 'genotype': genotype}

    def intersect_sample_variants(self, sample_pos):
        '''Find matching variant positions between samples'''
        intersect_pos = set(sample_pos[0]['positions'])

        for sample in sample_pos[1:]:
            intersect_pos &= set(sample['positions'])
        
        intersect_pos = list(intersect_pos)
        intersect_pos.sort()

        return intersect_pos

    def calc_baf(self, sample, chrom, positions):
        '''Returns baf and variant position'''
        allele_fractions = [] # Y axis
        variant_positions = [] # X axis

        for pos in positions:
            for rec in self.vcf.fetch(str(chrom), pos -1, pos):
                ad = rec.samples[sample]['AD']
                if ad and sum(ad) > 0:  # Avoid division by zero
                    baf = ad[1] / sum(ad)  # Alt / (Ref + Alt) 
                    allele_fractions.append(baf)
                    variant_positions.append(pos)
        
        return (allele_fractions, variant_positions)
    
    def plot_baf(self, plot_data, sample, chrom, start=1, end=None, genotype=None):
        # Extract plot data
        allele_fractions, variant_positions = plot_data

        samp_geno_label = []
        if isinstance(sample, list) and isinstance(genotype, list):
            for s, g in zip(sample, genotype):
                samp_geno_label.append(s)
                samp_geno_label.append(self.get_genotype(g, True)[:3])
        else:
            samp_geno_label.append(sample)
            samp_geno_label.append(self.get_genotype(genotype, True)[:3])

        samp_geno_label = '_'.join(samp_geno_label)

        if (start == 1 or start is None) and end is None:
            location = f"chr{chrom}"
        else:
            location = f"chr{chrom}.{start}-{end}"

        # Plotting the Allele Fraction Scatter Plot
        plt.figure(figsize=(10, 6))
        plt.scatter(variant_positions, allele_fractions, s=5, alpha=0.3, color='purple', edgecolors='purple')
        plt.axhline(0.667, color='black', linestyle='dashed', linewidth=1)
        plt.axhline(0.5, color='black', linestyle='dashed', linewidth=1)
        plt.axhline(0.337, color='black', linestyle='dashed', linewidth=1)

        # Set custom Y-axis ticks
        custom_ticks = [0, 0.25, 0.337, 0.5, 0.667, 0.75, 1]
        plt.yticks(custom_ticks, labels=[str(tick) for tick in custom_ticks])

        plt.xlabel("Genomic Position (bp)")
        plt.ylabel("Allele Fraction (Alt / (Ref + Alt))")
        plt.title(f"BAF Plot: {samp_geno_label} {location}")
        plt.ylim(-0.1, 1.1)  # AF values range between 0 and 1
        # plt.legend()

        # Save the plot as an image file
        output_filename = f"{samp_geno_label}_{location}_BAF.png"
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')  # Save as high-quality PNG
        plt.close()  # Close the plot to prevent it from displaying in some environments

        print(f"Plot saved as {output_filename}")
    
    def single_sample_baf(self, sample, chrom, start=1, end=None, genotype=None):
        '''Generate single sample baf plot'''
        var_pos = self.get_variants(sample, str(chrom), start, end, genotype)
        plot_data = self.calc_baf(sample, str(chrom), var_pos['positions'])
        self.plot_baf(plot_data, sample, str(chrom), start, end, genotype)

    def joint_call_bcf(self, samples, chrom, start=1, end=None, genotypes=None):
        '''Generate joint call baf based on individual sample genotypes'''
        samples_var_pos = []
        for sample, genotype in zip(samples, genotypes):
            samples_var_pos.append(self.get_variants(sample, chrom, start, end, genotype))
        shared_positions = self.intersect_sample_variants(samples_var_pos)
        plot_data = self.calc_baf(samples[0], str(chrom), shared_positions)
        self.plot_baf(plot_data, samples, str(chrom), start, end, genotypes)

    def run_single_plots(self):
        '''Automatically generate different genotype plots for each sample'''
        genotypes = ['0/1', '1/1', 'all']

        for sample in self.samples:
            for genotype in genotypes:
                # gt = self.get_genotype(genotype, False)
                self.single_sample_baf(sample, self.chrom, self.start, self.end, genotype)
    
    def run_joint_call_plots(self):
        '''Automatically generate joint call baf plots based on the number of samples given for typical scenarios'''
        genotypes_combos = self.genotype_combinations(len(self.samples))
        
        for genotypes_combo in genotypes_combos:
            samples_var_pos = []
            for sample, genotype in zip(self.samples, genotypes_combo):
                samples_var_pos.append(self.get_variants(sample, str(self.chrom), self.start, self.end, genotype))
            shared_positions = self.intersect_sample_variants(samples_var_pos)
            plot_data = self.calc_baf(self.samples[0], str(self.chrom), shared_positions)
            self.plot_baf(plot_data, self.samples, str(self.chrom), self.start, self.end, genotypes_combo)
    
    def genotype_combinations(self, no_samples):
        '''Return genotype combinations to automatically generate standard plots'''
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
            return None

        
def main():
    parser = argparse.ArgumentParser(description="")
    # Arguments
    parser.add_argument('-v', '--vcfFile', type=str, required=True, help="VCF file")
    parser.add_argument('-l', '--location', type=str, required=True, help="Genomic location either chr or chr:start-end")
    parser.add_argument('-s', '--samples', type=str, required=False, help="List of sample ids, starting with proband separated by spaces")
    parser.add_argument('-g', '--genotypes', type=str, required=False, help="List of genotypes matching the order of sample ids")
    parser.add_argument('-f', '--no_filtering', action='store_true', required=False, help="Accept varaints with other non-PASS filters, default=False")
    parser.add_argument('-vq', '--min_qual', type=int, required=False, default=30, help="Min variant quality score, default=30")
    parser.add_argument('-dp', '--min_dp', type=int, required=False, default=10, help="Min variant read depth, default=10")
    parser.add_argument('-gq', '--min_gq', type=int, required=False, default=20, help="Min genotype quality score, default=20")

    # Parse args
    args = parser.parse_args()

    # Create instance of class
    AF = Allele_Fraction(
        vcfFile=args.vcfFile,
        location=args.location,
        samples=args.samples,
        genotypes=args.genotypes,
        no_filtering=args.no_filtering,
        min_qual=args.min_qual,
        min_dp=args.min_dp,
        min_gq=args.min_gq
        )

if __name__ == '__main__':
    main()