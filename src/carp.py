import argparse, datetime, re

from .baf import Baf
from .dosage import Dosage
from .plots import Plots
from .utility import Utility

class Carp:

    def __init__(self, args):
        self.now = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        self.date = datetime.datetime.now().strftime("%Y-%m-%d")
        self.chrs = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
        self.mode = args.mode
        self.proband_id = args.proband_id
        self.prefix = args.prefix
        self.utility = Utility(self.prefix, self.proband_id, self.now)
        
        self.inDir = self.utility.verify_dir(args.inDir)
        self.outDir = self.utility.verify_dir(args.outDir)

        self.logger = self.utility.logger
        
        # do we want to have different min defaults for WES and WGS data?
        self.filters = {
            "no_filtering": args.no_filtering,
            "min_qual": args.min_qual,
            "min_dp": args.min_dp,
            "min_gq": args.min_gq,
            "min_mq": args.min_mq,
            "min_qd": args.min_qd,
        }
        
        self.capture = self.utility.set_capture(self.proband_id)
        vcf_regex = (
            rf"^(?!.*gnomad_filtered)"
            rf"(?:WGS_EX\d{{7}}-)*{re.escape(self.proband_id)}"
            rf"(?:(?:WGS_EX\d{{7}}|TwEx\d*_EX\d{{7}})-)*{re.escape(self.proband_id)}"
            rf"(?:-(?:WGS_EX\d{{7}}|TwEx\d*_EX\d{{7}}))*\.vcf\.gz$"
        )
        self.vcfFile = self.utility.find_file(self.inDir, vcf_regex)
        self.baf = Baf(self.logger, self.filters, self.vcfFile, self.chrs)
        self.samples = Baf.get_samples(self.baf.vcf, args.samples)
        self.genotypes = Baf.genotypes(args.genotypes)
        self.location = args.location
        self.chrom, self.start, self.end = Baf.get_location(args.location)
        # dosage
        self.cytoFile = args.cyto
        self.iscaFile = args.isca

    def run(self):
        """main function to direct this application"""
        # print filters being used
        for filter, value in self.filters.items():
            print(f"\t{filter}: {value}")
        
        # generate baf plots only
        if self.mode == "baf":
            print("Running baf mode")
            if str(self.location) == "all":
                print("Provide a chromosome (and optional coordinates) when running BAF mode")
                raise SystemExit
            
            if self.genotypes is None:
                if len(self.samples) == 1:
                    print(f"Automatically generating single sample BAF plots for {self.samples}.")
                    self.baf.run_single_plots(self.samples, self.chrom, self.start, self.end, self.outDir)
                else:
                    print(f"Automatically generating single sample and joint BAF plots for {self.samples}.")
                    self.baf.run_single_plots(self.samples, self.chrom, self.start, self.end, self.outDir)
                    self.baf.run_joint_call_plots(self.proband_id, self.samples, self.chrom, self.start, self.end, self.outDir)
            else:
                if len(self.samples) == 1 and len(self.genotypes) == 1:
                    print(f"Generating single BAF plot for {self.samples} with {self.genotypes} genotype.")
                    self.baf.run_single_plots(self.samples, self.chrom, self.start, self.end, self.outDir, self.genotypes)
                elif len(self.samples) == len(self.genotypes):
                    print(f"Generating joint call BAF plots for {self.samples} with {self.genotypes} genotypes.")
                    self.baf.run_joint_call_plots(self.proband_id, self.samples, self.chrom, self.start, self.end, self.outDir, self.genotypes)
                else:
                    print("Number of samples and genotypes do not match")  
                    raise SystemExit                     
        # generate ideogram plots
        elif self.mode == "ideogram":
            print("Running ideogram mode")
            read_depth_files = {}
            cnvs_files = {}            
            for sample in self.samples:
                print(f"Finding CNV files for {sample}")
                readDepthFile = self.utility.find_file(self.inDir, rf"^{re.escape(sample)}.*\.20000.data$")
                cnvsFile = self.utility.find_file(self.inDir, rf'cnvs_{re.escape(sample)}.*\.20000$')
                read_depth_files[sample] = readDepthFile
                cnvs_files[sample] = cnvsFile
                print(f"Loaded {cnvsFile} and {readDepthFile}")
            
            sample_genome_baf = {}
            dosage_dict = {}
            for sample in self.samples:
                if self.location == 'all':
                    # whole genome ideograms
                    print(f"Loading genome-wide BAF for {sample}")
                    bafData = self.baf.get_genome_wide_baf(sample)
                else:
                    # single chromosome ideograms
                    print(f"Loading Chromosome {self.location} BAF for {sample}")
                    bafData = self.baf.get_plot_data(sample, self.location, 1, None, 'all')

                sample_dosage = Dosage(
                    readDepthFile=read_depth_files[sample],
                    cnvsFile=cnvs_files[sample],
                    cytobandsFile=self.cytoFile,
                    iscaFile=self.iscaFile,
                    sample=sample,
                    prefix=self.prefix,
                    chrs=self.chrs,
                    noiseCutoff=0.3
                )

                dosage_dict[sample] = sample_dosage
                sample_genome_baf[sample] = bafData

            if self.location == "all":
                print("Generating whole genome ideogram PDF report")
                pdf_name = self.outDir / f"{self.prefix}_genome_ideogram_{self.date}.pdf"
                Plots.pdf_report(pdf_name, self.samples, sample_genome_baf, dosage_dict, self.chrs, self.proband_id, self.capture)
            else:
                print(f"Generating ideogram for Chromosome {self.location}")
                ideogramPlot = Plots.make_chromosome_ideograms(self.samples, sample_genome_baf, dosage_dict, self.location, self.capture)
                outname = self.outDir / f"{self.prefix}_chr{self.location}_ideogram_{self.date}.png"
                ideogramPlot.savefig(outname)
                ideogramPlot.close('all')                                     

        elif self.mode == "dosage":
            raise NotImplementedError("Dosage analysis not implemented yet")
        else:
            raise RuntimeError("Unrecognisable mode")

def main():
    parser = argparse.ArgumentParser(description="")
    # Required arguments
    parser.add_argument('-m', '--mode', type=str, required=True, help="Valid values: baf, ideogram, dosage")
    parser.add_argument('-l', '--location', type=str, required=True, help="Genomic location either chr or chr:start-end")
    parser.add_argument('-p', '--proband_id', type=str, required=True, help="proband ID")
    parser.add_argument('--prefix', type=str, required=True, help="Prefix for output file names")   
    # Optional arguments
    ## Input/Output
    parser.add_argument('-i', '--inDir', type=str, required=False, help="Input file location.")
    parser.add_argument('-o', '--outDir', type=str, required=False, help="Output directory for plots.")
    parser.add_argument('-s', '--samples', type=str, required=False, help="List of sample ids, starting with proband separated by spaces")
    parser.add_argument('-g', '--genotypes', type=str, required=False, help="List of genotypes matching the order of sample ids")
    parser.add_argument('-v', '--vcfFile', type=str, required=False, help="VCF file")
    parser.add_argument('--cyto', type=str, required=False, default="resources/hg38_cytoBand.txt", help="Path to cytobands file")
    parser.add_argument('--isca', type=str, required=False, default="resources/web_ClinGen_region_curation_list_GRCh38_20250425.tsv", help="Path to ISCA regions file")
    ## Variant filters
    parser.add_argument('-f', '--no_filtering', action='store_true', required=False, help="Accept varaints with other non-PASS filters (QD>2,MQ>40), default=False")
    parser.add_argument('-vq', '--min_qual', type=int, required=False, default=30, help="Min variant quality score, default=30")
    parser.add_argument('-dp', '--min_dp', type=int, required=False, default=10, help="Min variant read depth, default=10")
    parser.add_argument('-gq', '--min_gq', type=int, required=False, default=20, help="Min genotype quality score, default=20")
    parser.add_argument('-mq', '--min_mq', type=int, required=False, default=40, help="Min mapping quality score, default=40")
    parser.add_argument('-qd', '--min_qd', type=int, required=False, default=2, help="Min qual-by-depth score, default=2")
   
    # Parse args
    args = parser.parse_args()

    # Initialise the CARP
    carp = Carp(args)
    # Run CARP
    carp.run()
    
if __name__ == "__main__":
    main()