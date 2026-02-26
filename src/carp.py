import argparse, datetime, re

from pathlib import Path

from src.baf import Baf
from src.plots import Plots
from src.utility import Utility

class Carp:

    def __init__(self, args):
        self.now = datetime.datetime.now().strftime("%Y%m%d-%H%M%S") 
        self.mode = args.mode
        self.location = None
        self.proband_id = args.proband_id
        self.prefix = args.prefix
        
        self.inDir = self.utility.verify_dir(args.inDir)
        self.outDir = self.utility.verify_dir(args.outDir)

        self.utility = Utility(self.outDir, self.prefix, self.proband_id, self.now)
        self.logger = self.utility.logger
        
        # do we want to have different min defaults for WES and WGS data
        self.filters = {
            "no_filter": args.no_filter,
            "min_qual": args.min_qual,
            "min_dp": args.min_dp,
            "min_gq": args.min_gq,
            "min_mq": args.min_mq,
            "min_qd": args.min_qd,
        }
        

        self.capture = self.utility.set_capture(self.proband_id)
        self.vcfFile = self.utility.find_file(self.inDir, rf"^(?!.*gnomad_filtered).*{re.escape(self.proband_id)}.*\.vcf\.gz$")
        # self.samples = self.

        self.baf = Baf(self.logger, self.filters, self.vcfFile)



    def run(self):
        """main function to direct this application"""

        if self.mode == "baf":
            pass
        elif self.mode == "ideogram":
            pass
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
    

if __name__ == "__main__":
    main()