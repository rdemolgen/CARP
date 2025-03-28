from allele_fraction import Allele_Fraction
from savvycnv_dosage import Sample_Dosage
import pandas as pd
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from datetime import datetime
import math, argparse, re, os, sys

def make_genome_ideogram(AF, genome_baf, Dosage):
    '''  
        Plot dosage,baf ideograms for the whole genome and return the plt object
        Allows you to dump in a PDF file or save as an image 
    '''
    n_rows, n_cols = 2, 1
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(24, 12), layout="compressed", sharex=True)
    axes[0].set_ylim([0, 2])
    axes[0].set_yticks(np.arange(0, 2, 0.25))
    Dosage.plot_ideogram_ax(axes[0])
    # Set custom Y-axis ticks
    custom_ticks = [0, 0.25, 0.337, 0.5, 0.667, 0.75, 1]
    axes[1].set_yticks(custom_ticks, labels=[str(tick) for tick in custom_ticks])
    AF.plot_baf_ideogram(axes[1], genome_baf)
    plt.tight_layout()
    return plt

def make_chromosome_ideograms(AF, sample_genome_baf, dosage_dict, chr):
    '''  
        Plot dosage,baf ideograms for each chromosome and return the plt object
        Allows you to dump in a PDF file or save as an image 
    '''
    num_subfigs = len(AF.samples)
    fig = plt.figure(figsize=(24, 12))
    fig.suptitle(f"Chromosome {chr}", x=0.125)
    if num_subfigs > 1:
        subfigs = fig.subfigures(num_subfigs, 1, wspace=0, hspace=-0.1)
        for outerind, subfig in enumerate(subfigs.flat):
            sample = AF.samples[outerind]
            axs = subfig.subplots(2, 1, sharex=True)
            dosage_dict[sample].plot_chr_ideogram_ax(axs[0], dosage_dict[sample].grouped_read_depth.get_group(chr))
            try:
                AF.plot_baf(sample_genome_baf[sample].get_group(chr), sample, chr, ax=axs[1], genotype='all')
            except AttributeError:
                AF.plot_baf(sample_genome_baf[sample], sample, chr, ax=axs[1], genotype='all')
    else:
        sample = AF.samples[0]
        axs = fig.subplots(2, 1, sharex=True)
        dosage_dict[sample].plot_chr_ideogram_ax(axs[0], dosage_dict[sample].grouped_read_depth.get_group(chr))
        try:
            AF.plot_baf(sample_genome_baf[sample].get_group(chr), sample, chr, ax=axs[1], genotype='all')
        except AttributeError:
            AF.plot_baf(sample_genome_baf[sample], sample, chr, ax=axs[1], genotype='all')
    plt.subplots_adjust(wspace=0, hspace=0)
    return plt

def pdf_report(pdf_name, AF, sample_genome_baf, dosage_dict, chrs, proband):
    '''  
        Generate a PDF report containing dosage,baf ideograms for the whole genome (proband)
        and for each chromosome (all family members)
    '''
    with PdfPages(pdf_name) as pdf:
        plt = make_genome_ideogram(AF, sample_genome_baf[proband], dosage_dict[proband])
        pdf.savefig()
        plt.close()
        for c in chrs:
            plt = make_chromosome_ideograms(AF, sample_genome_baf, dosage_dict, c)
            pdf.savefig()
            plt.close()

# move the cytobands stuff to it's own section to avoid repeating with every Class object

def find_file(regex_pattern):
    '''  
        Find a file based on the provied glob_string, e.g. 'WGS_EX1234567*data'
        Error if the file can't be found or if there are duplicates of the file
    '''
    reg_files = [f for f in os.listdir('../test_files/') if re.match(regex_pattern, f)]

    if len(reg_files) > 1:
        print("Error: More than one file found matching the regex pattern " + regex_pattern)
        pass
    elif len(reg_files) == 0:
        print("Error: No file found matching the glob string" + regex_pattern)
        pass
    else:
        file_path = os.path.join('../test_files/', reg_files[0])
        return file_path

def main():
    parser = argparse.ArgumentParser(description="")
    # Arguments
    parser.add_argument('-v', '--vcfFile', type=str, required=False, help="VCF file")
    parser.add_argument('-l', '--location', type=str, required=True, help="Genomic location either chr or chr:start-end")
    parser.add_argument('-s', '--samples', type=str, required=False, help="List of sample ids, starting with proband separated by spaces")
    parser.add_argument('-g', '--genotypes', type=str, required=False, help="List of genotypes matching the order of sample ids")
    parser.add_argument('-o', '--outDir', type=str, required=False, help="Output directory for plots.")
    parser.add_argument('-f', '--no_filtering', action='store_true', required=False, help="Accept varaints with other non-PASS filters (QD>2,MQ>40), default=False")
    parser.add_argument('-vq', '--min_qual', type=int, required=False, default=30, help="Min variant quality score, default=30")
    parser.add_argument('-dp', '--min_dp', type=int, required=False, default=10, help="Min variant read depth, default=10")
    parser.add_argument('-gq', '--min_gq', type=int, required=False, default=20, help="Min genotype quality score, default=20")
    parser.add_argument('-mq', '--min_mq', type=int, required=False, default=40, help="Min mapping quality score, default=40")
    parser.add_argument('-qd', '--min_qd', type=int, required=False, default=2, help="Min qual-by-depth score, default=2")
    parser.add_argument('-fam', '--family', type=str, required=True, help="family")
    parser.add_argument('-p', '--proband_id', type=str, required=True, help="proband ID")
    parser.add_argument('-m', '--mode', type=str, required=True, help="Valid values: baf, ideogram, dosage")

    # Parse args
    args = parser.parse_args()

    if args.mode != 'dosage':
        print(f'Running in {args.mode} mode')

        print(f"Finding family VCF files for {args.proband_id}")
        vcf_pat = rf'^(?!.*gnomad_filtered).*{re.escape(args.proband_id)}.*\.vcf\.gz$'
        vcfFile = find_file(vcf_pat)
        print(f"Loaded {vcfFile}")

        if args.mode == 'baf' and str(args.location) == 'all':
            sys.exit('Please enter chromosome or coordinates when running in BAF mode (e.g. )')
        # populate an Allele_Fraction object in baf and ideogram mode
        # baf mode will auto-generate the baf plots. ideogram mode will just populate the class
        AF = Allele_Fraction(
            vcfFile=vcfFile,
            location=args.location,
            samples=args.samples,
            genotypes=args.genotypes,
            no_filtering=args.no_filtering,
            min_qual=args.min_qual,
            min_dp=args.min_dp,
            min_gq=args.min_gq,
            min_mq=args.min_mq,
            min_qd=args.min_qd,
            outDir=args.outDir,
            mode=args.mode
            )

    if args.mode == 'ideogram':
        # Generate ideogram plots
        sample_genome_baf = {}
        dosage_dict = {}
        for sample in AF.samples:
            if args.location == 'all':
                # whole genome ideograms
                print(f"Loading genome-wide BAF for {sample}")
                baf = AF.get_genome_wide_baf(sample)
            else:
                # single chromosome ideograms
                print(f"Loading Chromosome {args.location} BAF for {sample}")
                baf = AF.get_plot_data(sample, args.location, 1, None, 'all')

            print(f"Finding CNV files for {sample}")
            cnvs_pat = rf'^cnvs_{re.escape(sample)}.*\.50000$'
            readDepthFile = find_file(cnvs_pat)
            data_pat = rf'^{re.escape(sample)}.*\.50000.data$'
            cnvsFile = find_file(data_pat)
            print(f"Loaded {cnvsFile} and {readDepthFile}")

            Dosage = Sample_Dosage(
                readDepthFile=readDepthFile,
                cnvsFile=cnvsFile,
                cytobandsFile='../hg38_cytoBand.txt',
                sample=sample,
                family=args.family,
                noiseCutoff=0.3
                )

            dosage_dict[sample] = Dosage
            sample_genome_baf[sample] = baf

        now=datetime.now()
        date=now.strftime('%Y-%m-%d')

        if args.location == 'all':
            chrs = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
            pdf_name = f"{args.family}_genome_ideogram_{date}.pdf"
        else:
            chrs = [args.location]
            pdf_name = f"{args.family}_chr{args.location}_ideogram_{date}.pdf"

        print("Generating whole genome ideogram PDF report")
        pdf_report(pdf_name, AF, sample_genome_baf, dosage_dict, chrs, args.proband_id)

    elif args.mode == 'dosage':
        print(f'Running in {args.mode} mode')

if __name__ == '__main__':
    main()