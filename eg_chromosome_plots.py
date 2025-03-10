import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import glob, math
from matplotlib.backends.backend_pdf import PdfPages

def read_sample_files(cnvs_file, read_depth_file):
    cnvs = pd.read_csv(cnvs_file, sep='\t', header=None, names=["chrom", "start", "end", "val", "type"])
    read_depth = pd.read_csv(read_depth_file, sep='\t', header=None, names=["chrom", "bin", "dosage", "stdev", "unnorm_dosage"])
    read_depth['stdev_pos'] = read_depth['stdev'] + 1
    read_depth['stdev_neg'] = 1 - read_depth['stdev']
    return cnvs, read_depth

def limit_to_chr(read_depth, chr):
    chr_read_depth = read_depth[read_depth["chrom"] == chr]
    return chr_read_depth

def plot_sample(plot_cyto, cytobands, cnvs, chr_read_depth, chr, ax):
    ax.set_ylim([0, 2])
    ax.set_yticks([0, 0.5, 1, 1.5, 2])
    ax.set_yticklabels(["0", "0.5", "1", "1.5", "2"])
    ax.grid(False)
    # plot cytobands
    if plot_cyto:
        for _, row in cytobands.iterrows():
            if row['chrom'] == chr:
                ax.axvspan(row['start'], row['end'], facecolor=row['value'])
    # plot CNVs
    for _, row in cnvs.iterrows():
        if row["chrom"] == chr:
            # deletions
            if row['type'] == 1:
                ax.axvspan(row['start'], row['end'], ymin=0.5, ymax=1, facecolor='hotpink')
            # duplications
            elif row['type'] == 3:
                ax.axvspan(row['start'], row['end'], ymin=0, ymax=0.5, facecolor='deepskyblue')
    # plot read depth
    ax.plot(chr_read_depth['bin'],chr_read_depth['dosage'], color='black')
    ax.fill_between(chr_read_depth['bin'], chr_read_depth['stdev_neg'],chr_read_depth['stdev_pos'], color='springgreen')
    return ax

def format_cytobands(cytobands):
    # read and format hg38_cytoBand.txt
    cytobands['chrom'] = cytobands['chrom'].astype(str).str.replace("chr", "")
    cytobands = cytobands.astype({"chrom": str})
    # extract the centromeres to plot on the overview page
    df_acen = cytobands[cytobands['value'] == "acen"]
    # make centromeres grey
    df_acen['value'] = df_acen['value'].str.replace("acen", "#E8E8E8")
    # now remove acen (centromeres), stalk (short arms of acrocentric chromosomes) and gvar (heterochromatin - pericentric or telomeric) regions
    cytobands = cytobands[~cytobands['value'].str.contains("acen|stalk|gvar", na=False)]
    cytobands['value'] = cytobands['value'].str.replace("gneg", "#F8F8F8")
    cytobands = cytobands.replace({'value' : { 'gpos25' : '#E8E8E8', 'gpos50' : '#E8E8E8', 'gpos75' : '#E8E8E8', 'gpos100' : '#E8E8E8'}})
    return cytobands, df_acen

def plot_report(chrs, cnv_chrs, samples, sample_cnvs, sample_read_depths, df_cytobands, df_acens, pdf_name):
    with PdfPages(pdf_name) as pdf:
        n_rows, n_cols = 6, 4  # 6 rows, 4 columns
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 12), layout="compressed", sharey=True)  # Adjust figure size for readability
        axes = axes.flatten()  # Flatten for easy iteration
        fig.suptitle('proband_id')
        axes[0].set_ylim([0, 2])
        axes[0].set_yticks(np.arange(0, 2, 0.25))
        for idx, c in enumerate(chrs):
            # get chromosome size to use as x-axis max
            chr_max = math.ceil(df_cytobands[df_cytobands['chrom'] == c]['end'].max().item()/10_000_000)*10_000_000
            # plt.subplots_adjust(wspace=0.4,hspace=0.4)
            ax = axes[idx]
            cnvs = sample_cnvs[0]
            read_depth = sample_read_depths[0]
            chr_read_depth = limit_to_chr(read_depth, c)
            ax = plot_sample(True, df_acens, cnvs, chr_read_depth, c, ax)
            ax.set_title("Chromosome " + c)
            del_patch = mpatches.Patch(color='hotpink', label='Deletion')
            dup_patch = mpatches.Patch(color='deepskyblue', label='Duplication')
            noise_patch = mpatches.Patch(color='springgreen', label='Noise')
        plt.figlegend(handles=[del_patch, dup_patch, noise_patch], ncols=3, loc='upper right')
        plt.tight_layout()
        pdf.savefig()
        plt.close()
        for c in cnv_chrs:
            # get chromosome size to use as x-axis max
            chr_max = math.ceil(df_cytobands[df_cytobands['chrom'] == c]['end'].max().item()/10_000_000)*10_000_000
            # set up sub-plots.
            num_subplots = len(samples)
            fig, axs = plt.subplots(num_subplots, figsize=(16, 6), sharex=True)
            # plt.subplots_adjust(wspace=0.4,hspace=0.4)
            fig.suptitle('Chromosome ' + c)
            # x-axis is shared
            axs[0].set_xlim([0, chr_max])
            axs[0].set_xticks(np.arange(0, chr_max, 10000000))
            for idx, s in enumerate(samples):
                ax = axs[idx]
                cnvs = sample_cnvs[idx]
                read_depth = sample_read_depths[idx]
                chr_read_depth = limit_to_chr(read_depth, c)
                ax = plot_sample(True, df_cytobands, cnvs, chr_read_depth, c, ax)
                ax.set_title(s, loc="left")
            del_patch = mpatches.Patch(color='hotpink', label='Deletion')
            dup_patch = mpatches.Patch(color='deepskyblue', label='Duplication')
            noise_patch = mpatches.Patch(color='springgreen', label='Noise')
            plt.figlegend(handles=[del_patch, dup_patch, noise_patch], ncols=3, loc='upper right')
            plt.tight_layout()
            pdf.savefig()
            plt.close()

def main():
    chrs = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
    cnv_chrs=["1", "2", "18"]
    samples = ["","",""]
    family = ""
    # read in cytobands file and reformat
    cytobands = pd.read_csv("hg38_cytoBand.txt", sep='\t', header=None, names=["chrom", "start", "end", "band", "value"])
    df_cytobands, df_acens = format_cytobands(cytobands)
    # read in sample data
    # initialise lists to store sample dataframes
    cnvs_list = []
    read_depth_list = []
    for s in samples:
        cnvs_file = glob.glob(s + "_*cnvs")[0]
        read_depth_file = glob.glob(s + "_*readDepth")[0]
        print(cnvs_file, read_depth_file)
        cnvs, read_depth = read_sample_files(cnvs_file, read_depth_file)
        cnvs_list.append(cnvs)
        read_depth_list.append(read_depth)
    pdf_name = family + "_cnvs_xyax" + ".pdf"
    plot_report(chrs, cnv_chrs, samples, cnvs_list, read_depth_list, df_cytobands, df_acens, pdf_name)

if __name__=="__main__":
    main()