
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from matplotlib.axes import Axes
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path
from typing import Optional, Tuple
from .utility import Utility

class Plots:

    def __init__(self):
        pass

    @staticmethod
    def plot_baf_ideogram(ax: Axes, grouped_depth: dict, capture: str) -> Axes:
        """
            Return sub plots of baf per chromosome
        """
        colours = ['#4477AA', '#EE6677']
        # dot size
        dot_size = 0.002 if capture == 'genome' else 4
        # lables for chromosomes
        for num, (name, group) in enumerate(grouped_depth):
            # plot each group (chromosome) and colour using the color pallete
            group.plot(kind='scatter', x='genome_coordinate', y='allele_fraction',color=colours[num % len(colours)], ax=ax, legend=None, s=dot_size, rasterized=True)
        
        return ax

    @staticmethod
    def make_genome_ideogram(genome_baf: dict, dosage: dict, capture: str) -> plt:
        '''  
            Plot dosage,baf ideograms for the whole genome and return the plt object
            Allows you to dump in a PDF file or save as an image 
        '''
        n_rows, n_cols = 2, 1
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(24, 12), layout="compressed", sharex=True)
        dosage_data = dosage.grouped_read_depth.max()
        ymin, ymax = Plots.get_ylims(dosage_data["dosage"])
        axes[0].set_ylim(ymin, ymax)
        axes[0].set_yticks(np.arange(ymin, ymax, 0.25))
        Plots.plot_ideogram_ax(axes[0], capture, dosage)
        # Set custom Y-axis ticks
        custom_ticks = [0, 0.25, 0.337, 0.5, 0.667, 0.75, 1]
        axes[1].set_yticks(custom_ticks, labels=[str(tick) for tick in custom_ticks])
        Plots.plot_baf_ideogram(axes[1], genome_baf, capture)
        plt.tight_layout()

        return plt

    @staticmethod
    def make_chromosome_ideograms(samples: list, sample_genome_baf: dict, dosage_dict: dict, chr: str, capture: str) -> plt:
        '''  
            Plot dosage,baf ideograms for each chromosome and return the plt object
            Allows you to dump in a PDF file or save as an image 
        '''
        num_subfigs = len(samples)
        fig = plt.figure(figsize=(24, 12))
        fig.suptitle(f"Chromosome {chr}", x=0.03, y=0.99)
        # fig.patches.extend([mpatches.Rectangle((0.125,0.98),0.865,0.01,
        #                               facecolor='none', edgecolor='black',
        #                               transform=fig.transFigure, figure=fig)])
        if num_subfigs > 1:
            subfigs = fig.subfigures(num_subfigs, 1, wspace=0, hspace=-0.1)
            for outerind, subfig in enumerate(subfigs.flat):
                sample = samples[outerind]
                # check if the data is just for the one chr, or genoem data grouped by chr
                if isinstance(sample_genome_baf[sample], tuple):
                    chr_plot_data = sample_genome_baf[sample]
                else:
                    try:
                        chr_plot_data = sample_genome_baf[sample].get_group(chr)
                    except KeyError:
                        chr_plot_data = pd.DataFrame(columns=['chrom','position','allele_fraction','genome_coordinate'])
                axs = subfig.subplots(2, 1, sharex=True)
                Plots.plot_chr_ideogram_ax(sample, axs[0], dosage_dict[sample].grouped_read_depth.get_group(chr), capture, chr, outerind, dosage_dict[sample].cnvs, dosage_dict[sample].regions)
                Plots.plot_baf(chr_plot_data, sample, chr, capture, ax=axs[1], genotype='all')
        else:
            sample = samples[0]
            if isinstance(sample_genome_baf[sample], tuple):
                chr_plot_data = sample_genome_baf[sample]            
            else:
                chr_plot_data = sample_genome_baf[sample].get_group(chr)
            axs = fig.subplots(2, 1, sharex=True)
            Plots.plot_chr_ideogram_ax(sample, axs[0], dosage_dict[sample].grouped_read_depth.get_group(chr), capture, chr, 0, dosage_dict[sample].cnvs, dosage_dict[sample].regions)
            # dosage_dict[proband].cnvs_track(axs[0], chr)
            Plots.plot_baf(chr_plot_data, sample, chr, capture, ax=axs[1], genotype='all')

        # plt.subplots_adjust(wspace=0, hspace=0.1, top=0.2, left=0.08, right=0.1)
        plt.subplots_adjust(wspace=0, hspace=0, bottom=0.12, right=0.99)
        return plt

    @staticmethod
    def plot_ideogram_ax(ax: Axes, capture: str, dosage: dict) -> Axes:
        """
            Plot all subplots
        """
        colours = ['#4477AA', '#EE6677']
        # dot size
        dot_size = 5 if capture == 'genome' else 7
        # lables for chromosomes
        x_labels = []
        x_labels_pos = []
        ymax = dosage.grouped_read_depth.max()
        for num, (name, group) in enumerate(dosage.grouped_read_depth):
            # plot noise first
            ax.fill_between(group['genome_coordinate'], group['stdev_neg'],group['stdev_pos'], color='#CCBB44')
            # plot each group (chromosome) and colour using the color pallete
            group.plot(kind='scatter', x='genome_coordinate', y='dosage',color=colours[num % len(colours)], ax=ax, legend=None, s=dot_size, rasterized=True)
            x_labels.append(name)
            x_labels_pos.append((group['genome_coordinate'].iloc[-1] - (group['genome_coordinate'].iloc[-1] - group['genome_coordinate'].iloc[0])/2))
        ax.set_xlim([0, len(dosage.grouped_read_depth)])
        ax.set_xticks(x_labels_pos)
        ax.set_xticklabels(x_labels)
        ax.set_title(dosage.prefix + ', noise cut-off = ' + str(dosage.noiseCutoff))
        return ax

    @staticmethod
    def plot_chr_ideogram_ax(sample: str, ax: Axes, dosage_data: dict, capture: str, chr: str, outerind: int, cnvs: pd.DataFrame, regions: pd.DataFrame) -> Axes:
        """
            Plot chromosome ideogram subplots
        """
        # dot size
        dot_size = 5 if capture == 'genome' else 8
        ax.fill_between(dosage_data['bin_end'], dosage_data['stdev_neg'],dosage_data['stdev_pos'], color='#CCBB44')
        ax.scatter(dosage_data['bin_end'], dosage_data['dosage'], s=dot_size, rasterized=True)
        ax.set_ylabel("Dosage")
        ymin, ymax = Plots.get_ylims(dosage_data['dosage'])
        ax.set_ylim(ymin, ymax)  # dosage values range between 0 and 2
        ax.set_xlim(left=0)
        # ax.set_title(self.sample, loc='left')
        # ax.set_title(self.sample)
        ax.text(-0.05, 0, sample, rotation='horizontal',
                ha='right', va='center', transform=ax.transAxes,
                fontsize=12)
        Plots.cnvs_track(ax, chr, cnvs)
        if outerind == 0:
            Plots.isca_track(ax, chr, regions)
        return ax

    @staticmethod
    def pdf_report(pdf_name: str, samples: list, sample_genome_baf: dict, dosage_dict: dict, chrs: list, proband: str, capture: str):
        '''  
            Generate a PDF report containing dosage,baf ideograms for the whole genome (proband)
            and for each chromosome (all family members)
        '''
        with PdfPages(pdf_name) as pdf:
            plt = Plots.make_genome_ideogram(sample_genome_baf[proband], dosage_dict[proband], capture)
            pdf.savefig()
            plt.close()
            for c in chrs:
                plt = Plots.make_chromosome_ideograms(samples, sample_genome_baf, dosage_dict, c, capture)
                pdf.savefig()
                plt.close()

    @staticmethod
    def plot_baf(plot_data: tuple, sample: list, chrom: str, capture: Optional[str]="genome", ax=None, start: Optional[int]=1, end: Optional[int]=None, genotype: Optional[str]=None, outDir: Optional[Path]=None) -> plt:
        """
            Generate a single baf plot as a .png file
        """
        # Extract plot data
        try:
            allele_fractions, variant_positions = plot_data
        except (TypeError, ValueError, KeyError):
            variant_positions = plot_data['position']
            allele_fractions = plot_data['allele_fraction']

        samp_geno_label = []
        if isinstance(sample, list) and isinstance(genotype, list):
            for s, g in zip(sample, genotype):
                samp_geno_label.append(s)
                samp_geno_label.append(Utility.get_genotype(g, True)[:3])
        else:
            samp_geno_label.append(sample)
            samp_geno_label.append(Utility.get_genotype(genotype, True)[:3])
        samp_geno_label = '_'.join(samp_geno_label)

        if (start == 1 or start is None) and end is None:
            location = f"chr{chrom}"
        else:
            location = f"chr{chrom}.{start}-{end}"

        # dot size
        dot_size = 0.9 if capture == 'genome' else 10
        # Plotting the Allele Fraction Scatter Plot
        # Set custom Y-axis ticks
        custom_ticks = [0, 0.25, 0.337, 0.5, 0.667, 0.75, 1]
        if ax == None:
            plt.figure(figsize=(10, 6))
            pltAx = plt
            pltAx.yticks(custom_ticks, labels=[str(tick) for tick in custom_ticks])
            pltAx.xlabel("Genomic Position (bp)")
            pltAx.ylabel("Allele Fraction (Alt / (Ref + Alt))")
            pltAx.ylim(-0.1, 1.1)  # AF values range between 0 and 1
        else:
            pltAx = ax
            pltAx.set_yticks(custom_ticks, labels=[str(tick) for tick in custom_ticks])
            # pltAx.set_xlabel("Genomic Position (bp)")
            pltAx.set_ylabel("Allele Fraction")
            pltAx.set_ylim(-0.1, 1.1)  # AF values range between 0 and 1
            pltAx.set_xlim(left=0)
        pltAx.scatter(variant_positions, allele_fractions, s=dot_size, color='purple', edgecolors='purple', rasterized=True)
        pltAx.axhline(0.667, color='black', linestyle='dashed', linewidth=1)
        pltAx.axhline(0.5, color='black', linestyle='dashed', linewidth=1)
        pltAx.axhline(0.337, color='black', linestyle='dashed', linewidth=1)
        # pltAx.title(f"BAF Plot: {samp_geno_label} {location}")
        # plt.legend()
        if ax == None:
            # Save the plot as an image file
            output_filename = outDir / f"{samp_geno_label}_{location}_BAF.png"
            pltAx.savefig(output_filename, dpi=300, bbox_inches='tight')  # Save as high-quality PNG
            pltAx.close()  # Close the plot to prevent it from displaying in some environments
            print(f"Plot saved as {output_filename}")
        else:
            return pltAx

    @staticmethod
    def get_ylims(dosage: float) -> Tuple[float, float]:
        ''' Fix y-axis limits from 0-2 unless there is a dosage >2.2,
        then fix between the 0 and the max dosage '''
        ymin = -0.1
        if max(dosage) > 2.2:
            ymax = round(max(dosage), 2) + 0.1
        else:
            ymax = 2.2
        return ymin, ymax

    @staticmethod
    def cnvs_track(ax: Axes, chr: str, cnvs: pd.DataFrame) -> Axes:
        """
            Annotate plots with cnv track
        """
        track_y0 = 1.025
        track_h = 0.1
        track_box = mpatches.Rectangle((0, 1), width=1, height=0.1, facecolor='none', edgecolor='none', transform=ax.transAxes, clip_on=False)
        ax.add_patch(track_box)
        # transform y coordinates to plot regions using data coordinates
        patch_y0, patch_h = Plots.transform_y_point(ax, track_y0, track_h)
        cnv_count = 0
        for _, row in cnvs.iterrows():
            if row["chrom"] == chr:
                cnv_count += 1
                patch_col = 'red' if row['type'] == 'Deletion' else 'deepskyblue'
                cnv_rect = mpatches.Rectangle((row['start'], patch_y0), width=row['size'], height=patch_h, color=patch_col, clip_on=False)
                ax.add_patch(cnv_rect)
        text_y0 = track_y0 + (track_h / 2)
        cnvs_txt = f'CNVs ({cnv_count})'
        ax.text(-0.002, text_y0, cnvs_txt, fontsize=12, verticalalignment='center', horizontalalignment='right', transform=ax.transAxes)
        return ax

    @staticmethod
    def isca_track(ax: Axes, chr: str, regions: pd.DataFrame) -> Axes:
        """
            Annotate plots with cnv track
        """
        track_y0 = 1.15
        track_h = 0.05
        track_box = mpatches.Rectangle((0, track_y0), width=1, height=track_h, facecolor='none', edgecolor='none', transform=ax.transAxes, clip_on=False)
        ax.add_patch(track_box)
        # transform y coordinates to plot regions using data coordinates
        patch_y0, patch_h = Plots.transform_y_point(ax, track_y0, track_h)
        region_count = 0
        for _, row in regions.iterrows():
            if row["chrom"] == chr:
                region_count += 1
                region_rect = mpatches.Rectangle((row['start'], patch_y0), width=row['end']-row['start'], height=patch_h, color='limegreen', clip_on=False)
                ax.add_patch(region_rect)
        text_y0 = track_y0 + (track_h / 2)
        isca_txt = f'ISCA regions ({region_count})'
        ax.text(-0.002, text_y0, isca_txt, fontsize=12, verticalalignment='center', horizontalalignment='right', transform=ax.transAxes)
        return ax

    @staticmethod
    def transform_y_point(ax: Axes, relative_y0: float, relative_height: float) -> tuple[float, float]:
        """
            For adding tracks to figures. Given a relative y-axis (y0) coordinate, and a relative height convert to data coordinates
        """
        # get the display coordinate for the axes-relative y position, then use invert to get the data coordinate
        disp_y0 = ax.transAxes.transform((0, relative_y0))[1]  # Just need y
        data_y0 = ax.transData.inverted().transform((0, disp_y0))[1]
        # Same for height: get display coordinate of track_y0 + track_h, subtract
        disp_y1 = ax.transAxes.transform((0, relative_y0 + relative_height))[1]
        data_y1 = ax.transData.inverted().transform((0, disp_y1))[1]
        patch_h = data_y1 - data_y0
        return data_y0, patch_h