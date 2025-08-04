import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from natsort import index_natsorted
import glob
from datetime import datetime
from collections import OrderedDict

pd.options.mode.chained_assignment = None  # default='warn'

class Sample_Dosage():
    ''' Create dataframes and plot axes for dosage on a per sample level across the whole genome '''

    def __init__(self, readDepthFile, cnvsFile, cytobandsFile, sample, family, noiseCutoff):
        self.readDepthFile = readDepthFile
        self.cnvsFile = cnvsFile
        self.cytobandsFile = cytobandsFile
        self.sample = sample
        self.family = family
        self.noiseCutoff = noiseCutoff
        self.chrs = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]

        # load data
        self.read_depth = self.load_read_depth()
        self.cytobands, self.acens = self.format_cytobands()
        self.cnvs = self.load_cnvs()
        self.regions = self.load_regions()

        # process data
        self.grouped_read_depth = self.process_read_depth()
        # self.ideogram_ax = self.plot_ideogram_ax()

    def load_read_depth(self):
        read_depth = pd.read_csv(self.readDepthFile, sep='\t', header=None, names=["chrom", "bin_start", "bin_end", "dosage", "stdev", "unnorm_dosage", "del_phred", "dup_phred"])
        # calculate 1 SD either way
        read_depth['stdev_pos'] = read_depth['stdev'] + 1
        read_depth['stdev_neg'] = 1 - read_depth['stdev']
        # remove decoy, unlocalised, and unplaced contigs
        read_depth = read_depth[~read_depth.chrom.str.contains("random|decoy|Un|^M", regex=True)]
        # sort by chromosome
        read_depth = read_depth.sort_values(by="chrom", key=lambda x: np.argsort(index_natsorted(read_depth["chrom"])))
        return read_depth

    def load_cnvs(self):
        cnvs_df = pd.read_csv(self.cnvsFile, sep='\t', index_col=False, names=['chrom','start','end','type','bins','total_width_bins','phred','phred_by_bin','dosage','sample'])
        cnvs_df['size'] = cnvs_df['end'] - cnvs_df['start']
        cnvs_df['chrom'] = cnvs_df['chrom'].astype('str')
        return cnvs_df
    
    def load_regions(self):
        regions_file = 'web_ClinGen_region_curation_list_GRCh38_20250425.tsv'
        regions_df = pd.read_csv(regions_file, sep="\t", index_col=False)
        return regions_df

    def format_cytobands(self):
        # read and format hg38_cytoBand.txt
        cytobands = pd.read_csv(self.cytobandsFile, sep='\t', header=None, names=["chrom", "start", "end", "band", "value"])
        cytobands['chrom'] = cytobands['chrom'].astype(str).str.replace("chr", "")
        cytobands = cytobands.astype({"chrom": str})
        # extract the centromeres to plot on the overview page
        df_acen = cytobands[cytobands['value'] == "acen"]
        # make centromeres grey
        df_acen['value'] = df_acen['value'].str.replace("acen", "#E8E8E8")
        # now remove acen (centromeres), stalk (short arms of acrocentric chromosomes) and gvar (heterochromatin - pericentric or telomeric) regions from the main cytobands file
        cytobands = cytobands[~cytobands['value'].str.contains("acen|stalk|gvar", na=False)]
        cytobands['value'] = cytobands['value'].str.replace("gneg", "#F8F8F8")
        cytobands = cytobands.replace({'value' : { 'gpos25' : '#E8E8E8', 'gpos50' : '#E8E8E8', 'gpos75' : '#E8E8E8', 'gpos100' : '#E8E8E8'}})
        return cytobands, df_acen

    @staticmethod
    def get_genome_coordinates(dataframe, col_name):
        chr_pos = {}
        for name,group in dataframe.groupby('chrom', sort=False):
            chr_pos[name]=group[col_name].max().item()
        chr_pos_ordered = OrderedDict(chr_pos)
        # get the cumulative position across the genome
        for k,v in chr_pos_ordered.items():
            if k == '1':
                chr_pos_ordered[k] = chr_pos_ordered[k]
            else:
                chr_pos_ordered[k] = chr_pos_ordered[k] + Sample_Dosage.previous_value(chr_pos_ordered, k)
        return chr_pos_ordered

    @staticmethod
    def apply_genome_coordinates(row, ordered_dict, col_name):
        if row.get('chrom') == '1':
            genome_coord = row.get(col_name)
        else:
            prev_max = Sample_Dosage.previous_value(ordered_dict, row.get('chrom'))
            genome_coord = row.get(col_name) + prev_max
        return genome_coord

    @staticmethod
    def previous_value(dictionary, current_key):
        # Get the list of keys from the OrderedDict
        keys = list(dictionary.keys())
        # Get an index of the current key and offset it by -1
        index = keys.index(current_key) - 1
        # return the previous key's value
        return dictionary[keys[index]]

    def remove_centromeres(self, read_depth, centromeres, chrs):
        for c in chrs:
            acen_start = centromeres[centromeres['chrom'] == c]['start'].min().item()
            acen_end = centromeres[centromeres['chrom'] == c]['end'].max().item()
            read_depth = read_depth[ ((read_depth['chrom'] == c) & ((read_depth['bin_start'] < acen_start) | (read_depth['bin_start'] >= acen_end))) | 
                                                ((read_depth['chrom'] != c))]
        return read_depth

    def limit_noise(self, read_depth, noise_cutoff):
        high_noise = read_depth[read_depth['stdev'] >= noise_cutoff]
        clean = read_depth[read_depth['stdev'] < noise_cutoff]
        return high_noise, clean

    @staticmethod
    def group_and_get_cumulative(read_depth, genome_coord_dict, col_name):
        # get the cumulative coordinate over the whole genome
        read_depth['genome_coordinate'] = read_depth.apply(lambda row: Sample_Dosage.apply_genome_coordinates(row, genome_coord_dict, col_name), axis=1)
        # add an index for sequential bin, irrespective of chromosome
        # sorted_df['ind'] = range(len(sorted_df))
        # group by chromosome for plotting
        grouped_df = read_depth.groupby('chrom', sort=False)
        return grouped_df

    def process_read_depth(self):
        # get the cumulative position across the genome. for compatibility with BAF coordinates
        cumulative_genome_position = self.get_genome_coordinates(self.read_depth, 'bin_end')
        # remove centromeres
        depth_sans_acens = self.remove_centromeres(self.read_depth, self.acens, self.chrs)
        # remove bins that have noise above our cut-off for CNV calling
        high_noise, low_noise = self.limit_noise(depth_sans_acens, self.noiseCutoff)
        # sort and group
        grouped_read_depth = self.group_and_get_cumulative(low_noise, cumulative_genome_position, 'bin_end')
        return grouped_read_depth

    def plot_ideogram_ax(self, ax, capture):
        colours = ['#4477AA', '#EE6677']
        # dot size
        dot_size = 5 if capture == 'genome' else 7
        # lables for chromosomes
        x_labels = []
        x_labels_pos = []
        ymax = self.grouped_read_depth.max()
        print(ymax)
        for num, (name, group) in enumerate(self.grouped_read_depth):
            # plot noise first
            ax.fill_between(group['genome_coordinate'], group['stdev_neg'],group['stdev_pos'], color='#CCBB44')
            # plot each group (chromosome) and colour using the color pallete
            group.plot(kind='scatter', x='genome_coordinate', y='dosage',color=colours[num % len(colours)], ax=ax, legend=None, s=dot_size, rasterized=True)
            x_labels.append(name)
            x_labels_pos.append((group['genome_coordinate'].iloc[-1] - (group['genome_coordinate'].iloc[-1] - group['genome_coordinate'].iloc[0])/2))
        ax.set_xlim([0, len(self.grouped_read_depth)])
        ax.set_xticks(x_labels_pos)
        ax.set_xticklabels(x_labels)
        ax.set_title(self.family + ', noise cut-off = ' + str(self.noiseCutoff))
        return ax

    def plot_chr_ideogram_ax(self, ax, dosage_data, capture, chr, outerind):
        # dot size
        dot_size = 5 if capture == 'genome' else 8
        ax.fill_between(dosage_data['bin_end'], dosage_data['stdev_neg'],dosage_data['stdev_pos'], color='#CCBB44')
        ax.scatter(dosage_data['bin_end'], dosage_data['dosage'], s=dot_size, rasterized=True)
        ax.set_ylabel("Dosage")
        ymin, ymax = self.get_ylims(dosage_data['dosage'])
        ax.set_ylim(ymin, ymax)  # dosage values range between 0 and 2
        ax.set_xlim(left=0)
        # ax.set_title(self.sample, loc='left')
        # ax.set_title(self.sample)
        ax.text(-0.05, 0, self.sample, rotation='horizontal',
                ha='right', va='center', transform=ax.transAxes,
                fontsize=12)
        self.cnvs_track(ax, chr)
        if outerind == 0:
            self.isca_track(ax, chr)
        return ax

    def get_ylims(self, dosage):
        ''' Fix y-axis limits from 0-2 unless there is a dosage >2.2,
         then fix between the 0 and the max dosage '''
        ymin = -0.1
        if max(dosage) > 2.2:
            ymax = round(max(dosage), 2) + 0.1
        else:
            ymax = 2.2
        return ymin, ymax

    def cnvs_track(self, ax, chr):
        track_y0 = 1.025
        track_h = 0.1
        track_box = mpatches.Rectangle((0, 1), width=1, height=0.1, facecolor='none', edgecolor='none', transform=ax.transAxes, clip_on=False)
        ax.add_patch(track_box)
        # transform y coordinates to plot regions using data coordinates
        patch_y0, patch_h = self.transform_y_point(ax, track_y0, track_h)
        cnv_count = 0
        for _, row in self.cnvs.iterrows():
            if row["chrom"] == chr:
                cnv_count += 1
                patch_col = 'red' if row['type'] == 'Deletion' else 'deepskyblue'
                cnv_rect = mpatches.Rectangle((row['start'], patch_y0), width=row['size'], height=patch_h, color=patch_col, clip_on=False)
                ax.add_patch(cnv_rect)
        text_y0 = track_y0 + (track_h / 2)
        cnvs_txt = f'CNVs ({cnv_count})'
        ax.text(-0.002, text_y0, cnvs_txt, fontsize=12, verticalalignment='center', horizontalalignment='right', transform=ax.transAxes)
        return ax

    def isca_track(self, ax, chr):
        track_y0 = 1.15
        track_h = 0.05
        track_box = mpatches.Rectangle((0, track_y0), width=1, height=track_h, facecolor='none', edgecolor='none', transform=ax.transAxes, clip_on=False)
        ax.add_patch(track_box)
        # transform y coordinates to plot regions using data coordinates
        patch_y0, patch_h = self.transform_y_point(ax, track_y0, track_h)
        region_count = 0
        for _, row in self.regions.iterrows():
            if row["chrom"] == chr:
                region_count += 1
                region_rect = mpatches.Rectangle((row['start'], patch_y0), width=row['end']-row['start'], height=patch_h, color='limegreen', clip_on=False)
                ax.add_patch(region_rect)
        text_y0 = track_y0 + (track_h / 2)
        isca_txt = f'ISCA regions ({region_count})'
        ax.text(-0.002, text_y0, isca_txt, fontsize=12, verticalalignment='center', horizontalalignment='right', transform=ax.transAxes)
        return ax

    def transform_y_point(self, ax, relative_y0, relative_height):
        ''' For adding tracks to figures. Given a relative y-axis (y0) coordinate, and a relative height convert to data coordinates '''
        # get the display coordinate for the axes-relative y position, then use invert to get the data coordinate
        disp_y0 = ax.transAxes.transform((0, relative_y0))[1]  # Just need y
        data_y0 = ax.transData.inverted().transform((0, disp_y0))[1]
        # Same for height: get display coordinate of track_y0 + track_h, subtract
        disp_y1 = ax.transAxes.transform((0, relative_y0 + relative_height))[1]
        data_y1 = ax.transData.inverted().transform((0, disp_y1))[1]
        patch_h = data_y1 - data_y0
        return data_y0, patch_h

# def main():
#     # Create instance of class
#     Dosage = Sample_Dosage(
#         readDepthFile='WGS_EX2410666_22KLFTLT3.coverageBinner.50000.data',
#         cnvsFile='cnvs_WGS_EX2410666_22KLFTLT3.50000',
#         cytobandsFile='hg38_cytoBand.txt',
#         sample='WGS_EX2410666',
#         family='F08782',
#         noiseCutoff=0.3
#         )

# if __name__ == '__main__':
#     main()