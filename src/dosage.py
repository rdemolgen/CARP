import numpy as np
import pandas as pd

from collections import OrderedDict
from natsort import index_natsorted
from pathlib import Path
from typing import Tuple

class Dosage:
    """
        Create dataframes and plot axes for dosage on a per sample level across the whole genome
    """

    def __init__(self, sample, prefix, readDepthFile, cnvsFile, cytobandsFile, iscaFile, chrs, noiseCutoff=0.3):
        # self.logger = logger
        self.sample = sample
        self.prefix = prefix
        self.noiseCutoff = noiseCutoff
        self.chrs = chrs
        # load data
        self.readDepth = self.load_read_depth(readDepthFile)
        self.cytobands, self.acens = self.format_cytobands(cytobandsFile)
        self.cnvs = self.load_cnvs(cnvsFile)
        self.regions = self.load_regions(iscaFile)

        # process data
        self.grouped_read_depth = self.process_read_depth()

    def load_read_depth(self, readDepthFile: Path) -> pd.DataFrame:
        read_depth = pd.read_csv(readDepthFile, sep='\t', header=None, names=["chrom", "bin_start", "bin_end", "dosage", "stdev", "unnorm_dosage", "del_phred", "dup_phred"])
        # calculate 1 SD either way
        read_depth['stdev_pos'] = read_depth['stdev'] + 1
        read_depth['stdev_neg'] = 1 - read_depth['stdev']
        # remove decoy, unlocalised, and unplaced contigs
        read_depth = read_depth[~read_depth["chrom"].astype(str).str.contains("random|decoy|Un|^M", regex=True)]
        # sort by chromosome
        read_depth = read_depth.sort_values(by="chrom", key=lambda x: np.argsort(index_natsorted(read_depth["chrom"])))
        return read_depth

    def load_cnvs(self, cnvsFile: Path) -> pd.DataFrame:
        cnvs_df = pd.read_csv(cnvsFile, sep='\t', index_col=False, names=['chrom','start','end','type','bins','total_width_bins','phred','phred_by_bin','dosage','sample'])
        cnvs_df['size'] = cnvs_df['end'] - cnvs_df['start']
        cnvs_df['chrom'] = cnvs_df['chrom'].astype('str')
        return cnvs_df
    
    def load_regions(self, iscaFile: Path) -> pd.DataFrame:
        regions_df = pd.read_csv(iscaFile, sep="\t", index_col=False)
        return regions_df

    def format_cytobands(self, cytobandsFile: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
        # read and format hg38_cytoBand.txt
        cytobands = pd.read_csv(cytobandsFile, sep='\t', header=None, names=["chrom", "start", "end", "band", "value"])
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
                chr_pos_ordered[k] = chr_pos_ordered[k] + Dosage.previous_value(chr_pos_ordered, k)
        return chr_pos_ordered

    @staticmethod
    def apply_genome_coordinates(row, ordered_dict, col_name):
        if row.get('chrom') == '1':
            genome_coord = row.get(col_name)
        else:
            prev_max = Dosage.previous_value(ordered_dict, row.get('chrom'))
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
        read_depth['genome_coordinate'] = read_depth.apply(lambda row: Dosage.apply_genome_coordinates(row, genome_coord_dict, col_name), axis=1)
        # add an index for sequential bin, irrespective of chromosome
        # sorted_df['ind'] = range(len(sorted_df))
        # group by chromosome for plotting
        grouped_df = read_depth.groupby('chrom', sort=False)
        return grouped_df

    def process_read_depth(self):
        # get the cumulative position across the genome. for compatibility with BAF coordinates
        cumulative_genome_position = self.get_genome_coordinates(self.readDepth, 'bin_end')
        # remove centromeres
        depth_sans_acens = self.remove_centromeres(self.readDepth, self.acens, self.chrs)
        # remove bins that have noise above our cut-off for CNV calling
        high_noise, low_noise = self.limit_noise(depth_sans_acens, self.noiseCutoff)
        # sort and group
        grouped_read_depth = self.group_and_get_cumulative(low_noise, cumulative_genome_position, 'bin_end')
        return grouped_read_depth

