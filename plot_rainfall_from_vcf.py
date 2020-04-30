import numpy as np
import pandas as pd
from pysam import VariantFile
import matplotlib.pyplot as plt
import seaborn as sns

from nucleic import DNA
from nucleic import SnvSpectrum, Notation
import palettable

class ContigCoordinates:
    def __init__(self,vcf_in):
        self.contig_dct = self.get_contig_length(vcf_in)
        self.length_adj = self.get_contig_coord_offset()

    def get_contig_length(self, vcf_in, skip_alt_contig=True, skip_mito=True):
        '''Get contig lengths from VCF header.

        Parameters
        ----------
        vcf_in
        A pysam VariantFile object.

        skip_alt_contig
        Ignore alternative contigs.

        skip_mito
        Ignore mitochondrial contigs.
        '''
        contig_dct = {}
        for c in vcf_in.header.records:
          c = str(c)
          if 'contig' in c:
            if skip_alt_contig:
              if '_' in c:
                continue
            if skip_mito:
              if 'chrM' in c:
                continue

            contig = c.split('=')[2].split(',')[0]
            length = int(c.split('=')[3].split('>')[0])
            contig_dct[contig] = length

        return contig_dct
    def get_contig_coord_offset(self):
        '''Get offset coordinates for contigs.

        Parameters
        ----------
        vcf_in
        A pysam VariantFile object.

        skip_alt_contig
        Ignore alternative contigs.

        skip_mito
        Ignore mitochondrial contigs.
        '''
        length_adj = {}
        total = 0
        for c in self.contig_dct.keys():
            length_adj[c] = total
            total += self.contig_dct[c]
        return length_adj

def vcf_to_dist_table(vcf):

    vcf_in = VariantFile(vcf)

    coord_obj = ContigCoordinates(vcf_in)
    contig_dct = coord_obj.contig_dct
    length_adj = coord_obj.length_adj

    chroms = []
    starts = []
    coords = []
    muts = []

    for rec in vcf_in:
        chrom = rec.chrom
        if chrom not in length_adj.keys():
            continue
        for alt in rec.alts:
            snv =  DNA(rec.ref).to(alt).with_pyrimidine_ref()
            snv_str = str(snv).split('[')[1].split(']')[0]
            if len(snv_str) > 3:
                continue
            chroms.append(chrom)
            starts.append(rec.start)
            coords.append(rec.start+length_adj[chrom])
            muts.append(str(snv_str))

    mut_table = pd.DataFrame([chroms,
                      muts,
                      starts,
                      coords]).transpose()
    mut_table.columns = ['chr','mut','start','coord']

    shifted_list = [0] + mut_table['coord'].to_list()[:-1]
    mut_table['dist'] = mut_table['coord'] - shifted_list
    mut_table['maxed_dist'] = np.clip(mut_table['dist'],0,5000)

    return mut_table, length_adj, contig_dct

def plot_rainfall_from_vcf(vcf,sample=None):

    if not sample:
        sample = vcf

    mut_table,length_adj, contig_dct = vcf_to_dist_table(vcf)
    mut_table = mut_table.sort_values(by='mut')

    colors = palettable.colorbrewer.diverging.Spectral_6.hex_colors[::-1]
    subs = ['C→A', 'C→G', 'C→T', 'T→A', 'T→C', 'T→G']
    color_dct = dict(zip(subs,colors))

    plt.figure(figsize=(14,4),dpi=300)

    sns.scatterplot(x='coord',
                y='maxed_dist',
                data=mut_table,
                palette=color_dct,
                edgecolor='dimgray',
                s=200,
                hue='mut')

    plt.xticks(list(length_adj.values()),
               list(length_adj.keys()),
               rotation=45,
               fontsize=12,
               fontweight='bold')


    for c in length_adj:
        plt.axvline(x=length_adj[c],
                    color='darkgrey')

    plt.xlim(0,length_adj['chrY'] + contig_dct['chrY'])

    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()

    plt.legend(bbox_to_anchor=(1.0,0.85),
                handles=handles[1:],
                labels=labels[1:],
                prop={'size': 16,},
                markerscale=2,
                frameon=False)

    plt.ylabel('Genomic Distance',
               fontweight='bold',
               fontsize=16)

    plt.xlabel('Genomic Coordinate',
               fontweight='bold',
               fontsize=16)

    plt.title(sample,
              fontsize=20,
              fontweight='bold')

    return plt.gca()
