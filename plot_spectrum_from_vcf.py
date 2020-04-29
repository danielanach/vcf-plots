import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rc
from nucleic import DNA
from nucleic import SnvSpectrum, Notation
import numpy as np
import os.path as path
import palettable
import pandas as pd
from pyfaidx import Fasta
from pysam import VariantFile
import scipy.stats as stats
import seaborn as sns
from substitution_spectrum import SubstitutionSpectrum

def plot_spectrum_by_vaf(vcf, sample_name, n_bins=10, title=''):
    """Plots the substitution spectrum by VAF .

    Parameters
    ----------
    vcf : str
        VCF file path.

    sample_name : str
        Name of sample.

    n_bins : int
        Number of bins to split variants by VAF into. Default = 10.

    title : str
        Title to put on the plot.
    """

    vcf = VariantFile(vcf)
    c, bins = np.histogram([0],bins=n_bins,range=(0,1))
    spectrum = SubstitutionSpectrum(samples=bins)

    mut_bins = {k:[] for k in spectrum.mut_count_dict.keys()}

    for var in vcf:

        chrom = var.chrom
        start = var.start
        ref = var.ref

        is_indel = False

        if len(ref) > 1:
            is_indel = True
            continue

        for i in range(len(var.alts)):
            alt = var.alts[i]
            if len(alt) > 1:
                is_indel = True
                continue
            af = var.info['AF'][i]

            snv =  DNA(ref).to(alt).with_pyrimidine_ref()
            snv_str = str(snv).split('[')[1].split(']')[0]

            if not is_indel:
                mut_bins[snv_str].append(float(af))

    for mut in spectrum.mut_count_dict.keys():
        spectrum.mut_count_dict[mut] = list(np.histogram(mut_bins[mut],bins=n_bins,range=(0,1))[0])

    bin_names = ['{0:.1f}'.format(b) for b in bins]
    plot_spectrum_with_hist(spectrum,bin_names[1:],title=title)

def plot_spectrum_by_sample(vcf_lst, sample_list, title=''):
    """Plots the substitution spectrum by VAF .

    Parameters
    ----------
    vcf : str
        VCF file path.

    sample_list : list(str)
        list of sample names.

    n_bins : int
        number of bins to split variants by VAF into
    """

    spectrum = SubstitutionSpectrum(samples = sample_list)

    for v in vcf_lst:

        vcf = VariantFile(v)
        samples = sample_list

        snvs = []

        for var in vcf:

            chrom = var.chrom
            start = var.start
            ref = var.ref

            is_indel = False

            if len(ref) > 1:
                is_indel = True
                continue

            for alt in var.alts:

                if len(alt) > 1:
                    is_indel = True
                    continue

                snv =  DNA(ref).to(alt).with_pyrimidine_ref()
                snv_str = str(snv).split('[')[1].split(']')[0]

                if not is_indel:
                    snvs.append(snv_str)

        for mut in spectrum.mut_count_dict.keys():
            spectrum.mut_count_dict[mut].append(snvs.count(mut))

    plot_spectrum_with_hist(spectrum,sample_list,title=title)


def plot_spectrum(spectrum, samples,  title):
    """Plots the substitution spectrum .

    Parameters
    ----------
    spectrum : SubstitutionSpectrum
        Contains substituion type and fraction contained.
    """

    mut_dict = spectrum.normalized_muts()
    names = samples

    sample_counts = spectrum.counts_by_category()

    colors = palettable.colorbrewer.diverging.Spectral_6.hex_colors[::-1]
    sns.set_style('white')

    r = list(range(len(samples)))
    barWidth = 0.85

    sns.set_style('white')
    fig = plt.figure(figsize=(4,4),dpi=250,facecolor='w')
    ax = fig.add_subplot(1, 1, 1)

    ax.bar(r,
            mut_dict['T→G'],
            color=colors[5],
            edgecolor='white',
            width=barWidth,
            label='T→G')

    ax.bar(r,
            mut_dict['T→C'],
            bottom=mut_dict['T→G'],
            color=colors[4],
            edgecolor='white',
            width=barWidth,
            label='T→C')

    ax.bar(r,
            mut_dict['T→A'],
            bottom=[i+j for i,j in zip(mut_dict['T→C'],
                                       mut_dict['T→G'])],
            color=colors[3],
            edgecolor='white',
            width=barWidth,
            label='T→A')

    ax.bar(r,
            mut_dict['C→T'],
            bottom=[i+j+k for i,j,k in zip(mut_dict['T→A'],
                                           mut_dict['T→C'],
                                           mut_dict['T→G'])],
            color=colors[2],
            edgecolor='white',
            width=barWidth,
            label='C→T')

    ax.bar(r,
            mut_dict['C→G'],
            bottom=[i+j+k+l for i,j,k,l in zip(mut_dict['C→T'],
                                               mut_dict['T→A'],
                                               mut_dict['T→C'],
                                               mut_dict['T→G'])],
            color=colors[1],
            edgecolor='white',
            width=barWidth,
            label='C→G')

    ax.bar(r,
            mut_dict['C→A'],
            bottom=[i+j+k+l+m for i,j,k,l,m in zip(mut_dict['C→G'],
                                                   mut_dict['C→T'],
                                                   mut_dict['T→A'],
                                                   mut_dict['T→C'],
                                                   mut_dict['T→G'])],
            color=colors[0],
            edgecolor='white',
            width=barWidth,
            label='C→A')

    plt.xticks(r, names,rotation=90)

    ax.set_xlabel('Sample',fontweight='bold',fontsize=18)
    ax.set_ylabel('Fraction of mutations',fontweight='bold',fontsize=18)

    handles, labels = ax.get_legend_handles_labels()

    ax.legend(handles = handles[::-1], labels=labels[::-1], loc='upper left', bbox_to_anchor=(1,0.7), ncol=1)
    ax.set_title(title, fontsize=10)
    return spectrum, ax


def plot_spectrum_with_hist(spectrum, samples,  title):
    """Plots the substitution spectrum .

    Parameters
    ----------
    spectrum : SubstitutionSpectrum
        Contains substituion type and fraction contained.
    """

    mut_dict = spectrum.normalized_muts()
    names = samples

    sample_counts = spectrum.counts_by_category()

    colors = palettable.colorbrewer.diverging.Spectral_6.hex_colors[::-1]

    r = list(range(len(samples)))
    barWidth = 0.85
    mpl.style.use('fivethirtyeight')
    sns.set_style('white')
    fig = plt.figure(figsize=(6,8),dpi=250,facecolor='w')

    ax_0 = fig.add_subplot(2, 1, 1)

    ax_0.bar(r,sample_counts,color='grey')
    ax_0.set_xticklabels([])
    ax_0.axes.get_xaxis().set_visible(False)
    ax_0.set_ylabel('# of mutations',fontweight='bold',fontsize=18)
    ax_0.set_title(title, fontweight='bold', fontsize=20)

    ax = fig.add_subplot(2, 1, 2,sharex=ax_0)

    ax.bar(r,
            mut_dict['T→G'],
            color=colors[5],
            edgecolor='white',
            width=barWidth,
            label='T→G')

    ax.bar(r,
            mut_dict['T→C'],
            bottom=mut_dict['T→G'],
            color=colors[4],
            edgecolor='white',
            width=barWidth,
            label='T→C')

    ax.bar(r,
            mut_dict['T→A'],
            bottom=[i+j for i,j in zip(mut_dict['T→C'],
                                       mut_dict['T→G'])],
            color=colors[3],
            edgecolor='white',
            width=barWidth,
            label='T→A')

    ax.bar(r,
            mut_dict['C→T'],
            bottom=[i+j+k for i,j,k in zip(mut_dict['T→A'],
                                           mut_dict['T→C'],
                                           mut_dict['T→G'])],
            color=colors[2],
            edgecolor='white',
            width=barWidth,
            label='C→T')

    ax.bar(r,
            mut_dict['C→G'],
            bottom=[i+j+k+l for i,j,k,l in zip(mut_dict['C→T'],
                                               mut_dict['T→A'],
                                               mut_dict['T→C'],
                                               mut_dict['T→G'])],
            color=colors[1],
            edgecolor='white',
            width=barWidth,
            label='C→G')

    ax.bar(r,
            mut_dict['C→A'],
            bottom=[i+j+k+l+m for i,j,k,l,m in zip(mut_dict['C→G'],
                                                   mut_dict['C→T'],
                                                   mut_dict['T→A'],
                                                   mut_dict['T→C'],
                                                   mut_dict['T→G'])],
            color=colors[0],
            edgecolor='white',
            width=barWidth,
            label='C→A')

    plt.xticks(r, names,rotation=90)

    ax.set_xlabel('Sample',fontweight='bold',fontsize=18)
    ax.set_ylabel('Fraction of mutations',fontweight='bold',fontsize=18)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles = handles[::-1], labels=labels[::-1], loc='upper left', bbox_to_anchor=(1,0.7), ncol=1)

    plt.subplots_adjust(wspace=0, hspace=0.05)
    plt.show()

    return spectrum, ax
