# -*- coding: utf-8 -*-
"""
Created on Mon May 10 20:42:42 2021

@author: djross
"""
import os
import glob

import pandas as pd
import numpy as np

from Bio import Align 
from Bio import SeqIO
from Bio import AlignIO

import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from matplotlib.ticker import MultipleLocator

import seaborn as sns
sns.set()
# set global default style:
sns.set_style("white")
sns.set_style("ticks", {'xtick.direction':'in', 'xtick.top':True, 'ytick.direction':'in', 'ytick.right':True, })
#sns.set_style({"axes.labelsize": 20, "xtick.labelsize" : 16, "ytick.labelsize" : 16})

plt.rcParams['axes.labelsize'] = 20
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16

plt.rcParams['legend.fontsize'] = 14
plt.rcParams['legend.edgecolor'] = 'k'

sanger_channels = ["DATA9", "DATA10", "DATA11", "DATA12"]


def mott_trimming_fr(record, trim=0.05):
    """
    This method trims the low-quality reads from the ends of each read
    using Mott's algorithm in both directions and returns the idices
    that indicate the most conservative trimming on each end
    
    Parameters
    ----------
    record : Biopython sequence object to be trimmed
        
    trim : float
        A value between 0 and 1 indicating the error probability thershhold
        used for the sequence trimming

    Returns
    -------
    the indices for the first and last sequence element that should be kept
    (or np.nan if the sequence quality never exceeds the threshold)
        
    """
    
    # Forward direction
    qual = np.array(record.letter_annotations['phred_quality'])
    score_list = trim - (10 ** (qual / -10.0))
    summed_score_list = []
    for x in score_list:
        if len(summed_score_list)>0:
            summed_score = summed_score_list[-1] + x
        else:
            summed_score = x
        if summed_score < 0:
            summed_score = 0
        summed_score_list.append(summed_score)
    summed_list1 = np.array(summed_score_list)
    
    try:
        f_start1 = np.where(summed_list1>0)[0][0]
        f_end1 = np.where(summed_list1==max(summed_list1))[0][0]
    except IndexError:
        f_start1 = np.nan
        f_end1 = np.nan
    
    # Reverse direction
    record_r = record.reverse_complement()
    qual = np.array(record_r.letter_annotations['phred_quality'])
    score_list = trim - (10 ** (qual / -10.0))
    summed_score_list = []
    for x in score_list:
        if len(summed_score_list)>0:
            summed_score = summed_score_list[-1] + x
        else:
            summed_score = x
        if summed_score < 0:
            summed_score = 0
        summed_score_list.append(summed_score)
    summed_list2 = np.array(summed_score_list)
    
    try:
        f_end2 = len(record_r) - 1 - np.where(summed_list2>0)[0][0]
        f_start2 = len(record_r) - 1 - np.where(summed_list2==max(summed_list2))[0][0]
    except IndexError:
        f_start2 = np.nan
        f_end2 = np.nan
    
    if np.any(np.isnan([f_start1, f_end1, f_start2, f_end2])):
        return np.nan, np.nan, np.nan, np.nan
    else:
        return max(f_start1, f_start2), min(f_end1, f_end2)
    
def num_matches(align1):
    """
    This method calculates the number of matched bases and the number of
    mismathced bass in a pairwise alignment

    Parameters
    ----------
    align1 : Bio.Align.PairwiseAlignment
        The alignment

    Returns
    -------
    match_count : int
        number of matched bases.
    mismatch_count : int
        number of mismatched bases.

    """
    align_str = f'{align1}'
    align_str = align_str.split('\n')
    match_str = align_str[1]
    match_count = match_str.count('|')
    mismatch_count = match_str.count('.')
    
    return match_count, mismatch_count
    
def align_sanger(record1, record2, trim=0.01, verbose=True):
    """
    This method performs a pairwise alignment of two very similar reads.
    
    The method first trims the low-quality reads from the ends of each read
    (if trim != None), and then aligns record1 with the reverse complement of
    record 2
    
    Parameters
    ----------
    record1, record2 : Biopython sequence objects to be aligned
        record1 and record2 must have with phred_quality annotation for trimming
        e.g. created from raw Sanger data via 'record1 = SeqIO.read(x, "abi")'
        
    trim : float
        A value between 0 and 1 indicating the error probability thershhold
        used for the sequence trimming. If trim == None, the input sequences
        are not trimmed
    
    verbose : Boolean
        If true, method prints info about the resulting alignments

    Returns
    -------
    the resulting alignments, as a Bio.Align.PairwiseAlignments object
        (an iterator over pairwise alignments)
        
    """
    
    # Aligner settings chosen for matching nearly identical sequences
    #     e.g., foward and reverse Sanger reads, 
    #     or Sanger reads to reference with only a few expected mutations
    aligner = Align.PairwiseAligner()
    aligner.match_score = 5
    aligner.mismatch_score = -9
    aligner.mode = 'global'
    aligner.target_internal_open_gap_score = -12
    aligner.target_internal_extend_gap_score = -3
    aligner.query_internal_open_gap_score = -12
    aligner.query_internal_extend_gap_score = -3
    
    if trim is not None:
        trim_ind1 = mott_trimming_fr(record1, trim=trim)
        trim_ind2 = mott_trimming_fr(record2, trim=trim)
        new_record1 = record1[trim_ind1[0]:trim_ind1[1]+1]
        new_record2 = record2[trim_ind2[0]:trim_ind2[1]+1]
    else:
        new_record1 = record1
        new_record2 = record2
    
    alignments = aligner.align(new_record1.seq, new_record2.seq.reverse_complement())
    if verbose: print(f'{len(alignments)} alignment(s) found with score: {alignments.score}')
    
    
    align1 = alignments[0]
    if verbose: print(f'{align1.aligned}')
    align_str = f'{align1}'
    align_str = align_str.split('\n')
    if len(align1.aligned[0])>1:
        target_str = align_str[0]
        target_str = target_str.strip('-')
        if '--' in target_str:
            if verbose: print('multi-base gap in forward sequence')
            #x = target_str.find('--')
            #print(f'{x}, {target_str}')
        elif '-' in target_str:
            count = target_str.count('-')
            if verbose: print(f'{count} gap(s) in forward sequence')
            #x = target_str.find('-')
            #print(f'{x}, {target_str}')
            
        query_str = align_str[2]
        query_str = query_str.strip('-')
        if '--' in query_str:
            if verbose: print('multi-base gap in reverse sequence')
            #print(f'{query_str}')
        elif '-' in query_str:
            count = query_str.count('-')
            if verbose: print(f'{count} gap(s) in reverse sequence')
            #print(f'{query_str}')
    else:
        if verbose: print('no gaps')
        
    match_str = align_str[1]
    count = match_str.count('|')
    if verbose: print(f'{count} matches in alignment')
    count = match_str.count('.')
    if verbose: print(f'{count} mismatches in alignment')
    #print(align1)
    return alignments

def is_good_sanger(alignment, min_matches=20, max_mismatch_ratio=0.02):
    """
    This method decides whether or not a pairwise alignment should be counted
    as a good result for matched forward and reverse Sanger sequences

    Parameters
    ----------
    alignment : Bio.Align.PairwiseAlignment
        The alignment, resulting from a pairwise alignment of mathced forward 
        and reverse Sanger reads.
    min_matches : int, optional
        The minumum number of matches for a good alignment. The default is 20.
    max_mismatch_ratio : float, optional
        The maximum value of the ratio (number of matches)/(number of mismatches)
        for a good alignment. The default is 0.02.

    Returns
    -------
    Boolean
        True if the alignement is good.

    """
    x, y = num_matches(alignment)
    return (x >= min_matches) and (y/x <= max_mismatch_ratio)


def slice_sanger(sequence, b0=0, b1=None):
    """
    This method slices Sanger data, including the chromatorgram data
    and peak locations for plotting

    Parameters
    ----------
    sequence : Biopython sequence object
        The sequence to be sliced.
    b0 : int, optional
        The first index for the slice. The default is 0.
    b1 : int, optional
        The end index for the slice. The default is None.

    Returns
    -------
    new_seq : Biopython sequence object
        The sliced sequence.

    """
    # If b0 and/or b1 are None, still trim chromatorgram data
    if b1 == None:
        b1 = len(sequence)
    
    # Built-in slicing handles sequenc and quality information
    new_seq = sequence[b0: b1]
    new_seq.annotations["abif_raw"] = {}
    
    # chromatogram data and peak locations needs to be added on
    raw_data = sequence.annotations["abif_raw"]
    peak_locations = np.array(raw_data['PLOC1'])
    
    # left and right edges of each chromatogram peak
    peak_left = peak_locations
    peak_left = (peak_left[:-1] + peak_left[1:])/2
    peak_right = np.append( peak_left, peak_locations[-1] + 1/2*(peak_left[-1] - peak_left[-2]) )
    peak_left = np.insert( peak_left, 0, peak_locations[0] - 1/2*(peak_left[1] - peak_left[0]) )
    peak_left = np.array([int(np.round(x)) if x>0 else 0 for x in peak_left])
    peak_right = np.array([int(np.round(x)) if x>0 else 0 for x in peak_right])
    
    # new start and end of chromatogram
    chrom_start = peak_left[b0]
    chrom_end = peak_right[b1-1]
    
    peak_locations = tuple(peak_locations[b0: b1] - chrom_start)
    new_seq.annotations["abif_raw"]['PLOC1'] = peak_locations
    
    for ch in sanger_channels:
        new_seq.annotations["abif_raw"][ch] = raw_data[ch][chrom_start: chrom_end+1]
        
    return new_seq


def sanger_reverse_complement(sequence):
    """
    This method returns the reverse-complement of Sanger data, 
    including the chromatorgram data and peak locations for plotting

    Parameters
    ----------
    sequence : Biopython sequence object
        The sequence to be reverse-complemented.

    Returns
    -------
    new_seq : Biopython sequence object
        The reverse-complement sequence.

    """
    
    # Built-in slicing handles sequence and quality information
    new_seq = sequence.reverse_complement()
    new_seq.annotations["abif_raw"] = {}
    
    # chromatogram data and peak locations needs to be added on
    raw_data = sequence.annotations["abif_raw"]
    
    for ch, new_ch in zip(sanger_channels, sanger_channels[::-1]):
        new_seq.annotations["abif_raw"][new_ch] = raw_data[ch][: : -1]
    
    peak_locations = np.array(raw_data['PLOC1'])
    
    peak_locations = tuple(len(raw_data[ch]) - 1 - peak_locations[: : -1])
    new_seq.annotations["abif_raw"]['PLOC1'] = peak_locations
        
    return new_seq


def plot_sanger(sequence, start_base, end_base, ax, 
                ax2=None, 
                offset=0, 
                quality_only=False, 
                letters_on_top=False,
                is_trimmed=False):
    # start_base and end_base are given in biology notation, i.e. first base of sequence is "1" (not "0")
    
    # chromatogram data
    raw_data = sequence.annotations["abif_raw"]
    channel_data = [raw_data[x] for x in sanger_channels]
    peak_locations = np.array(raw_data['PLOC1'])
    pal = sns.color_palette('dark')
    
    channel_colors = [pal[0]] + [pal[1]] + [pal[-1]] + [pal[8]]
    plot_sequence = sequence
    
    # left and right edges of each chromatogram peak
    peak_left = peak_locations
    peak_left = (peak_left[:-1] + peak_left[1:])/2
    peak_right = np.append( peak_left, peak_locations[-1] + 1/2*(peak_left[-1] - peak_left[-2]) )
    peak_left = np.insert( peak_left, 0, peak_locations[0] - 1/2*(peak_left[1] - peak_left[0]) )
    peak_left = np.array([int(np.round(x)) if x>0 else 0 for x in peak_left])
    peak_right = np.array([int(np.round(x)) if x>0 else 0 for x in peak_right])
    
    # Called sequence
    dna_seq = str(plot_sequence.seq)
    
    # Quality scores
    qual = plot_sequence.letter_annotations['phred_quality']
    
    # Plot the quality score behind everything else 
    base_positions = [i for i in range(start_base, end_base+1)]
    base_calls = dna_seq[start_base-1:end_base]
    ax.set_yticks([])
    if (end_base - start_base) < 200:
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        linewidth = None
    else:
        linewidth = 0.5
    if letters_on_top:
        ax.tick_params(labelbottom=True, labeltop=False)
        v_text = 'bottom'
        y_text = 1.02
    else:
        ax.tick_params(labelbottom=False, labeltop=True)
        v_text = 'top'
        y_text = -0.03
    qual_x = np.array( [ np.array([x - 0.5]*2) for x in base_positions ] ).flatten()
    qual_x = np.append(qual_x, [base_positions[-1]+0.5]*2)
    qual_y = np.array( [ np.array([x]*2) for x in qual[start_base-1:end_base] ] ).flatten()
    qual_y = np.append(qual_y, 0)
    qual_y = np.insert(qual_y, 0, 0)
    if is_trimmed:
        color = sns.color_palette()[3]
        f_alpha = 0.05
        l_alpha = .3
    else:
        color = sns.color_palette()[2]
        f_alpha = 0.1
        l_alpha = 1
    ax.fill_between(qual_x+offset, qual_y, color=color, alpha=f_alpha, zorder=-2)
    ax.plot(qual_x[1:-1]+offset, qual_y[1:-1], color=color, alpha=l_alpha, zorder=-1, linewidth=linewidth)
    
    ax.autoscale(enable=True)
    ylim = ax.get_ylim()
    ax.set_ylim(0, ylim[1])
    
    if not quality_only:
        ax2.set_yticks([])
        trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
        if is_trimmed:
            alpha = 0.5
        else:
            alpha = 1
        for center, base in zip(base_positions, base_calls):
            ax.text(center+offset, y_text, base, horizontalalignment='center', verticalalignment=v_text,
                   fontname="Courier New", size=20, transform=trans, alpha=alpha)
            left = peak_left[center-1]
            right = peak_right[center-1]
            for y_data, c in zip(channel_data, channel_colors):
                y = y_data[left:right+1]
                x = np.linspace(center-0.5, center+0.5, len(y)) + offset
                ax2.plot(x, y, color=c, alpha=alpha);
    
    
    