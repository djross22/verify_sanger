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
from Bio.SeqRecord import SeqRecord

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
    
def align_sanger(record1, record2, verbose=True):
    """
    This method performs a pairwise alignment of two very similar reads.
    
    The method first trims the low-quality reads from the ends of each read
    (if trim != None), and then aligns record1 with record 2
    
    Parameters
    ----------
    record1, record2 : Biopython sequence objects to be aligned
        record1 and record2 must have with phred_quality annotation for trimming
        e.g. created from raw Sanger data via 'record1 = SeqIO.read(x, "abi")'
    
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
    
    alignments = aligner.align(record1.seq, record2.seq)
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
    
    # Find mapping between indices for consensus and record1
    f_qual = []
    f_ind = []
    i = 0
    input_qual = record1.letter_annotations['phred_quality']
    start_gap = True
    for ch in align_str[0]:
        if ch == '-':
            if start_gap:
                q = 0
                ind = 'none' #-1
            else:
                if i < len(input_qual):
                    q = (input_qual[i-1] + input_qual[i])/2
                    ind = 'gap' #(i-1 + i)/2
                else:
                    q = 0
                    ind = 'none' #-1
        else:
            q = input_qual[i]
            ind = i
            i += 1
            start_gap = False
        f_qual.append(q)
        f_ind.append(ind)
    f_ind = np.array(f_ind, dtype=object)
    
    # Find mapping between indices for consensus and record2
    r_qual = []
    r_ind = []
    i = 0
    input_qual = record2.letter_annotations['phred_quality']
    start_gap = True
    for ch in align_str[2]:
        if ch == '-':
            if start_gap:
                q = 0
                ind = 'none' #-1
            else:
                if i < len(input_qual):
                    q = (input_qual[i-1] + input_qual[i])/2
                    ind = 'gap' #(i-1 + i)/2
                else:
                    q = 0
                    ind = 'none' #-1
        else:
            q = input_qual[i]
            ind = i
            i += 1
            start_gap = False
        r_qual.append(q)
        r_ind.append(ind)
    r_ind = np.array(r_ind, dtype=object)
    
    # Find consensus (using quality scores), coverage, 
    #     and the indices for any mismatches
    consensus_seq = ''
    consensus_qual = []
    mismatch_ind = []
    coverage = []
    for i, (ch1, ch2, q1, q2) in enumerate(zip(align_str[0], align_str[2], f_qual, r_qual)):
        if ch1==ch2:
            s = ch1
            q = q1 + q2
            c = 2
        else:
            if q1>=q2:
                s = ch1
                q = q1
                #if q2!=0: mismatch_ind.append(i)
            else:
                s = ch2
                q = q2
                #if q1!=0: mismatch_ind.append(i)
            c = 0 if s == '-' else 1
            mismatch_ind.append(i)
        consensus_seq += s
        consensus_qual.append(q)
        coverage.append(c)
    
    mismatch_ind = np.array(mismatch_ind)
    temp = align_str[1].strip('-')
    if len(temp)>0:
        mismatch_ind = mismatch_ind[mismatch_ind>=align_str[1].find(temp[0])]
        mismatch_ind = mismatch_ind[mismatch_ind<=align_str[1].rfind(temp[-1])]
    else:
        mismatch_ind = np.array([])
        
    # Make the consensus sequence into a Bio.SeqRecord object
    consensus_seq = SeqRecord(consensus_seq)
    consensus_seq.letter_annotations['phred_quality'] = consensus_qual
    consensus_seq.letter_annotations['coverage'] = coverage
    
    # Save results as properties of the alignment 
    align1.consensus_seq = consensus_seq
    align1.mismatch_ind = mismatch_ind
    align1.f_ind = f_ind
    align1.r_ind = r_ind
    align1.record1 = record1
    align1.record2 = record2
    align1.align_str = align_str
    
    return align1


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
        
    # This method doesn't support negative number indices, 
    #     or b0, b1>=len(sequence)
    if (b0<0) or (b1>=len(sequence)):
        raise ValueError(f'b0 = {b0} but it must be between 0 and the seqeunce length - 1 ({len(sequence)-1})')
    if (b1<0) or (b1>=len(sequence)):
        raise ValueError(f'b1 = {b1} but it must be between 0 and the seqeunce length - 1 ({len(sequence)-1})')
    
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
                letters_on_bottom=True,
                letters_on_top=False,
                is_trimmed=False,
                include_coverage=False):
    # start_base and end_base are given in biology notation, i.e. first base of sequence is "1" (not "0")
    
    if not quality_only:
        # chromatogram data
        raw_data = sequence.annotations["abif_raw"]
        channel_data = [raw_data[x] for x in sanger_channels]
        peak_locations = np.array(raw_data['PLOC1'])
        pal = sns.color_palette('dark')
        
        channel_colors = [pal[0]] + [pal[1]] + [pal[-1]] + [pal[8]]
        
        # left and right edges of each chromatogram peak
        peak_left = peak_locations
        peak_left = (peak_left[:-1] + peak_left[1:])/2
        peak_right = np.append( peak_left, peak_locations[-1] + 1/2*(peak_left[-1] - peak_left[-2]) )
        peak_left = np.insert( peak_left, 0, peak_locations[0] - 1/2*(peak_left[1] - peak_left[0]) )
        peak_left = np.array([int(np.round(x)) if x>0 else 0 for x in peak_left])
        peak_right = np.array([int(np.round(x)) if x>0 else 0 for x in peak_right])
    
    # Called sequence
    dna_seq = str(sequence.seq)
    
    # Quality scores
    qual = sequence.letter_annotations['phred_quality']
    
    # Plot the quality score behind everything else 
    base_positions = [i for i in range(start_base, end_base+1)]
    base_calls = dna_seq[start_base-1:end_base]
    ax.set_yticks([])
    
    if (end_base - start_base) < 200:
        if letters_on_bottom or letters_on_top:
            ax.xaxis.set_minor_locator(MultipleLocator(1))
        linewidth = None
    else:
        linewidth = 0.5
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
    ax.fill_between(qual_x+offset, qual_y, color=color, alpha=f_alpha, zorder=-20)
    ax.plot(qual_x[1:-1]+offset, qual_y[1:-1], color=color, alpha=l_alpha, zorder=-19, linewidth=linewidth)
    if include_coverage:
        cover = sequence.letter_annotations['coverage']
        y = np.array( [ np.array([x]*2) for x in cover[start_base-1:end_base] ] ).flatten()
        y = np.append(y, 0)
        y = np.insert(y, 0, 0)
        y = y * max(qual_y)/max(y) / 3
        color = sns.color_palette()[0]
        ax.fill_between(qual_x+offset, y, color=color, alpha=f_alpha, zorder=-10)
        ax.plot(qual_x[1:-1]+offset, y[1:-1], color=color, alpha=l_alpha, zorder=-9, linewidth=linewidth)
    
    ax.autoscale(enable=True)
    ylim = ax.get_ylim()
    ax.set_ylim(0, ylim[1])
    
    ax2.set_yticks([])
    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    if is_trimmed:
        alpha = 0.5
    else:
        alpha = 1
    ax.tick_params(labelbottom=True, labeltop=True)
    if letters_on_bottom:
        ax.tick_params(labelbottom=False)
    if letters_on_top:
        ax.tick_params(labeltop=False)
        
    for center, base in zip(base_positions, base_calls):
        if letters_on_bottom:
            ax.text(center+offset, -0.03, base, horizontalalignment='center', verticalalignment='top',
                   fontname="Courier New", size=20, transform=trans, alpha=alpha)
        if letters_on_top:
            ax.text(center+offset, 1.02, base, horizontalalignment='center', verticalalignment='bottom',
                   fontname="Courier New", size=20, transform=trans, alpha=alpha)
    
    if not quality_only:
        for center, base in zip(base_positions, base_calls):
            left = peak_left[center-1]
            right = peak_right[center-1]
            for y_data, c in zip(channel_data, channel_colors):
                y = y_data[left:right+1]
                x = np.linspace(center-0.5, center+0.5, len(y)) + offset
                ax2.plot(x, y, color=c, alpha=alpha);
    
    
    