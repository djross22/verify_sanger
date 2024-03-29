# -*- coding: utf-8 -*-
"""
Created on Mon May 10 20:42:42 2021

@author: djross
"""
import os
import glob
import copy
import numbers

import pandas as pd
import numpy as np

from Bio import Align
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord

from verify_sanger.PlottableAlignment import PlottableAlignment
from verify_sanger.PlottableRecord import PlottableRecord

import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from matplotlib.ticker import MultipleLocator
import matplotlib.patches as patches

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

# sanger_channels correspond to the ABI data traces for the letters: G, A, T, C
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
    

def num_gaps(align1, gap_seq=None):
    """
    This method calculates the number of internal gaps in a pairwise alignment

    Parameters
    ----------
    align1 : Bio.Align.PairwiseAlignment
        The alignment
        
    gap_seq : string
        string which indicates which part of the alignment to use for 
        counting gaps

    Returns
    -------
    gap_count : int
        number of internal gaps.

    """
    align_str = f'{align1}'
    align_str = align_str.split('\n')
    if gap_seq == 'seq1':
        gap_ind = 0
    elif gap_seq == 'seq2':
        gap_ind = 2
    else:
        gap_ind = 1
    gap_str = align_str[gap_ind]
    gap_str = gap_str.strip('-')
    gap_count = gap_str.count('-')
    
    return gap_count    

    
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
    
def align_sanger(record1, record2, verbose=True, find_consensus=True, 
                 default_quality1=0, include_stop=False, ungap=False):
    """
    This method performs a pairwise alignment of two very similar reads.
    
    Parameters
    ----------
    record1, record2 : Biopython sequence objects to be aligned
        record2 must have with phred_quality annotation if find_consensus is True
        e.g. created from raw Sanger data via 'record1 = SeqIO.read(x, "abi")'
    
    verbose : Boolean
        If true, method prints info about the resulting alignments
    
    find_consensus : Boolean
        If true, method uses quality scores to detirmine a consensus sequence
    
    default_quality1 : float
        The quality score assumed for record1 if it does not have a 
        'phred_quality' entry in the letter_annotations dictionary, 
        e.g., if record1 iS a reference sequence
    
    ungap : Boolean
        If true, method removes gap characters from input sequences before 
        running alignment

    Returns
    -------
    align1 : the resulting alignment, as a Bio.Align.PairwiseAlignment object
        the zeroth element from the iterator from the pairwise alignment call
        
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
    if include_stop:
        aligner.alphabet += '*'
    
    if ungap:
        record1 = ungap_seqrecord(record1)
        record2 = ungap_seqrecord(record2)
    alignments = aligner.align(record1.seq, record2.seq)
    if verbose: 
        try:
            print(f'{len(alignments)} alignment(s) found with score: {alignments.score}')
        except OverflowError:
            print(f'Overflow number of alignment(s) found with score: {alignments.score}')
    
    align1 = alignments[0]
    if verbose: print(f'{align1.aligned}')
    align_str = f'{align1}'
    align_str = align_str.split('\n')
    if len(align1.aligned[0])>1:
        
        target_str = ''
        for s in align_str:
            if 'target' in s:
                s = s.replace('target', '')
                s = s.strip()
                s = s[s.find(' ')+1:]
                if ' ' in s:
                    s = s[:s.find(' ')]
                target_str += s
        target_str = target_str.strip('-')
        
        if '--' in target_str:
            if verbose: print('multi-base gap in sequence no. 1')
            #x = target_str.find('--')
            #print(f'{x}, {target_str}')
        elif '-' in target_str:
            count = target_str.count('-')
            if verbose: print(f'{count} gap(s) in sequence no. 1')
            #x = target_str.find('-')
            #print(f'{x}, {target_str}')
            
        query_str = align_str[2]
        
        query_str = ''
        for s in align_str:
            if 'query' in s:
                s = s.replace('query', '')
                s = s.strip()
                s = s[s.find(' ')+1:]
                if ' ' in s:
                    s = s[:s.find(' ')]
                query_str += s
        
        query_str = query_str.strip('-')
        if '--' in query_str:
            if verbose: print('multi-base gap in sequence no. 2')
            #print(f'{query_str}')
        elif '-' in query_str:
            count = query_str.count('-')
            if verbose: print(f'{count} gap(s) in sequence no. 2')
            #print(f'{query_str}')
    else:
        if verbose: print('no gaps')
        
    match_str = align_str[1]
    count = match_str.count('|')
    if verbose: print(f'{count} matches in alignment')
    count = match_str.count('.')
    if verbose: print(f'{count} mismatches in alignment')
    
    # Remake list of quality scores for record1, with gaps based on alignment
    f_qual = []
    i = 0
    if find_consensus:
        if 'phred_quality' in record1.letter_annotations.keys():
            input_qual = record1.letter_annotations['phred_quality']
        else:
            input_qual = [default_quality1]*len(record1)
    else:
        input_qual = [0]*len(record1)
    start_gap = True
    for ch in align_str[0]:
        if ch == '-':
            if start_gap:
                q = 0
            else:
                if i < len(input_qual):
                    q = (input_qual[i-1] + input_qual[i])/2
                else:
                    q = 0
        else:
            q = input_qual[i]
            i += 1
            start_gap = False
        f_qual.append(q)
    
    # Remake list of quality scores for record2, with gaps based on alignment
    r_qual = []
    i = 0
    if find_consensus:
        input_qual = record2.letter_annotations['phred_quality']
    else:
        input_qual = [0]*len(record2)
    start_gap = True
    for ch in align_str[2]:
        if ch == '-':
            if start_gap:
                q = 0
            else:
                if i < len(input_qual):
                    q = (input_qual[i-1] + input_qual[i])/2
                else:
                    q = 0
        else:
            q = input_qual[i]
            i += 1
            start_gap = False
        r_qual.append(q)
    
    # Find consensus (using quality scores), coverage, 
    #     and the indices for any mismatches
    # The mismatch_ind list is a list of mismatches and gaps between the 
    #     two input sequences. It is useful to have even when the consensus
    #     sequence is not.
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
        
    if find_consensus:
        # Make the consensus sequence into a PlottableRecord object
        consensus_seq = SeqRecord(Seq(consensus_seq))
        consensus_seq.letter_annotations['phred_quality'] = consensus_qual
        consensus_seq = PlottableRecord(consensus_seq)
        consensus_seq.coverage = coverage
    else:
        consensus_seq = None
    
    # Save results as a PlottableAlignment
    align1 = PlottableAlignment(input_alignment=align1, 
                                consensus_seq=consensus_seq,
                                mismatch_ind=mismatch_ind,
                                record1=record1,
                                record2=record2)
    
    return align1


def zoom_out_plot(align1, title=None, 
                  seq1_label='Sequence 1:', seq2_label='Sequence 2:'):
        
    plt.rcParams["figure.figsize"] = [12, 4]
    fig, axs = plt.subplots(3, 1)
    if title is not None:
        fig.suptitle(title, size=20, x=0.1,
                     verticalalignment='bottom',
                     horizontalalignment='left')
    axs[0].get_shared_x_axes().join(*axs)
    #axs[0].get_shared_y_axes().join(*axs)
    ax2 = [ ax.twinx()  for ax in axs ]
    
    f_seq = align1.record1
    r_seq = align1.record2
    consensus_seq = align1.consensus_seq
    
    plot_sanger(consensus_seq, axs[0], ax2=ax2[0], 
                include_chromatograms=False, include_coverage=True, letters_on_bottom=False)
    
    plot_sanger(f_seq, axs[1], ax2=ax2[1], 
                include_chromatograms=False, letters_on_bottom=False)
    
    plot_sanger(r_seq, axs[2], ax2=ax2[2], 
                include_chromatograms=False, letters_on_bottom=False)
    
    axs[0].tick_params(labelbottom=False)
    axs[1].tick_params(labelbottom=False)
    axs[1].tick_params(labeltop=False)
    axs[2].tick_params(labeltop=False)
    
    for ax, s_label in zip(axs, ['Consensus:', seq1_label, seq2_label]):
        ax.text(-0.01, 0.5, s_label, horizontalalignment='right', verticalalignment='center',
                size=20, transform=ax.transAxes)
    
    shift = -0.02
    for i, ax in enumerate(axs):
        box = ax.get_position()
        box.y0 = box.y0 - shift*i
        box.y1 = box.y1 - shift*i
        ax.set_position(box)
        
    return fig, axs
        

def zoom_in_plot(align1, index, zoom_span=10, title=None, verbose=False,
                 seq1_label='Sequence 1:', seq2_label='Sequence 2:',
                 include_chromatograms=True, compare_to_ref=False,
                 anchor_feature=None):
    
    if isinstance(index, numbers.Integral):
        index = slice(index, index+1)
    zoom = slice(index.start-zoom_span, index.stop+zoom_span+1)
    
    if verbose:
        print(align1.align_str[0][zoom])
        print(align1.align_str[1][zoom])
        print(align1.align_str[2][zoom])
    
    plt.rcParams["figure.figsize"] = [12, 5]
    fig, axs = plt.subplots(3, 1)
    if title is not None:
        fig.suptitle(title, size=20, x=0.1,
                     verticalalignment='bottom',
                     horizontalalignment='left')
    axs[0].get_shared_x_axes().join(*axs)
    axs[0].get_shared_y_axes().join(*axs)
    ax2 = [ ax.twinx()  for ax in axs ]
    
    f_seq = align1.record1
    r_seq = align1.record2
    consensus_seq = align1.consensus_seq
    
    if anchor_feature is not None:
        anchor_offset = 0
        
        anchor_feat = None
        for feat in align1.record1.features:
            if feat.qualifiers['label'][0] == anchor_feature:
                anchor_feat = feat
                break
        if anchor_feat is not None:
            anchor_offset = -anchor_feat.location.start.position
        else:
            anchor_offset = 0
    else:
        anchor_offset = 0
    
    plot_sanger(consensus_seq, axs[0], zoom.start+1, zoom.stop, 
                ax2=ax2[0], offset=anchor_offset,
                include_chromatograms=False, include_coverage=True)
    
    plot_sanger(f_seq, axs[1], zoom.start+1, zoom.stop, ax2=ax2[1], offset=anchor_offset,
                letters_on_top=True, include_chromatograms=include_chromatograms,
                ref_seq_plot=compare_to_ref)
    
    plot_sanger(r_seq, axs[2], zoom.start+1, zoom.stop, ax2=ax2[2], offset=anchor_offset,
                letters_on_top=True, include_chromatograms=include_chromatograms)
    
    for ax, s_label in zip(axs, ['Consensus:', seq1_label, seq2_label]):
        ax.text(-0.01, 0.5, s_label, horizontalalignment='right', 
                verticalalignment='center',
                size=20, transform=ax.transAxes)
        
    shift = 0.06
    for i, ax in enumerate(axs):
        box = ax.get_position()
        box.y0 = box.y0 - shift*i
        box.y1 = box.y1 - shift*i
        ax.set_position(box)
        ylim = ax.get_ylim()
        ax.set_ylim(ylim)
        for x in range(index.start, index.stop):
            rect = patches.Rectangle((x+0.5+anchor_offset, ylim[0]), 1, ylim[1]-ylim[0], 
                                     linewidth=1, edgecolor='k', facecolor='r', alpha=0.2, zorder=-100)
            ax.add_patch(rect)
        
    if compare_to_ref:
        axs[1].set_axis_off()
        ax2[1].set_axis_off()
        shift1 = 0.135
        shift2 = shift1 - 0.08
        box = axs[1].get_position()
        box.y0 = box.y0 + shift1 + shift2
        box.y1 = box.y1 + shift2
        axs[1].set_position(box)
        
        shift3 = 0.05
        box = axs[2].get_position()
        box.y0 = box.y0 + shift1 + shift2 + shift3
        box.y1 = box.y1 + shift1 + shift2 + shift3
        axs[2].set_position(box)
    
    return fig, axs


def make_type_color_dict(record1):
    type_list = np.unique([x.type for x in record1.features])
    return {t:c for t, c in zip(type_list, sns.color_palette('pastel')) }


def plot_features(record1, ax, x_offset=0):
    for feat in record1.features:
        p0 = feat.location.start.position + x_offset
        p1 = feat.location.end.position + x_offset
        c = make_type_color_dict(record1)[feat.type]
        rect = patches.Rectangle((p0+0.5, -0.5), p1-p0, 1, 
                                 linewidth=1, edgecolor='k', facecolor=c, alpha=1, zorder=20)
        ax.add_patch(rect)
        if (p1-p0)>=len(record1)/5:
            label = feat.qualifiers['label'][0]
            ax.text((p1+p0)/2, -0.1, label, horizontalalignment='center', verticalalignment='center',
                     size=20, transform=ax.transData, alpha=1, zorder=30)
        ax.plot([0.5 + x_offset, len(record1)+0.5 + x_offset], [0,0], 'k', zorder=10);
        ax.set_axis_off()


def compare_to_ref_plot(align1, title=None, 
                        seq1_label='Reference:', seq2_label='Sanger:',
                        anchor_feature=None):
    
    plt.rcParams["figure.figsize"] = [12, 3]
    fig, axs = plt.subplots(2, 1)
    if title is not None:
        fig.suptitle(title, size=20, x=0.1,
                     verticalalignment='bottom',
                     horizontalalignment='left')
    axs[0].get_shared_x_axes().join(*axs)
    #axs[0].get_shared_y_axes().join(*axs)
    
    #f_seq = align1.record1
    r_seq = align1.record2
    #consensus_seq = align1.consensus_seq
        
    #for b, x in zip(f_block, f_offset):
    #    vs.plot_sanger(f_seq, b[0]+1, b[1]+1, axs[0], ax2=None, offset=x, 
    #                include_chromatograms=False, letters_on_bottom=False, include_quality=False)
    if anchor_feature is not None:
        anchor_feat = None
        for feat in align1.record1.features:
            if feat.qualifiers['label'][0] == anchor_feature:
                anchor_feat = feat
                break
        if anchor_feat is not None:
            anchor_offset = -anchor_feat.location.start.position
        else:
            anchor_offset = 0
    else:
        anchor_offset = 0
        
    plot_features(align1.record1, axs[0], x_offset=anchor_offset)
    
    plot_sanger(r_seq, axs[1], ax2=None, offset=anchor_offset, 
                include_chromatograms=False, letters_on_bottom=False)
            
    # Add red rectangles above Sanger quality plot to show locations of 
    #     mismatches and gaps
    ylim = axs[1].get_ylim()
    rect_h = (ylim[1] - ylim[0])*0.15
    rect_y0 = ylim[1]
    for ind in align1.mismatch_ind:
        rect_x0 = ind + 0.5 + anchor_offset
        rect = patches.Rectangle((rect_x0, rect_y0), 1, rect_h, 
                                 linewidth=2, edgecolor='r', facecolor='r', 
                                 alpha=1, zorder=100)
        axs[1].add_patch(rect)
    axs[1].set_ylim(ylim[0], ylim[1]+rect_h*2)
    
    axs[0].tick_params(labelbottom=False)
    axs[1].tick_params(labelbottom=True)
    axs[1].tick_params(labeltop=False)
    
    for ax, s_label in zip(axs, [seq1_label, seq2_label]):
        ax.text(-0.01, 0.5, s_label, horizontalalignment='right', verticalalignment='center',
                size=20, transform=ax.transAxes)
    
    shift = 0.2
    box = axs[0].get_position()
    box.y0 = box.y0 + shift
    axs[0].set_position(box)
    box = axs[1].get_position()
    box.y0 = box.y0 + shift
    box.y1 = box.y1 + shift
    axs[1].set_position(box)
        
    return fig, axs


def is_good_sanger(alignment, min_matches=20, max_mismatch_ratio=0.02,
                   max_gap_ratio=0.015):
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
        The maximum value of the ratio: (number of mismatches)/(number of matches)
        for a good alignment. The default is 0.02.
    max_gap_ratio : float, optional
        The maximum value of the ratio: (number of gaps)/(number of matches)
        for a good alignment. Calculated separately for each sequence in the
        alignment. The default is 0.015.

    Returns
    -------
    ret_val : Boolean
        True if the alignement is good.

    """
    
    x, y = num_matches(alignment)
    ret_val = x >= min_matches
    if ret_val:
        ret_val = ret_val & (y/x <= max_mismatch_ratio)
        z1 = num_gaps(alignment, gap_seq='seq1')
        z2 = num_gaps(alignment, gap_seq='seq2')
        ret_val = ret_val & (z1/x <= max_gap_ratio)
        ret_val = ret_val & (z2/x <= max_gap_ratio)
    
    return ret_val


def sanger_reverse_complement(record1):
    """
    This method returns the reverse-complement of Sanger data, 
    including the chromatorgram data and peak locations for plotting

    Parameters
    ----------
    record1 : Biopython sequence record object
        The sequence to be reverse-complemented.

    Returns
    -------
    new_record : Biopython sequence record object
        The reverse-complement sequence.

    """
    
    # Built-in slicing handles sequence and quality information
    new_record = record1.reverse_complement()
    new_record.annotations["abif_raw"] = {}
    
    # chromatogram data and peak locations needs to be added on
    raw_data = record1.annotations["abif_raw"]
    #TODO: transfer other properties/attributes
    
    for ch, new_ch in zip(sanger_channels, sanger_channels[::-1]):
        new_record.annotations["abif_raw"][new_ch] = raw_data[ch][: : -1]
    
    peak_locations = np.array(raw_data['PLOC1'])
    
    peak_locations = tuple(len(raw_data[ch]) - 1 - peak_locations[: : -1])
    new_record.annotations["abif_raw"]['PLOC1'] = peak_locations
        
    return new_record


def ungap_seqrecord(record1, inplace=False):
    base_array = np.array([x for x in f'{record1.seq}'])
    remove_list = np.where(base_array=='-')[0]
    
    if inplace:
        new_record = record1
    else:
        new_record = copy.deepcopy(record1)
        
    new_anotations = {}
    for key in new_record.letter_annotations.keys():
        new_val = np.delete(new_record.letter_annotations[key], remove_list)
        new_anotations[key] = new_val
    
    new_record.letter_annotations = {}
    
    new_record.seq = new_record.seq.ungap('-')
    new_record.letter_annotations = new_anotations
    
    if inplace:
        return
    else:
        return new_record


def plot_sanger(record, ax, 
                start_base=0, 
                end_base=None, 
                ax2=None, 
                offset=0, 
                include_chromatograms=True, 
                letters_on_bottom=True,
                letters_on_top=False,
                is_trimmed=False,
                include_coverage=False,
                include_quality=True,
                ref_seq_plot=False):
    
    if end_base is None:
        end_base = len(record)
    
    # Don't show y-axis ticks
    ax.set_yticks([])
    if ax2 is not None:
        ax2.set_yticks([])
        ax2.get_shared_x_axes().join(ax, ax2)
    
    if ref_seq_plot:
        include_quality = False
        include_coverage = False
        include_chromatograms = False
        letters_on_bottom = False
        letters_on_top = False
        letters_in_middle = True
    else:
        letters_in_middle = False
        
    # start_base and end_base are given in biology notation, i.e. first base of sequence is "1" (not "0")
    # Also, end_base is the last plotted base, not the last+1
    if start_base<1:
        start_base = 1
    if end_base>len(record):
        end_base = len(record)
    
    # Plot the chromatogram data
    if include_chromatograms :
        chrom_data = record.chromatogram_plot_data(start_base-1, end_base)
        
        pal = sns.color_palette('dark')
        channel_colors = [pal[0]] + [pal[1]] + [pal[-1]] + [pal[8]]
        
        for y, c in zip(chrom_data[1:], channel_colors):
            ax2.plot(chrom_data[0], y, color=c);
    
    if (end_base - start_base) < 200:
        if letters_on_bottom or letters_on_top:
            ax.xaxis.set_minor_locator(MultipleLocator(1))
        linewidth = None
    else:
        linewidth = 0.5
        
    # Plot the quality score (behind everything else) 
    if include_quality:
        qual_data = record.quality_plot_data(start_base-1, end_base)
        qual_x = qual_data[0]
        qual_y = qual_data[1]
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
    
    # Plot the sequence coverage (behind everything else)
    if include_coverage:
        cover_data = record.coverage_plot_data(start_base-1, end_base)
        cover_x = cover_data[0]
        cover_y = cover_data[1]
        
        cover_y = cover_y * max(qual_y)/max(cover_y) / 3
        color = sns.color_palette()[0]
        ax.fill_between(cover_x+offset, cover_y, color=color, alpha=f_alpha, zorder=-10)
        ax.plot(cover_x[1:-1]+offset, cover_y[1:-1], color=color, alpha=l_alpha, zorder=-9, linewidth=linewidth)
    
    ax.autoscale(enable=True)
    ylim = ax.get_ylim()
    ax.set_ylim(0, ylim[1])
    
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
    if letters_in_middle:
        ax.tick_params(labeltop=False)
        ax.tick_params(labelbottom=False)
        
    # Label x-axis with called sequence
    dna_seq = str(record.seq)
    base_calls = dna_seq[start_base-1:end_base]
    base_positions = [i+1 for i in range(start_base-1, end_base)]
    for center, base in zip(base_positions, base_calls):
        if letters_on_bottom:
            ax.text(center+offset, -0.03, base, horizontalalignment='center', verticalalignment='top',
                   fontname="Courier New", size=20, transform=trans, alpha=alpha)
        if letters_on_top:
            ax.text(center+offset, 1.02, base, horizontalalignment='center', verticalalignment='bottom',
                   fontname="Courier New", size=20, transform=trans, alpha=alpha)
        if letters_in_middle:
            ax.text(center+offset, 0.5, base, horizontalalignment='center', verticalalignment='center',
                   fontname="Courier New", size=20, transform=trans, alpha=alpha)
    
    
def translate_with_gaps(record1):
    seq1 = record1.seq.ungap(gap='-')
    new_record = SeqRecord(seq1)
    return new_record.translate()


def get_sequence_vs_reference(ref_alignment, ref_feature):
    
    # Find the feature in the reference SeqRecord from the ref_alignment
    reference_feat = None
    for feat in ref_alignment.record1.features:
        if feat.qualifiers['label'][0] == ref_feature:
            reference_feat = feat
            break
        
    if reference_feat is None:
        return None, None
    
    # These are the start and end positions of the insert CDS in the aligned reference sequence 
    ref_start = reference_feat.location.start.position
    ref_end = reference_feat.location.end.position
    
    test_dna = ref_alignment.consensus_seq[ref_start: ref_end]
    
    test_aminos = translate_with_gaps(test_dna)
    
    return test_dna, test_aminos


def find_mutations_vs_reference(ref_alignment, ref_feature, verbose=True):
    
    # Find the feature in the reference SeqRecord from the ref_alignment
    reference_feat = None
    for feat in ref_alignment.record1.features:
        if feat.qualifiers['label'][0] == ref_feature:
            reference_feat = feat
            break
        
    if reference_feat is None:
        return None, None
    
    # These are the start and end positions of the insert CDS in the aligned reference sequence 
    ref_start = reference_feat.location.start.position
    ref_end = reference_feat.location.end.position
    
    ref_dna = ref_alignment.record1[ref_start: ref_end]
    test_dna = ref_alignment.consensus_seq[ref_start: ref_end]
    
    ref_aminos = translate_with_gaps(ref_dna)
    test_aminos = translate_with_gaps(test_dna)
    
    # The same aligner settings used for Sanger sequence alignments
    #     also seem to work for amino acid alignments.
    amino_alignment = align_aminos(ref_aminos, test_aminos, verbose=verbose)
    
    substitution_codes, indel_codes = mutation_codes_from_alignment(amino_alignment)
                    
    return amino_alignment, substitution_codes, indel_codes


def align_aminos(ref_aminos, test_aminos, verbose=True):
    # The same aligner settings used for Sanger sequence alignments
    #     also seem to work for amino acid alignments.
    amino_alignment = align_sanger(ref_aminos, test_aminos, 
                                   find_consensus=False, include_stop=True,
                                   verbose=verbose)
    
    return amino_alignment


def mutation_codes_from_alignment(align1):
    substitution_codes = []
    indel_codes = []
    for ind in align1.mismatch_ind:
        wt_amino = align1.align_str[0][ind]
        new_amino = align1.align_str[2][ind]
        if align1.align_str[1][ind] == '.':
            # Substitution
            substitution_codes.append(f'{wt_amino}{ind+1}{new_amino}')
        elif align1.align_str[1][ind] == '-':
            # Insertion or deletion
            #TODO: remove duplicate entries in indel_codes for 
            #     multi-base/multi-residue indels
            if wt_amino == '-':
                # Insertion
                p1 = ind
                ch = align1.align_str[0][p1]
                while ch =='-':
                    p1 -= 1
                    ch = align1.align_str[0][p1]
                p2 = ind
                ch = align1.align_str[0][p2]
                while ch =='-':
                    p2 += 1
                    ch = align1.align_str[0][p2]
    
                wt_amino1 = f'{align1.record1.seq}'[p1]
                wt_amino2 = f'{align1.record1.seq}'[p2]
            
                indel_codes.append(f'{wt_amino1}{p1+1}_{wt_amino2}{p2+1}ins{new_amino}')
            elif new_amino == '-':
                # Deletion
                p1 = ind
                ch = align1.align_str[2][p1]
                while ch =='-':
                    p1 -= 1
                    ch = align1.align_str[2][p1]
                p1 += 1
                p2 = ind
                ch = align1.align_str[2][p2]
                while ch =='-':
                    p2 += 1
                    ch = align1.align_str[2][p2]
                p2 -= 1
    
                wt_amino1 = f'{align1.record1.seq}'[p1]
                wt_amino2 = f'{align1.record1.seq}'[p2]
                
                if p1==p2:
                    indel_codes.append(f'{wt_amino1}{p1+1}del')
                else:
                    indel_codes.append(f'{wt_amino1}{p1+1}_{wt_amino2}{p2+1}del')
                    
    return substitution_codes, indel_codes


def mutation_codes_from_name(name, added_codes=[]):
    name = name[name.find('(')+1:]
    name = name[:name.rfind(')')]
    codes = name.split('/') + added_codes
    codes = [x[:x.find('-')] if '-' in x else x for x in codes]
    pos = [int(x[1:-1]) for x in codes]
    df = pd.DataFrame({'codes':codes, 'pos':pos})
    df.sort_values(by='pos', inplace=True)
    codes = list(df.codes)
    return codes

