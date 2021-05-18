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
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord

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
    
    The method first trims the low-quality reads from the ends of each read
    (if trim != None), and then aligns record1 with record 2
    
    Parameters
    ----------
    record1, record2 : Biopython sequence objects to be aligned
        record1 and record2 must have with phred_quality annotation for trimming
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
    if include_stop:
        aligner.alphabet += '*'
    
    if ungap:
        record1 = ungap_seqrecord(record1)
        record2 = ungap_seqrecord(record2)
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
            if verbose: print('multi-base gap in sequence no. 1')
            #x = target_str.find('--')
            #print(f'{x}, {target_str}')
        elif '-' in target_str:
            count = target_str.count('-')
            if verbose: print(f'{count} gap(s) in sequence no. 1')
            #x = target_str.find('-')
            #print(f'{x}, {target_str}')
            
        query_str = align_str[2]
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
    
    # Find mapping between indices for consensus and record1
    f_qual = []
    f_ind = []
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
    if find_consensus:
        input_qual = record2.letter_annotations['phred_quality']
    else:
        input_qual = [0]*len(record2)
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
        # Make the consensus sequence into a Bio.SeqRecord object
        consensus_seq = SeqRecord(Seq(consensus_seq))
        consensus_seq.letter_annotations['phred_quality'] = consensus_qual
        consensus_seq.letter_annotations['coverage'] = coverage
    
    # Save results as properties of the alignment 
    if find_consensus:
        align1.consensus_seq = consensus_seq
    align1.mismatch_ind = mismatch_ind
    align1.f_ind = f_ind
    align1.r_ind = r_ind
    align1.record1 = record1
    align1.record2 = record2
    align1.align_str = align_str
    
    return align1


def make_blocks(block, breaks):
    if len(breaks)>0:
        temp = []
        next_x = 0
        for x in breaks:
            y = block[next_x: x]
            if len(y)>0:
                temp.append([y[0], y[-1]])
            else:
                temp.append([block[0], block[-1]])
            next_x = x + 1
        y = block[next_x: ]
        if len(y)>0:
            temp.append([y[0], y[-1]])
        else:
            temp.append([block[0], block[-1]])
    else:
        if len(block)>0:
            temp = [ [block[0], block[-1]] ]
        else:
            temp = None
    return temp


def zoom_out_plot(align1, title=None, 
                  seq1_label='Sequence 1:', seq2_label='Sequence 2:'):
    f_block = align1.f_ind
    f_block = f_block[f_block!='none']
    f_breaks = np.where(f_block=='gap')[0]
    r_block = align1.r_ind
    r_block = r_block[r_block!='none']
    r_breaks = np.where(r_block=='gap')[0]
    f_block = make_blocks(f_block, f_breaks)
    r_block = make_blocks(r_block, r_breaks)
    f_offset = [ np.where(align1.f_ind==x[0])[0][0]-x[0] for x in f_block ]
    r_offset = [ np.where(align1.r_ind==x[0])[0][0]-x[0] for x in r_block ]
    
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
    
    plot_sanger(consensus_seq, 1, len(consensus_seq), axs[0], ax2=ax2[0], offset=0,
                include_chromatograms=False, include_coverage=True, letters_on_bottom=False)
    
    for b, x in zip(f_block, f_offset):
        plot_sanger(f_seq, b[0]+1, b[1]+1, axs[1], ax2=ax2[1], offset=x, 
                    include_chromatograms=False, letters_on_bottom=False)
    
    for b, x in zip(r_block, r_offset):
        plot_sanger(r_seq, b[0]+1, b[1]+1, axs[2], ax2=ax2[2], offset=x, 
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
        

def zoom_in_plot(align1, zoom_ind, zoom_span=10, title=None, verbose=False,
                 seq1_label='Sequence 1:', seq2_label='Sequence 2:',
                 include_chromatograms=True, compare_to_ref=False):
    
    f_block = align1.f_ind[zoom_ind-zoom_span: zoom_ind+zoom_span+1]
    f_block = f_block[f_block!='none']
    f_breaks = np.where(f_block=='gap')[0]
    r_block = align1.r_ind[zoom_ind-zoom_span: zoom_ind+zoom_span+1]
    r_block = r_block[r_block!='none']
    r_breaks = np.where(r_block=='gap')[0]
    f_block = make_blocks(f_block, f_breaks)
    r_block = make_blocks(r_block, r_breaks)
    if (r_block is None) or (f_block is None):
        if verbose: print('No good sequence blocks to plot')
        return None, None
    f_block = [x for x in f_block if 'gap' not in x]
    r_block = [x for x in r_block if 'gap' not in x]
    f_offset = [ np.where(align1.f_ind==x[0])[0][0]-x[0] for x in f_block ]
    r_offset = [ np.where(align1.r_ind==x[0])[0][0]-x[0] for x in r_block ]
    if verbose:
        print(align1.align_str[0][zoom_ind-zoom_span:zoom_ind+zoom_span+1])
        print(align1.align_str[1][zoom_ind-zoom_span:zoom_ind+zoom_span+1])
        print(align1.align_str[2][zoom_ind-zoom_span:zoom_ind+zoom_span+1])
    
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
    
    plot_sanger(consensus_seq, zoom_ind-zoom_span+1, zoom_ind+zoom_span+1, axs[0], ax2=ax2[0], offset=0,
                include_chromatograms=False, include_coverage=True)
    
    for b, x in zip(f_block, f_offset):
        plot_sanger(f_seq, b[0]+1, b[1]+1, axs[1], ax2=ax2[1], offset=x,
                    letters_on_top=True, include_chromatograms=include_chromatograms,
                    ref_seq_plot=compare_to_ref)
    
    for b, x in zip(r_block, r_offset):
        plot_sanger(r_seq, b[0]+1, b[1]+1, axs[2], ax2=ax2[2], offset=x,
                    letters_on_top=True, include_chromatograms=include_chromatograms)
    
    for ax, s_label in zip(axs, ['Consensus:', seq1_label, seq2_label]):
        ax.text(-0.01, 0.5, s_label, horizontalalignment='right', verticalalignment='center',
                size=20, transform=ax.transAxes)
        
    shift = 0.06
    for i, ax in enumerate(axs):
        box = ax.get_position()
        box.y0 = box.y0 - shift*i
        box.y1 = box.y1 - shift*i
        ax.set_position(box)
        ylim = ax.get_ylim()
        ax.set_ylim(ylim)
        rect = patches.Rectangle((zoom_ind+0.5, ylim[0]), 1, ylim[1]-ylim[0], 
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
    f_block = align1.f_ind
    f_block = f_block[f_block!='none']
    f_breaks = np.where(f_block=='gap')[0]
    r_block = align1.r_ind
    r_block = r_block[r_block!='none']
    r_breaks = np.where(r_block=='gap')[0]
    f_block = make_blocks(f_block, f_breaks)
    r_block = make_blocks(r_block, r_breaks)
    f_offset = [ np.where(align1.f_ind==x[0])[0][0]-x[0] for x in f_block ]
    r_offset = [ np.where(align1.r_ind==x[0])[0][0]-x[0] for x in r_block ]
    
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
            anchor_offset = -np.mean(f_offset) - anchor_feat.location.start.position
        else:
            anchor_offset = 0
    else:
        anchor_offset = 0
        
    plot_features(align1.record1, axs[0], x_offset=np.mean(f_offset)+anchor_offset)
    
    for b, x in zip(r_block, r_offset):
        if (type(b[0]) is int) and (type(b[1]) is int):
            plot_sanger(r_seq, b[0]+1, b[1]+1, axs[1], ax2=None, offset=x+anchor_offset, 
                        include_chromatograms=False, letters_on_bottom=False)
            
    # Add red rectangles above Sanger quality plot to show locations of 
    #     mismatches and gaps
    ylim = axs[1].get_ylim()
    rect_h = (ylim[1] - ylim[0])*0.15
    rect_y0 = ylim[1]
    for ind in align1.mismatch_ind:
        rect_x0 = ind + 0.5 + np.mean(f_offset) + anchor_offset
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


def slice_sanger(record1, b0=0, b1=None):
    """
    This method slices Sanger data, including the chromatorgram data
    and peak locations for plotting

    Parameters
    ----------
    record1 : Biopython sequence record object
        The sequence to be sliced.
    b0 : int, optional
        The first index for the slice. The default is 0.
    b1 : int, optional
        The end index for the slice. The default is None.

    Returns
    -------
    new_record : Biopython sequence record object
        The sliced sequence.

    """
    # If b0 and/or b1 are None, still trim chromatorgram data
    if b1 == None:
        b1 = len(record1)
        
    # This method doesn't support negative number indices, 
    #     or b0, b1>=len(record1)
    if (b0<0) or (b1>=len(record1)):
        raise ValueError(f'b0 = {b0} but it must be between 0 and the seqeunce length - 1 ({len(record1)-1})')
    if (b1<0) or (b1>=len(record1)):
        raise ValueError(f'b1 = {b1} but it must be between 0 and the seqeunce length - 1 ({len(record1)-1})')
    
    # Built-in slicing handles sequenc and quality information
    new_record = record1[b0: b1]
    new_record.annotations["abif_raw"] = {}
    
    # chromatogram data and peak locations needs to be added on
    raw_data = record1.annotations["abif_raw"]
    peak_locations = np.array(raw_data['PLOC1'])
    #TODO: transfer other properties/attributes
    
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
    new_record.annotations["abif_raw"]['PLOC1'] = peak_locations
    
    for ch in sanger_channels:
        new_record.annotations["abif_raw"][ch] = raw_data[ch][chrom_start: chrom_end+1]
        
    return new_record


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


def ungap_seqrecord(record1):
    base_array = np.array([x for x in f'{record1.seq}'])
    remove_list = np.where(base_array=='-')[0]
    
    new_anotations = {}
    for key in record1.letter_annotations.keys():
        new_val = np.delete(record1.letter_annotations[key], remove_list)
        new_anotations[key] = new_val
    
    record1.letter_annotations = {}
    
    record1.seq = record1.seq.ungap('-')
    record1.letter_annotations = new_anotations
    
    return record1


def plot_sanger(sequence, start_base, end_base, ax, 
                ax2=None, 
                offset=0, 
                include_chromatograms=True, 
                letters_on_bottom=True,
                letters_on_top=False,
                is_trimmed=False,
                include_coverage=False,
                include_quality=True,
                ref_seq_plot=False):
    
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
    if start_base<1:
        start_base = 1
    if end_base>len(sequence):
        end_base = len(sequence)
    
    if include_chromatograms :
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
    
    # Plot the quality score and coverage behind everything else 
    base_positions = [i for i in range(start_base, end_base+1)]
    base_calls = dna_seq[start_base-1:end_base]
    ax.set_yticks([])
    qual_x = np.array( [ np.array([x - 0.5]*2) for x in base_positions ] ).flatten()
    qual_x = np.append(qual_x, [base_positions[-1]+0.5]*2)
    
    if (end_base - start_base) < 200:
        if letters_on_bottom or letters_on_top:
            ax.xaxis.set_minor_locator(MultipleLocator(1))
        linewidth = None
    else:
        linewidth = 0.5
        
    # Quality scores
    if include_quality:
        qual = sequence.letter_annotations['phred_quality']
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
    
    if ax2 is not None:
        ax2.set_yticks([])
        
    if include_chromatograms:
        for center, base in zip(base_positions, base_calls):
            left = peak_left[center-1]
            right = peak_right[center-1]
            for y_data, c in zip(channel_data, channel_colors):
                y = y_data[left:right+1]
                x = np.linspace(center-0.5, center+0.5, len(y)) + offset
                ax2.plot(x, y, color=c, alpha=alpha);
    
    
def translate_with_gaps(record1):
    seq1 = record1.seq.ungap(gap='-')
    new_record = SeqRecord(seq1)
    return new_record.translate()


def find_mutations_vs_reference(ref_alignment, ref_feature, verbose=True):
    
    # Find the feature in the reference SeqRecord from the ref_alignment
    reference_feat = None
    for feat in ref_alignment.record1.features:
        if feat.qualifiers['label'][0] == ref_feature:
            reference_feat = feat
            break
        
    if reference_feat is None:
        return None, None
    
    # These are the start and end positions of the insert CDS in the reference sequence 
    ref_b0 = reference_feat.location.start.position
    ref_b1 = reference_feat.location.end.position
    
    # These are the start and end positions of the insert CDS in the alignment/consensus sequence
    ref_start = np.where(ref_alignment.f_ind==ref_b0)[0][0]
    ref_end = np.where(ref_alignment.f_ind==ref_b1)[0][0]
    
    ref_dna = ref_alignment.record1[ref_b0: ref_b1]
    test_dna = ref_alignment.consensus_seq[ref_start: ref_end]
    
    ref_aminos = ref_dna.translate()
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
        wt_pos = align1.f_ind[ind]
        wt_amino = align1.align_str[0][ind]
        new_amino = align1.align_str[2][ind]
        if align1.align_str[1][ind] == '.':
            # Substitution
            substitution_codes.append(f'{wt_amino}{wt_pos+1}{new_amino}')
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
    
                wt_pos1 = align1.f_ind[p1]
                wt_pos2 = align1.f_ind[p2]
                wt_amino1 = f'{align1.record1.seq}'[wt_pos1]
                wt_amino2 = f'{align1.record1.seq}'[wt_pos2]
            
                indel_codes.append(f'{wt_amino1}{wt_pos1+1}_{wt_amino2}{wt_pos2+1}ins{new_amino}')
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
    
                wt_pos1 = align1.f_ind[p1]
                wt_pos2 = align1.f_ind[p2]
                wt_amino1 = f'{align1.record1.seq}'[wt_pos1]
                wt_amino2 = f'{align1.record1.seq}'[wt_pos2]
                
                if wt_pos1==wt_pos2:
                    indel_codes.append(f'{wt_amino1}{wt_pos1+1}del')
                else:
                    indel_codes.append(f'{wt_amino1}{wt_pos1+1}_{wt_amino2}{wt_pos2+1}del')
                    
    return substitution_codes, indel_codes


def mutation_codes_from_name(name, added_codes=[]):
    name = name[name.find('(')+1:]
    name = name[:name.rfind(')')]
    codes = name.split('/') + added_codes
    pos = [int(x[1:-1]) for x in codes]
    df = pd.DataFrame({'codes':codes, 'pos':pos})
    df.sort_values(by='pos', inplace=True)
    codes = list(df.codes)
    return codes

