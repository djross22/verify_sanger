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

def mott_trimming_fr(record, limit=0.05):
    """
    This method trims the low-quality reads from the ends of each read
    using Mott's algorithm in both directions and returns the idices
    that indicate the most conservative trimming on each end
    
    Parameters
    ----------
    record : Biopython sequence object to be trimmed
        
    limit : float
        A value between 0 and 1 indicating the error probability thershhold
        used for the sequence trimming

    Returns
    -------
    the indices for the first and last sequence element that should be kept
    (or np.nan if the sequence quality never exceeds the threshold)
        
    """
    
    # Forward direction
    qual = np.array(record.letter_annotations['phred_quality'])
    score_list = limit - (10 ** (qual / -10.0))
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
    score_list = limit - (10 ** (qual / -10.0))
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
        trim_ind1 = mott_trimming_fr(record1, limit=trim)
        trim_ind2 = mott_trimming_fr(record2, limit=trim)
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

