# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 19:42:27 2021

@author: djross
"""

import numbers
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# sanger_channels correspond to the ABI data traces for the letters: G, A, T, C
sanger_channels = ["DATA9", "DATA10", "DATA11", "DATA12"]


class PlottableRecord(SeqRecord):
    def __init__(self, input_record):
        super().__init__(
            seq=input_record.seq,
            id=input_record.id,
            name=input_record.name,
            description=input_record.description,
            dbxrefs=input_record.dbxrefs,
            features=input_record.features,
            annotations=input_record.annotations,
            letter_annotations=input_record.letter_annotations,
            )
        
        self.init_chrom_plot_data(input_record)
        self.coverage = [0 if s == '-' else 1 for s in self.seq]
        
    def init_chrom_plot_data(self, input_record):
        if 'abif_raw' in input_record.annotations.keys():
            raw_data = input_record.annotations["abif_raw"]
            channel_data = [raw_data[x] for x in sanger_channels]
            peak_locations = np.array(raw_data['PLOC1'])
            peak_left, peak_right = peak_edges(peak_locations)
            
            chrom_x = []
            chrom_g = []
            chrom_a = []
            chrom_t = []
            chrom_c = []
            for i, (left, right) in enumerate(zip(peak_left, peak_right)):
                for y_data, y_list in zip(channel_data, [chrom_g, chrom_a, 
                                                         chrom_t, chrom_c]):
                    y = y_data[left:right+1]
                    y_list.append(np.array(y))
                x = np.linspace(i+0.5, i+1.5, len(y))
                chrom_x.append(x)
                
            chrom_x = np.array(chrom_x)
            chrom_g = np.array(chrom_g)
            chrom_a = np.array(chrom_a)
            chrom_t = np.array(chrom_t)
            chrom_c = np.array(chrom_c)
            
            self._chrom_data = [chrom_x, chrom_g, chrom_a, 
                                chrom_t, chrom_c]
        else:
            self._chrom_data = [np.array([np.nan]*len(self))]
            
    
    def __getitem__(self, index):
        # If index is an integer, Instead of super() behavior (return single letter),
        #     return PlottableRecord object with a sequence of length 1.
        if isinstance(index, numbers.Integral):
            index = slice(index, index+1)
        new_record = SeqRecord(
            seq=self.seq,
            id=self.id,
            name=self.name,
            description=self.description,
            dbxrefs=self.dbxrefs,
            features=self.features,
            annotations=self.annotations,
            letter_annotations=self.letter_annotations,
            )
        new_record = new_record[index.start:index.stop:index.step]
        new_record = PlottableRecord(new_record)
        
        # Slicing the chromatogram data:
        new_data = [x[index.start:index.stop:index.step] for x in self._chrom_data]
        x_shift = [x[0]-i-0.5 for i, x in enumerate(new_data[0])]
        new_data[0] = new_data[0] - x_shift
        
        new_record._chrom_data = new_data
        
        return new_record
    
    
    def reverse_complement(self):
        new_record = super().reverse_complement()
        
        new_record = PlottableRecord(new_record)
        
        new_record.coverage = self.coverage[::-1]
        
        new_data = []
        
        data = self._chrom_data[0][::-1]
        x_shift = [x[0]-i-0.5 for i, x in enumerate(data)]
        new_data.append(data - x_shift)
        
        for data in self._chrom_data[:0:-1]:
            data = [d[::-1] for d in data[::-1]]
            new_data.append(data)
            
        new_record._chrom_data = new_data
        
        return new_record
    
    
    @classmethod
    def make_aligned_record(cls, input_record, aligned_seq):
        """
        The method is effectively an alternate constructor for making a new
        PlottableRecord instance based on a sequence alignment

        Parameters
        ----------
        input_record : sequence record, PlotttableRecord or Bio.SeqRecord.SeqRecord
            The sequence record used as input to the alignment.
        aligned_seq : DNA sequence, Bio.Seq.Seq or str
            The aligned sequence, with possible gaps.

        Raises
        ------
        ValueError
            aligned_seq must result from an alignment of input_record. So, 
            after ungapping, aligned_seq must be the same sequence as input_record.

        Returns
        -------
        new_record : PlottableRecord
            The new plottable sequence record, adjusted to match the alignment
            (i.e. gaps insterted into the quality and chromatogram data, etc.).

        """
        
        aligned_seq = str(aligned_seq)
        input_seq = str(input_record.seq)
        
        # Check that align_str contains the same sequence as input_record, with possible gaps
        if aligned_seq.replace('-', '') != input_seq.replace('-', ''):
            raise ValueError("After ungapping, aligned_seq must be the same sequence as input_record.")
        
        new_record = SeqRecord(seq=Seq(aligned_seq),
                               id=input_record.id,
                               name=input_record.name,
                               description=input_record.description,
                               dbxrefs=input_record.dbxrefs,
                               features=input_record.features,
                               annotations=input_record.annotations,
                               letter_annotations=None,
                               )
        
        new_record = PlottableRecord(new_record)
        
        # Loop over all associated lists/arrays that have an entry for 
        #     each base in the sequence
        for k, v in input_record.letter_annotations.items():
            new_record.letter_annotations[k] = align_data(v, input_seq, aligned_seq)
        
        #new_record.coverage = align_data(new_record.coverage, input_seq, aligned_seq, align_str)
        
        new_data = []
        
        data = input_record._chrom_data[0]
        data = np.array( align_data(data, input_seq, aligned_seq, fill_item=np.array([0, 1])) )
        x_shift = [x[0]-i-0.5 for i, x in enumerate(data)]
        new_data.append(data - x_shift)
        
        for data in input_record._chrom_data[1:]:
            data = np.array( align_data(data, input_seq, aligned_seq, fill_item=np.array([np.nan, np.nan])) )
            new_data.append(data)
            
        new_record._chrom_data = new_data
        
        new_record.seq = Seq(aligned_seq)
        
        return new_record
    
    
    def chromatogram_plot_data(self, start=0, end=None):
        if end == None:
            end = len(self)
            
        out_data = [x[start:end] for x in self._chrom_data]
        out_data = [np.array([item for sublist in data for item in sublist]) for data in out_data]
        return out_data
        
    
    
    def quality_plot_data(self, start=0, end=None):
        if end == None:
            end = len(self)
        
        qual_x = qual_cover_x(start, end)
        
        qual = self.letter_annotations['phred_quality']
        qual_y = np.array( [ np.array([x]*2) for x in qual[start:end] ] ).flatten()
        qual_y = np.append(qual_y, 0)
        qual_y = np.insert(qual_y, 0, 0)
        
        return [qual_x, qual_y]
    
    
    def coverage_plot_data(self, start=0, end=None):
        if end == None:
            end = len(self)
        
        cover_x = qual_cover_x(start, end)
        
        cover_y = np.array( [ np.array([x]*2) for x in self.coverage[start:end] ] ).flatten()
        cover_y = np.append(cover_y, 0)
        cover_y = np.insert(cover_y, 0, 0)
        
        return [cover_x, cover_y]
    
    
def qual_cover_x(start, end):
        x_data = np.array( [ np.array([x + 0.5]*2) for x in range(start, end) ] ).flatten()
        x_data = np.append(x_data, [end+0.5]*2)
        
        return x_data
    
def peak_edges(peak_locations):
    """
    Method to get the left and right edges of chromatogram peaks, based on the
    input array of peak centers

    Parameters
    ----------
    peak_locations : np.array
        1D array of positions of peak centers in Sagner chromatogram data.

    Returns
    -------
    peak_left : np.array
        1D array of positions of left edges of Sanger chromatogram peaks
    peak_right : np.array
        1D array of positions of right edges of Sanger chromatogram peaks

    """
    
    peak_left = peak_locations
    peak_left = (peak_left[:-1] + peak_left[1:])/2
    peak_right = np.append( peak_left, peak_locations[-1] + 1/2*(peak_left[-1] - peak_left[-2]) )
    peak_left = np.insert( peak_left, 0, peak_locations[0] - 1/2*(peak_left[1] - peak_left[0]) )
    peak_left = np.array([int(np.round(x)) if x>0 else 0 for x in peak_left])
    peak_right = np.array([int(np.round(x)) if x>0 else 0 for x in peak_right])
    
    return peak_left, peak_right


def align_data(data, input_seq, aligned_seq, fill_item=0):
    """
    Method to create new lists of data (letter annotations, chromatogram data)
        based on an alignment (i.e. insterting gaps, etc.)

    Parameters
    ----------
    data : array like
        the data to be aligned.
    input_seq : string
        The DNA sequence of the input record used for the alignment.
    aligned_seq : string
        The aligned sequence, with possible gaps.

    Returns
    -------
    new_data : list
        the resulting aligned data.

    """
    
    if len(data) != len(input_seq):
        raise ValueError(f"The length of data ({len(data)}) must equal the length of input_seq ({len(input_seq)})")
        
    new_data = []
    for b1 in aligned_seq:
        # gap
        if b1 == '-':
            new_data.append(fill_item)
        else:
            if len(data) == 0:
                raise ValueError(f'Unexpected non-gap base found in aligned_seq ({b1}) after end of data')
            elif b1 == input_seq[0]:
                new_data.append(data[0])
                data = data[1:]
                input_seq = input_seq[1:]
            else:
                raise ValueError(f'Unexpected missmatch between bases in aligned_seq ({b1}) and input_seq ({input_seq[0]})')
    
    return new_data
