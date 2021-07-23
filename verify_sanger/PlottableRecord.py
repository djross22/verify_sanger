# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 19:42:27 2021

@author: djross
"""

import numpy as np
from Bio.SeqRecord import SeqRecord

# sanger_channels correspond to the ABI data traces for the letters: G, A, T, C
sanger_channels = ["DATA9", "DATA10", "DATA11", "DATA12"]

class PlottableRecord(SeqRecord):
    def __init__(self, input_record):
        self.description = input_record.description
        self.id = input_record.id
        self.name = input_record.name
        self._seq = input_record.seq
        
        self._per_letter_annotations = input_record._per_letter_annotations
        
        #for k in input_record.letter_annotations.keys():
        #    self.letter_annotations[k] = input_record.letter_annotations[k]
        
        self.init_chrom_plot_data(input_record)
        
    def init_chrom_plot_data(self, input_record):
        if 'abif_raw' in input_record.annotations.keys():
            raw_data = input_record.annotations["abif_raw"]
            channel_data = [raw_data[x] for x in sanger_channels]
            peak_locations = np.array(raw_data['PLOC1'])
            peak_left, peak_right = peak_edges(peak_locations)
            
            self._chrom_x = []
            self._chrom_g = []
            self._chrom_a = []
            self._chrom_t = []
            self._chrom_c = []
            for i, (left, right) in enumerate(zip(peak_left, peak_right)):
                for y_data, y_list in zip(channel_data, [self._chrom_g, self._chrom_a, 
                                                         self._chrom_t, self._chrom_c]):
                    y = y_data[left:right+1]
                    y_list.append(np.array(y))
                x = np.linspace(i+0.5, i+1.5, len(y))
                self._chrom_x.append(x)
                
            self._chrom_x = np.array(self._chrom_x)
            self._chrom_g = np.array(self._chrom_g)
            self._chrom_a = np.array(self._chrom_a)
            self._chrom_t = np.array(self._chrom_t)
            self._chrom_c = np.array(self._chrom_c)
            
            self._chrom_data = [self._chrom_x, self._chrom_g, self._chrom_a, 
                                self._chrom_t, self._chrom_c]
            
    
    def chromatogram_plot_data(self, start=0, end=None):
        if end == None:
            end = len(self._chrom_data[0])
            
        out_data = [x[start:end] for x in self._chrom_data]
        out_data = [np.array([item for sublist in data for item in sublist]) for data in out_data]
        return out_data
        
    
    
    def quality_plot_data(self, start=0, end=None):
        if end == None:
            end = len(self._chrom_data[0])
        
        qual_x = qual_cover_x(start, end)
        
        qual = self.letter_annotations['phred_quality']
        qual_y = np.array( [ np.array([x]*2) for x in qual[start:end] ] ).flatten()
        qual_y = np.append(qual_y, 0)
        qual_y = np.insert(qual_y, 0, 0)
        
        return [qual_x, qual_y]
    
    
    def coverage_plot_data(self, start=0, end=None):
        if end == None:
            end = len(self._chrom_data[0])
        
        cover_x = qual_cover_x(start, end)
        cover = [0 if s == '-' else 1 for s in self.seq]
        
        cover_y = np.array( [ np.array([x]*2) for x in cover[start:end] ] ).flatten()
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
