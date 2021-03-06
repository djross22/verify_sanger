# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 19:14:38 2021

@author: djross
"""

from Bio.Align import PairwiseAlignment
from verify_sanger.PlottableRecord import PlottableRecord 

class PlottableAlignment(PairwiseAlignment):
    def __init__(self, input_alignment,
                 consensus_seq,
                 mismatch_ind,
                 record1,
                 record2):
        
        super().__init__(target=input_alignment.target, 
                         query=input_alignment.query, 
                         path=input_alignment.path, 
                         score=input_alignment.score,
                         )
        self.consensus_seq = consensus_seq
        self.mismatch_ind = mismatch_ind
        
        # Initially set record1 and record2 to None (so they exist), 
        #     then replace them a few lines down
        self.record1 = None
        self.record2 = None
        
        # Make sure inputs, record1 and record2, are PlottableRecords
        if not isinstance(record1, PlottableRecord):
            record1 = PlottableRecord(record1)
        if not isinstance(record2, PlottableRecord):
            record2 = PlottableRecord(record2)
        
        self.record1 = PlottableRecord.make_aligned_record(self, record1, target=True)
        self.record2 = PlottableRecord.make_aligned_record(self, record2, target=False)
    
    @property
    def align_str(self):
        return f'{self}'.split('\n')
    
    
    def map_coordinate(self, input_index, target=True):
        """
        maps a location in one of the input sequences used for the alignment
        to a locaiton in the corresponding output sequence

        Parameters
        ----------
        input_index : int
            The position within the input sequence to be mapped
        target : Boolean, optional
            If True, the input sequence is self.target; 
            if False, the input sequence is self.query. 
            The default is True.

        Returns
        -------
        output_index : int
            The index in the output sequence that aligns with 
            input_indes in the input sequence.

        """
        input_seq = self.target if target else self.query
        if target:
            if self.record1 is not None:
                output_seq = self.record1.seq
            else:
                output_seq = self.align_str[0]
        else:
            if self.record2 is not None:
                output_seq = self.record2.seq
            else:
                output_seq = self.align_str[2]
                
        return map_index_across_gaps(input_index, input_seq, output_seq)
    
    
def map_index_across_gaps(input_index, input_seq, output_seq):
    
    # Make sure the input_index is withing the input sequence
    if input_index < 0:
        raise ValueError("The input_index ({input_index}) must be non-negative")
    elif input_index > len(input_seq):
        raise ValueError(f"The input_index ({input_index}) must less than ar equal to the length of the input sequence ({len(input_seq)})")
    
    # After ungapping (removing dashes), the input_seq and output_seq should be the same
    if str(input_seq).replace('-', '') != str(output_seq).replace('-', ''):
        raise ValueError("After ungapping (removing dashes), the input_seq and output_seq must be the same")
        
    # If input_seq[input_index] is a gap, find last non-gap and return corresponding index+1
    input_gap_correction = 0
    if input_index < len(input_seq):
        while input_seq[input_index] == '-':
            if input_index == 0:
                return 0
            input_index -= 1
            input_gap_correction = 1
            
    
    i = 0
    j = 0
    while (i < len(input_seq)) and (j < len(output_seq)):
        b1 = input_seq[i]
        b2 = output_seq[j]
        if b1 == b2:
            if i == input_index:
                break
            else:
                i += 1
                j += 1
        elif b1 == '-':
            i += 1
        elif b2 == '-':
            j += 1
        else:
            # If there is a missmatch, one of the bases should be a dash; if not, something weird has happended 
            #TODO: raise an error here?
            print(f'Unexpected base missmatch in map_index_across_gaps(input_index={input_index})')
    
    j += input_gap_correction
        
    return j
        
        