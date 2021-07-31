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
        output_seq = self.record1.seq if target else self.record2.seq
        
        # Make sure the input_index is withing the input sequence
        if input_index < 0:
            raise ValueError("The input_index ({input_index}) must be non-negative")
        elif input_index >= len(input_seq):
            raise ValueError(f"The input_index ({input_index}) must less than the length of the input sequence ({len(input_seq)})")
        
        i = 0
        for j, b in enumerate(output_seq):
            if b == input_seq[i]:
                if i == input_index:
                    break
                else:
                    i += 1
            elif b != '-':
                # Each base, b, should match or be a dash; if not, something weird has happended 
                print(f'Unexpected base missmatch in map_coordinate(input_index={input_index}, target={target})')
        
        return j
        