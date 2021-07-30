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
        
        align_str = self.align_str
        self.record1 = PlottableRecord.make_aligned_record(record1, align_str[0])
        self.record2 = PlottableRecord.make_aligned_record(record2, align_str[2])
    
    @property
    def align_str(self):
        return f'{self}'.split('\n')
        