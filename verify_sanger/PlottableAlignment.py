# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 19:14:38 2021

@author: djross
"""

from Bio.Align import PairwiseAlignment

class PlottableAlignment(PairwiseAlignment):
    def __init__(self, input_alignment):
        super().__init__(target=input_alignment.target, 
                         query=input_alignment.query, 
                         path=input_alignment.path, 
                         score=input_alignment.score,
                         )