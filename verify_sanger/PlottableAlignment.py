# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 19:14:38 2021

@author: djross
"""

from Bio.Align import PairwiseAlignment

class PlottableAlignment(PairwiseAlignment):
    def __init__(self, input_alignment):
        self.target = input_alignment.target
        self.query = input_alignment.query
        self.path = input_alignment.path
        self.score = input_alignment.score