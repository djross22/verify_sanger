# -*- coding: utf-8 -*-
"""
Created on Mon May 10 2021

@author: djross
"""

from setuptools import setup, find_packages

setup(name='verify_sanger',
      version='1.0',
      description='Python package for automated verification of genes inserted into a plasmid or other construct using Sanger sequencing data (.ab1 files) and a reference sequence',
      packages=find_packages(),
      zip_safe=False)