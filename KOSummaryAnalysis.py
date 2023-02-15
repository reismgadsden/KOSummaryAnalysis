"""
KOSummaryAnalysis.py

This program will analyze the annotated .tsv/.csv as well as the
Trinity output file in order to pull the following statistics for
each KO
- Total number of transcripts
- Total number of unique gene ids
- Average of log fold change
- Minimum of log fold change
- Maximum of log fold change
- P-value of the highest log fold change

author - Reis Mercer Gadsden
version - 2023/02/01
institution - Appalachian State University

github -
"""

"""""""""""
" IMPORTS "
"""""""""""
import pandas
import numpy
import os


class KOSummaryAnalysis:

    def __init__(self, input_file=None):
        if input_file is None:
            print("Please provide an input file!")