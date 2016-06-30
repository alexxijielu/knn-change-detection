'''Utility functions for the opening and writing of matrices of feature vectors.

Copyright notice: You are free to modify, extend and distribute 
this code for non-commercial purposes only, provided that the 
author of the original code is mentioned as the original author of the code.

Author:     Alex Lu 
            Computer Science Department, University of Toronto
            alexlu@cs.toronto.edu
Last Revised:     May 2nd, 2016, 4:53 PM EST
'''

import csv
import numpy as np


def openGeneMatrix (fileName):
    '''Opens a gene matrix file and returns the feature labels, gene labels and gene matrix
    Input: Path of file to be opened (as a string)
    Output: feature labels, gene labels, and gene matrix'''
    file = open(fileName)
    list = csv.reader(file, delimiter='\t')
    matrix = np.array([row for row in list])
    
    genelist = matrix[1:, 0]
    headers = matrix[0, :]
    genematrix = matrix[1:, 1:]
    genematrix = genematrix.astype(np.float)
    
    file.close()
    
    return headers, genelist, genematrix


def packageGeneMatrix (fileName, headers, genelist, genematrix):
    '''Combines the feature labels, gene labels and gene matrix and writes to output
    Input: Path of file to be written to, feature labels, gene labels, and gene matrix
    Output: Writes to path of file'''
    file = open(fileName, "w")
    file.truncate()
    
    file.write('\t'.join(headers) + '\n')
    for i in range (0, genematrix.shape[0]):
        file.write(genelist[i] + '\t')
        for j in genematrix[i]:
            file.write((str(j)) + '\t')
        file.write('\n')
    
    print("Written to file + " + fileName)
        
    file.close()