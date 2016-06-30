'''Various functions for unsupervised kNN change detection. Includes both
a z-score and a modified z-score version.

Author:     Alex Lu 
            Computer Science Department, University of Toronto
            alexlu@cs.toronto.edu
Last Revised:     May 2nd, 2016, 4:53 PM EST
'''

import numpy as np
import sklearn.metrics.pairwise as skdist


def gaussianParams(k, geneMatrix, distMatrix):
    '''Generates a mean and variance vector for each gene in a gene matrix using the k closest genes using euclidean distance
    Inputs: normalization window k, gene matrix of size D x N to be normalized, matrix to calculate distances
    Output: mean and variance of k NN of genes'''
    dist = skdist.pairwise_distances(distMatrix.T, metric='euclidean')
    nearest = np.argsort(dist, axis=1)[:,1:(k+1)]
    
    means = np.zeros(geneMatrix.shape)
    for gene in range(0, nearest.shape[0]):
        for index in nearest[gene]:
            for feature in range(0, len(geneMatrix[index])):
                means[gene][feature] += geneMatrix[index][feature]
        means[gene] = means[gene] / k 
        if ((gene + 1) % 200 == 0 or (gene + 1) == geneMatrix.shape[0]):
            print ("Calculated means for %d out of %d genes." % ((gene + 1), geneMatrix.shape[0]))  
    
    variances = np.zeros(geneMatrix.shape)
    for gene in range(0, nearest.shape[0]):
        for index in nearest[gene]:
            for feature in range(0, len(geneMatrix[index])):
                variances[gene][feature] += np.square( (geneMatrix[index][feature] - means[index][feature]) )
        variances[gene] = variances[gene] / k 
        if ((gene + 1) % 200 == 0 or (gene + 1) == geneMatrix.shape[0]):
            print ("Calculated variances for %d out of %d genes." % ((gene + 1), geneMatrix.shape[0]))      
             
    return means, variances


def calculateZScores(geneMatrix, means, variances):
    '''Calculates z-score vectors for each gene given mean and variance vectors of kNN neighbors
    Inputs: gene matrix and corresponding mean and variances of kNN neighbors 
    Output: zscores'''
    zscores = np.zeros(geneMatrix.shape)  
    for gene in range(0, geneMatrix.shape[0]):
        for feature in range (0, len(geneMatrix[gene])):
            zscores[gene][feature] = (geneMatrix[gene][feature] - means[gene][feature]) / np.sqrt(variances[gene][feature])      
    return zscores


def modParams(k, geneMatrix, distMatrix):
    '''Generates a mean and MAD vector for each gene in a gene matrix using the k closest genes using euclidean distance
    Inputs: normalization window k, gene matrix of size D x N to be normalized, matrix to calculate distances
    Output: mean and MAD of k NN of genes'''
    dist = skdist.pairwise_distances(distMatrix.T, metric='euclidean')
    nearest = np.argsort(dist, axis=1)[:,1:(k+1)]
    
    means = np.zeros(geneMatrix.shape)
    for gene in range(0, nearest.shape[0]):
        for index in nearest[gene]:
            for feature in range(0, len(geneMatrix[index])):
                means[gene][feature] += geneMatrix[index][feature]
        means[gene] = means[gene] / k 
        if ((gene + 1) % 200 == 0 or (gene + 1) == geneMatrix.shape[0]):
            print ("Calculated means for %d out of %d genes." % ((gene + 1), geneMatrix.shape[0]))
    
    medians = np.zeros(geneMatrix.shape)
    MAD = np.zeros(geneMatrix.shape)
    for gene in range(0, nearest.shape[0]):
        neighbors = []
        for index in nearest[gene]:
            neighbors.append(geneMatrix[index])
            
        neighbors = np.array(neighbors)
        for feature in range(0, len(neighbors[0])):
            medians[gene][feature] = np.median(neighbors[:, feature])
            MAD[gene][feature] = np.median(abs(medians[gene][feature] - neighbors[:, feature]))
        if ((gene + 1) % 200 == 0 or (gene + 1) == geneMatrix.shape[0]):
            print ("Calculated MAD for %d out of %d genes." % ((gene + 1), geneMatrix.shape[0]))      
             
    return means, MAD

def calculateModZScores(geneMatrix, means, MAD):
    '''Calculates modified z-score vectors for each gene given mean and variance vectors of kNN neighbors
    Inputs: gene matrix and corresponding mean and variances of kNN neighbors 
    Output: zscores'''
    zscores = np.zeros(geneMatrix.shape)  
    for gene in range(0, geneMatrix.shape[0]):
        for feature in range (0, len(geneMatrix[gene])):
            zscores[gene][feature] = 0.6745 * (geneMatrix[gene][feature] - means[gene][feature]) / (MAD[gene][feature])      
    return zscores