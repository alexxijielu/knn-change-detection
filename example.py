import knnChangeDetect
import util

#unpackage the headers, genes, and features for the change matrix and the WT matrix
headers, genelist, genematrix = util.openGeneMatrix("\Delta.txt")
wtheaders, wtgenelist, wtgenematrix = util.openGeneMatrix("\WT.txt")

#generate mean and MAD for each gene, and then the modified z-scores
means, MAD = knnChangeDetect.modWeights(50, genematrix, wtgenematrix)
zscores = knnChangeDetect.calculateModZScores(genematrix, means, variances)

#print results
for index in range(0, zscores.shape[0]):
    print genelist[index], 
    print "\t",
    for score in zscores[index]:
        print score,
        print "\t",
    print ""