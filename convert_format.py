#!/usr/bin/env python3

from __future__ import division
import pandas as pd
import numpy as np
import csv
import os
import sys


def binary_to_sparse_format(X, Y, file):

    """
    Converts a samples times SNPs matrix and the corresponding samples times
    tasks response matrix to the sparse format required by the Safe Pattern
    Pruning algorithm by Nakagawa et al. and saves it as a csv file.
    (https://github.com/takeuchi-lab/SafePatternPruning/blob/master/itemLasso/data/dna)

    Parameters
    ----------
    file: str
        Name of the output csv file written by the function

    X: (n x p) pandas DataFrame or (n x p) numpy array
        samples times SNPs matrix

    Y: (n x T) pandas DataFrame or (n x T) numpy array
        Y must be the same type as X.
        samples times tasks response matrix

    """

    with open(file, "w") as fp:
        wr = csv.writer(fp, delimiter='\t')
        if isinstance(X, pd.core.frame.DataFrame):
            for i in range(X.shape[0]):
                tmp = [str(j + 1) + ':' + str(X.iloc[i, j])
                       for j in range(X.shape[1])
                       if X.iloc[i, j] != 0]
                tmp = Y.iloc[i, :].tolist() + tmp
                wr.writerow(tmp)
        if isinstance(X, np.ndarray):
            for i in range(X.shape[0]):
                tmp = [str(j + 1) + ':' + str(X[i, j])
                       for j in range(X.shape[1])
                       if X[i, j] != 0]
                tmp.insert(0, Y[i])
                wr.writerow(tmp)

x_file = sys.argv[1]
y_file = sys.argv[2]

output = sys.argv[3]

testX = np.genfromtxt(x_file, delimiter=",")[:,1:]
testY = np.genfromtxt(y_file, delimiter=",")[:,1:]
testY = testY.reshape(np.shape(testY)[0])

binary_to_sparse_format(testX, testY, output)

#def func(lam, *args):
#    return(np.exp(lam) / (np.exp(lam) - 1) - 1 / lam - args[0])
#
## Folder where the simulated data is saved
#path_to_data = "./data_preprocessed/"
#createFolder(path_to_data)
## Folder where figures are saved
#path_to_figures = "./figures/"
#createFolder(path_to_figures)
#
#
## Simulate datasets with iid variables drawn from a Bernoulli distribution
## with parameter q itself drawn from a uniform distribution on [0.1, 0.5]
## and columns ordered in decreasing order of vector size.
## (data is written in the format required by SPP)
## ------------------------------------------------------------------------
#rs = 0
#np.random.seed(rs)
#
#dims = [(300, 1000), (1000, 1000), (1000, 3000), (1000, 10000), (10000, 1000)]
#
#coef = np.sort(np.random.normal(0, 1, size=100))[::-1]
#ind = np.random.choice(range(1000), size=(len(coef), 2))
#
#for (n, p) in dims:
#    print n, p
#    X = np.zeros((n, p))
#    for j in range(p):
#        q = np.random.uniform(0.1, 0.5)
#        X[:, j] = np.random.binomial(n=1, p=q, size=n)
#    print "Bernoulli drawn"
#    # Order columns by decreasing order of vector size
#    X = X[:, np.argsort(-X.sum(axis=0))]
#    y = np.zeros((n))
#    for i in range(ind.shape[0]):
#        j = ind[i, 0]
#        k = ind[i, 1]
#        y += coef[i]*X[:, j]*X[:, k]
#    binary_to_sparse_format(
#        X=X, Y=y,
#        file=path_to_data + "Bernoulli_n" + str(n) + "_p" + str(p) +
#        "_qunif_coefnormal_rs" + str(rs) + "_nnzd.tsv")
#    # Write the index and coefficients of the ground truth features in another
#    # file
#    gt = np.concatenate([ind, coef.reshape((len(coef), 1))], axis=1)
#    np.savetxt(path_to_data + "GroundTruth_n" + str(n) + "_p" + str(p) +
#               "_qunif_coefnormal_rs" + str(rs) + "_nnzd.csv", X=gt,
#               delimiter=",")
#    print "file written"
#
#
## Simulate datasets with iid variables drawn from a Bernoulli distribution
## with parameter q drawn from a uniform distribution on [0.1, 0.5]
## and y defined by its cumsum which is itself parametrized by an exponential
## and columns ordered in decreasing order of vector size.
## (data is written in the format required by SPP)
## ------------------------------------------------------------------------
#rs = 0
#
## Choose the dimension of the simulated dataset
#dims = [(1000, 1000)]
## Choose the area under the cumulative sul curve of the response
#l_auc = [0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]
#
#l_lam = [fsolve(func, 1, args=(auc)) for auc in l_auc]
#l_handles = []
#colors = cm.rainbow(np.linspace(0, 1, len(l_auc)))
#plt.figure(figsize=(5.4, 4.4))
#for (n, p) in dims:
#    print n, p
#    X = np.zeros((n, p))
#    np.random.seed(rs)
#    for j in range(p):
#        q = np.random.uniform(0.1, 0.5)
#        X[:, j] = np.random.binomial(n=1, p=q, size=n)
#    print "Bernoulli drawn"
#    # Order columns by decreasing order of vector size
#    X = X[:, np.argsort(-X.sum(axis=0))]
#    for (a, auc) in enumerate(l_auc):
#        lam = l_lam[a]
#        cumsum = [1 / (1 - np.exp(-lam)) * (1 - np.exp(-lam * i))
#                  for i in np.arange(0, 1, 1. / ((n + 1) / 2))]
#        h, = plt.plot(cumsum, color=colors[a])
#        l_handles.append(h)
#        yp = np.array([cumsum[i] - cumsum[i - 1]
#                      for i in range(1, int(n / 2) + 1)])
#        y = np.append(yp, -yp)
#        binary_to_sparse_format(
#            X=X, Y=y,
#            file=path_to_data + "Bernoulli_n" + str(n) + "_p" + str(p) +
#            "_qunif_yexp_auc" + str(auc) + "_rs" + str(rs) + "_nnzd.tsv")
#        print "file written"
#plt.legend(l_handles, l_auc, title=r'$\kappa$', loc=4, fontsize=13,
#           frameon=False)
#plt.xlabel("Coordinate", fontsize=14)
#plt.ylabel("Cumulative sum", fontsize=14)
#plt.tick_params(axis='both', labelsize=14)
#plt.savefig(path_to_figures + "response_cumsum.pdf", format='pdf',
#            bbox_inches='tight')
