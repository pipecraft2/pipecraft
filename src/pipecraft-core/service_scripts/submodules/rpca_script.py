import sys
import os
import biom
import skbio
import numpy as np
import pandas as pd
from typing import Union
from deicode.matrix_completion import MatrixCompletion
from deicode.preprocessing import rclr
from scipy.linalg import svd

## This is a script for custom RPCA
## Input = RCLR-transformed matrix (rows = samples, cols = OTUs)



## Read the command-line arguments
input_file = sys.argv[1]   # input_file = "rclr.tsv"

## Output directory
try:
  output_dir = sys.argv[2]
except IndexError:
  output_dir = "Res"

## Opitonal arguments: number of PCA components and OptSpace iterations
## currently not implemented (using `auto_rpca`)
try:
  DEFAULT_RANK = int(sys.argv[3])
except IndexError:
  DEFAULT_RANK = int(3)

try:
  DEFAULT_ITERATIONS = int(sys.argv[4])
except IndexError:
  DEFAULT_ITERATIONS = int(5)



## Load input data
table = pd.io.parsers.read_csv(input_file, sep="\t", index_col=0)

## Convert class 'pandas.core.frame.DataFrame' into 'numpy.ndarray'
NDARRAY = table.select_dtypes(include = float).to_numpy()

## Preview input table
# print(type(table))
# print(table)

## Preview numpy array
# print(type(NDARRAY))
# print(NDARRAY)

def rpca(table, #biom.Table,
         n_components: Union[int, str] = DEFAULT_RANK,
         max_iterations: int = DEFAULT_ITERATIONS) -> (
        skbio.OrdinationResults,
        skbio.DistanceMatrix):
    """
    Runs RPCA on a constructed RCLR matrix
    """

    # get shape of table
    # n_features, n_samples = table.shape

    # OptSpace (RPCA) on robust-CLR (rclr) preprocessed data
    opt = MatrixCompletion(n_components=n_components,
                           max_iterations=max_iterations).fit(NDARRAY)

    # get new n-comp when applicable
    n_components = opt.s.shape[0]
    # get PC column labels for the skbio OrdinationResults
    rename_cols = ['PC' + str(i + 1) for i in range(n_components)]
    # get completed matrix for centering
    X = opt.sample_weights @ opt.s @ opt.feature_weights.T
    # center again around zero after completion
    X = X - X.mean(axis=0)
    X = X - X.mean(axis=1).reshape(-1, 1)
    # re-factor the data
    u, s, v = svd(X)
    # only take n-components
    u = u[:, :n_components]
    v = v.T[:, :n_components]
    # calc. the new variance using projection
    p = s**2 / np.sum(s**2)
    p = p[:n_components]
    s = s[:n_components]
    # save the loadings
    feature_loading = pd.DataFrame(v, index=table.columns, columns=rename_cols)
    sample_loading = pd.DataFrame(u, index=table.index, columns=rename_cols)
    # % var explained
    proportion_explained = pd.Series(p, index=rename_cols)
    # get eigenvalues
    eigvals = pd.Series(s, index=rename_cols)

    # if the n_components is two add PC3 of zeros
    # this is referenced as in issue in
    # <https://github.com/biocore/emperor/commit
    # /a93f029548c421cb0ba365b4294f7a5a6b0209ce>
    # discussed in DEICODE -- PR#29
    if n_components == 2:
        feature_loading['PC3'] = [0] * len(feature_loading.index)
        sample_loading['PC3'] = [0] * len(sample_loading.index)
        eigvals.loc['PC3'] = 0
        proportion_explained.loc['PC3'] = 0

    # save ordination results
    short_method_name = 'rpca_biplot'
    long_method_name = '(Robust Aitchison) RPCA Biplot'
    ord_res = skbio.OrdinationResults(
        short_method_name,
        long_method_name,
        eigvals.copy(),
        samples=sample_loading.copy(),
        features=feature_loading.copy(),
        proportion_explained=proportion_explained.copy())
    
    # save distance matrix
    dist_res = skbio.stats.distance.DistanceMatrix(
        opt.distance, ids=sample_loading.index)

    return ord_res, dist_res


def auto_rpca(table, # biom.Table,
              max_iterations: int = DEFAULT_ITERATIONS) -> (
        skbio.OrdinationResults,
        skbio.DistanceMatrix):
    """Runs RPCA but with auto estimation of the
       rank peramater.
    """
    ord_res, dist_res = rpca(table, n_components='auto')
    return ord_res, dist_res





ord_res, dist_res = auto_rpca(table)

os.makedirs(output_dir, exist_ok=True)

ord_res.write(os.path.join(output_dir, 'ordination.txt'))
dist_res.write(os.path.join(output_dir, 'distance-matrix.tsv'))

