Information on the files in this repository:

SOFTWARE:
knnChangeDetect.py: Python functions for kNN change detection. Requires numpy and scikit learn, as well as util.py.
util.py: Utility functions for unpackaging tab-delimited gene feature files. Requires numpy.
example.py: Example usage of the above two packages.

DATA:
Profiles.xlsx: Per-gene features (mean, trunucated mean, and median profiling) from single cell data (not included due to size limitations) extracted from the WT3 and RPD3_1 datasets in the CYCLoPs database using image analysis software by Handfield et al. 2015. To use these with the software provided, save individual pages of the spreadsheet as tab-delimited text files.
z_scores.xlsx: Processed data (modified z-scores) derived from features in Profiles.xlsx using knn change detection method.