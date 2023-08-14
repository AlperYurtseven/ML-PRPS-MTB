import pandas as pd
import numpy as np


def tetrachoric(x, y):
    """Calculate the tetrachoric correlation coefficient between two dichotomous variables.
       Args:
           x (array-like): An array of dichotomous data.
           y (array-like): An array of dichotomous data with the same length as x.

       Returns:
           float: The tetrachoric correlation coefficient between x and y.

       """
    a = np.sum(y[x == 1])  # number of y=1 given x=1
    b = len(x) - np.sum(y[x == 1])  # number of y=0 given x=1
    c = np.sum(y[x == 0])  # number of y=1 given x=0
    d = len(x) - np.sum(y[x == 0])  # number of y=0 given x=0
    if b == 0:
        b = 1
    if c == 0:
        c = 1
    odds = a * d / (b * c)
    corr = np.cos(np.pi / (1 + np.sqrt(odds)))
    return corr


def tetrachoric_matrix(data):
    """Create a correlation matrix of tetrachoric correlations for a DataFrame of dichotomous data.

        This function takes a pandas DataFrame where each column is an array of dichotomous data and calculates
        the tetrachoric correlation coefficient between all pairs of columns in the DataFrame. The resulting
        correlation matrix is a pandas DataFrame with variable names as column and row indices.

        Args:
            data (pandas.DataFrame): A DataFrame of dichotomous data. Indices should be the identifiers

        Returns:
            pandas.DataFrame: A correlation matrix of tetrachoric correlations for the input data.

        """

    corr_matrix = np.eye((data.shape[1]))
    row_index, column_index = np.tril_indices(corr_matrix.shape[0], -1)

    for row_ndx, col_ndx in zip(row_index, column_index):
        corr_matrix[row_ndx, col_ndx] = tetrachoric(data.iloc[:, row_ndx],
                                                    data.iloc[:, col_ndx])
        corr_matrix[col_ndx, row_ndx] = corr_matrix[row_ndx, col_ndx]

    # corr_matrix=np.abs(corr_matrix)
    corr_matrix = pd.DataFrame(corr_matrix, index=data.columns, columns=data.columns)

    return corr_matrix

# Usage example

#feature_df=pd.read_csv("all_features_tsv_file", sep="\t", usecols=[list_of_features_of_interest])
#feat_matrix=tetrachoric_matrix(feature_df)
#cmap=sns.cm.vlag
#sns.clustermap(feat_matrix, cmap=cmap)