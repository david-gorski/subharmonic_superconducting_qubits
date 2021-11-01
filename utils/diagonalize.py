import numpy as np

def diagonalize(matrix):
    # clearly only works if conditions for diagonalization are met
    # matrix = np.array(matrix)
    w,v = np.linalg.eig(matrix)
    # D = np.reshape(np.array([vec for vec in v]), matrix.shape)
    diagonalized_matrix = np.matrix(v).getH() * np.matrix(matrix) * np.matrix(v)
    return diagonalized_matrix