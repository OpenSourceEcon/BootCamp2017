# svd_image_compression.py
"""Volume 1A: SVD and Image Compression.
<Name>
<Class>
<Date>
"""

from scipy import linalg as la
import numpy as np
from matplotlib import pyplot as plt

# Problem 1
def truncated_svd(A,k=None):
    """Computes the truncated SVD of A. If r is None or equals the number
        of nonzero singular values, it is the compact SVD.
    Parameters:
        A: the matrix
        k: the number of singular values to use
    Returns:
        U - the matrix U in the SVD
        s - the diagonals of Sigma in the SVD
        Vh - the matrix V^H in the SVD
    """
    raise NotImplementedError("truncated_svd incomplete")

# Problem 2
def visualize_svd():
    """Plot each transformation associated with the SVD of A."""
    raise NotImplementedError("visualize_svd incomplete")

# Problem 3
def svd_approx(A, k):
    """Returns best rank k approximation to A with respect to the induced 2-norm.

    Inputs:
    A - np.ndarray of size mxn
    k - rank

    Return:
    Ahat - the best rank k approximation
    """
    raise NotImplementedError("svd_approx incomplete")

# Problem 4
def lowest_rank_approx(A,e):
    """Returns the lowest rank approximation of A with error less than e
    with respect to the induced 2-norm.

    Inputs:
    A - np.ndarray of size mxn
    e - error

    Return:
    Ahat - the lowest rank approximation of A with error less than e.
    """
    raise NotImplementedError("lowest_rank_approx incomplete")

# Problem 5
def compress_image(filename,k):
    """Plot the original image found at 'filename' and the rank k approximation
    of the image found at 'filename.'

    filename - jpg image file path
    k - rank
    """
    raise NotImplementedError("compress_image incomplete")
