# qr_lstsq_eigs.py
"""Volume 1A: QR 2 (Least Squares and Computing Eigenvalues).
<Name>
<Class>
<Date>"""

import numpy as np
from cmath import sqrt
from scipy import linalg as la
from matplotlib import pyplot as plt


# Problem 1
def least_squares(A, b):
    """Calculate the least squares solutions to Ax = b using QR decomposition.

    Inputs:
        A ((m,n) ndarray): A matrix of rank n <= m.
        b ((m, ) ndarray): A vector of length m.

    Returns:
        x ((n, ) ndarray): The solution to the normal equation.
    """
    raise NotImplementedError("Problem 1 Incomplete")


# Problem 2
def line_fit():
    """Load the data from housing.npy. Use least squares to calculate the line
    that best relates height to weight.

    Plot the original data points and the least squares line together.
    """
    raise NotImplementedError("Problem 2 Incomplete")


# Problem 3
def polynomial_fit():
    """Load the data from housing.npy. Use least squares to calculate
    the polynomials of degree 3, 6, 9, and 12 that best fit the data.

    Plot the original data points and each least squares polynomial together
    in individual subplots.
    """
    raise NotImplementedError("Problem 3 Incomplete")


def plot_ellipse(a, b, c, d, e):
    """Plot an ellipse of the form ax^2 + bx + cxy + dy + ey^2 = 1."""
    theta = np.linspace(0, 2*np.pi, 200)
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    A = a*(cos_t**2) + c*cos_t*sin_t + e*(sin_t**2)
    B = b*cos_t + d*sin_t
    r = (-B + np.sqrt(B**2 + 4*A))/(2*A)

    plt.plot(r*cos_t, r*sin_t, lw=2)
    plt.gca().set_aspect("equal", "datalim")

# Problem 4
def ellipse_fit():
    """Load the data from ellipse.npy. Use least squares to calculate the
    ellipse that best fits the data.

    Plot the original data points and the least squares ellipse together.
    """
    raise NotImplementedError("Problem 4 Incomplete")


# Problem 5
def power_method(A, N=20, tol=1e-12):
    """Compute the dominant eigenvalue of A and a corresponding eigenvector
    via the power method.

    Inputs:
        A ((n,n) ndarray): A square matrix.
        N (int): The maximum number of iterations.
        tol (float): The stopping tolerance.

    Returns:
        (foat): The dominant eigenvalue of A.
        ((n, ) ndarray): An eigenvector corresponding to the dominant
            eigenvalue of A.
    """
    raise NotImplementedError("Problem 5 Incomplete")


# Problem 6
def qr_algorithm(A, N=50, tol=1e-12):
    """Compute the eigenvalues of A via the QR algorithm.

    Inputs:
        A ((n,n) ndarray): A square matrix.
        N (int): The number of iterations to run the QR algorithm.
        tol (float): The threshold value for determining if a diagonal block
            is 1x1 or 2x2.

    Returns:
        ((n, ) ndarray): The eigenvalues of A.
    """
    raise NotImplementedError("Problem 6 Incomplete")
