__author__ = 'Alessandro Orchini'

# ! /usr/bin/python
"""Solve a (linear) eigenvalue problem. Returns left and right eigenvectors for adjoint calculations."""

from scipy import linalg as lg
import numpy as np


def solveEigProb(M):
    """
    :param M: linearized matrix (Jacobian)
    :return: d (eigenvalues), W(left eigenvectors, for adjoint calculations), V(right eigenvectors).
    """

    d, W, V = lg.eig(M, left=True)

    return [d, W, V]


def selectUnstable(d, W, V):
    """
    :param d: eigenvalues
    :param W: left eigenvectors
    :param V: right eigenvectors
    :return: eigenvalues with non-negative growth rates and their respective eigenvectors
    """

    # Select positive growth rate
    index = np.where(np.real(d) >= 0)[0]

    dTemp = np.zeros(len(index)) + 0j
    VTemp = np.zeros((V.shape[0], len(index))) + 0j
    WTemp = np.zeros((W.shape[0], len(index))) + 0j

    for ii in range(0, len(index)):
        dTemp[ii] = d[index[ii]]
        VTemp[:, ii] = V[:, index[ii]]
        WTemp[:, ii] = W[:, index[ii]]

    d = dTemp
    V = VTemp
    W = WTemp

    # Check there are no modes with zero frequency. Print a warning otherwise
    for ii in range(0, len(d)):
        if (np.imag(d[ii]) == 0):
            print('Warning: this mode has growth rate %f (non-negative) and zero frequency. '
                  'I will ignore it.' % np.real(d[ii]))

    # Select positive frequencies
    index = np.where(np.imag(d) > 0)[0]

    dTemp = np.zeros(len(index)) + 0j
    VTemp = np.zeros((V.shape[0], len(index))) + 0j
    WTemp = np.zeros((W.shape[0], len(index))) + 0j

    for ii in range(0, len(index)):
        dTemp[ii] = d[index[ii]]
        VTemp[:, ii] = V[:, index[ii]]
        WTemp[:, ii] = W[:, index[ii]]

    d = dTemp
    V = VTemp
    W = WTemp

    return [d, W, V]