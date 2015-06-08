__author__ = 'Alessandro Orchini'

# !/usr/bin/python
"""Functions that evaluate finite different derivatives or differentiation matrices

names: mehtodN_stencilM_return
    - method:
        -# FD = FiniteDifference
        -# NU = NonUniform (Finite Difference)
        -# CB = Chebyshev
    - N = order of derivative
    - stencil:
        -# CT = Central
        -# FW = Forward
        -# BW = Backward
        -# CB = Chebyshev
    - M = order of accuracy
    - return:
        -# D = Derivative
        -# DM = Differentiation Matrix
"""

import numpy as np


def FD1_CT2_D(f, dx):
    """
    Returns f first derivative, second order dx accurate, finite difference centered
    :param f: function
    :param dx: grid space
    :return: df/dx
    """

    dfdr = np.zeros([len(f), len(f)])

    # First point use forward scheme first order
    dfdr[0, 0] = -1.5
    dfdr[0, 1] = +2.0
    dfdr[0, 2] = -0.5

    for ii in range(1, len(f) - 1):
        dfdr[ii][ii - 1] = -0.5
        dfdr[ii][ii + 1] = +0.5

    # Last point use backward scheme first order
    dfdr[-1, -1] = +1.5
    dfdr[-1, -2] = -2.0
    dfdr[-1, -3] = +0.5

    return np.dot(dfdr / dx, f)


def FD1_FW2_D(f, dx):
    """
    Returns f first derivative, second order dx accurate, finite difference forward
    :param f: function
    :param dx: grid space
    :return: df/dx
    """

    dfdr = np.zeros([len(f), len(f)])

    for ii in range(0, len(f) - 2):
        dfdr[ii][ii] = -1.5
        dfdr[ii][ii + 1] = +2.0
        dfdr[ii][ii + 2] = -0.5

    # Vorlast point use centered scheme to keep second order
    dfdr[-2, -3] = -0.5
    dfdr[-2, -1] = +0.5

    # Last point use backward scheme first order
    dfdr[-1, -1] = +1.5
    dfdr[-1, -2] = -2.0
    dfdr[-1, -3] = +0.5

    return np.dot(dfdr / dx, f)


def FD1_BW2_D(f, dx):
    """
    Returns f first derivative, second order dx accurate, finite difference backward
    :param f: function
    :param dx: grid space
    :return: df/dx
    """

    dfdr = np.zeros([len(f), len(f)])

    # First point use forward scheme first order
    dfdr[0, 0] = -1.5
    dfdr[0, 1] = +2.0
    dfdr[0, 2] = -0.5

    # Second point use centered scheme to keep second order
    dfdr[1, 0] = -0.5
    dfdr[1, 2] = +0.5

    for ii in range(2, len(f)):
        dfdr[ii][ii] = +1.5
        dfdr[ii][ii - 1] = -2.0
        dfdr[ii][ii - 2] = +0.5

    return np.dot(dfdr / dx, f)


def FD2_CT2_D(f, dx):
    """
    Returns f second derivative, second order dx accurate, finite difference centered
    :param f: function
    :param dx: grid space
    :return: df/dx
    """

    dfdr = np.zeros([len(f), len(f)])

    # First point use forward scheme first order
    dfdr[0, 0] = +1.0
    dfdr[0, 1] = -2.0
    dfdr[0, 2] = +1.0

    for ii in range(1, len(f) - 1):
        dfdr[ii][ii - 1] = +1.0
        dfdr[ii][ii] = -2.0
        dfdr[ii][ii + 1] = +1.0

    # Last point use backward scheme first order
    dfdr[-1, -1] = +1.0
    dfdr[-1, -2] = -2.0
    dfdr[-1, -3] = +1.0

    return np.dot(dfdr / (dx * dx), f)


def NU1_BW2_D(f, dx):
    """
    Returns f first derivative, second order dx accurate, finite difference backward, nonuniform mesh
    :param f: function
    :param dx: vector of grid spacings
    :return: df/dx
    """
    if len(dx)+1 != len(f):
        print 'Error: dimensions of N and dx mismatch in nonuniform differentiation matrix calculation'
        return

    dfdr = np.zeros([len(f), len(f)])

    # First point, use forward scheme
    dfdr[0,0] = -(2.0*dx[0] + 1.0*dx[1]) / (dx[0] * (dx[0] + dx[1]))
    dfdr[0,1] = +(dx[0] + dx[1]) / (dx[0] * dx[1])
    dfdr[0,2] = -(dx[0]) / (dx[1] * (dx[0] + dx[1]))

    # Second point, use centered scheme
    dfdr[1,0] = -(dx[1]) / (dx[0] * (dx[1] + dx[0]))
    dfdr[1,1] = (dx[1] - dx[0]) / (dx[0] * dx[1])
    dfdr[1,2] = +(dx[0]) / (dx[1] * (dx[0] + dx[1]))

    for ii in range(2, len(f)):
        dfdr[ii, ii] = (2. * dx[ii - 1] + dx[ii - 2]) / (dx[ii - 1] * (dx[ii - 1] + dx[ii - 2]))
        dfdr[ii, ii - 1] = -(dx[ii - 2] + dx[ii - 1]) / (dx[ii - 2] * dx[ii - 1])
        dfdr[ii, ii - 2] = +(dx[ii - 1]) / ((dx[ii - 2] + dx[ii - 1]) * dx[ii - 2])

    return np.dot(dfdr,f)



def FD1_CT2_DM(N, dx):
    """
    Returns first differentiation matrix, second order dx accurate, finite difference centered
    :param N: matrix dimension
    :param dx: grid space
    :return: d/dx
    """

    dfdr = np.zeros([N, N])

    # First point use forward scheme
    dfdr[0, 0] = -1.5
    dfdr[0, 1] = +2.0
    dfdr[0, 2] = -0.5

    for ii in range(1, N - 1):
        dfdr[ii][ii - 1] = -0.5
        dfdr[ii][ii + 1] = +0.5

    # Last point use backward scheme
    dfdr[-1, -1] = +1.5
    dfdr[-1, -2] = -2.0
    dfdr[-1, -3] = +0.5

    return dfdr / dx


def FD1_FW2_DM(N, dx):
    """
    Returns first differentiation matrix, second order dx accurate, finite difference forward
    :param N: matrix dimension
    :param dx: grid space
    :return: d/dx
    """

    dfdr = np.zeros([N, N])

    for ii in range(0, N - 2):
        dfdr[ii][ii] = -1.5
        dfdr[ii][ii + 1] = +2.0
        dfdr[ii][ii + 2] = -0.5

    # Vorlast point use centered scheme
    dfdr[-2, -3] = -0.5
    dfdr[-2, -1] = +0.5

    # Last point use backward scheme
    dfdr[-1, -1] = +1.5
    dfdr[-1, -2] = -2.0
    dfdr[-1, -3] = +0.5

    return dfdr / dx


def FD1_BW2_DM(N, dx):
    """
    Returns first differentiation matrix, second order dx accurate, finite difference backward
    :param N: matrix dimension
    :param dx: grid space
    :return: d/dx
    """

    dfdr = np.zeros([N, N])

    # First point use forward scheme
    dfdr[0, 0] = -1.5
    dfdr[0, 1] = +2.0
    dfdr[0, 2] = -0.5

    # Second point use centered scheme
    dfdr[1, 0] = -0.5
    dfdr[1, 2] = +0.5

    for ii in range(2, N):
        dfdr[ii][ii] = +1.5
        dfdr[ii][ii - 1] = -2.0
        dfdr[ii][ii - 2] = +0.5

    return dfdr / dx


def FD2_CT2_DM(N, dx):
    """
    Returns second differentiation matrix, second order dx accurate, finite difference centered
    :param N: matrix dimension
    :param dx: grid space
    :return: d/dx
    """

    dfdr = np.zeros([N, N])

    # First point use forward scheme
    dfdr[0, 0] = +1.0
    dfdr[0, 1] = -2.0
    dfdr[0, 2] = +1.0

    for ii in range(1, N - 1):
        dfdr[ii][ii - 1] = +1.0
        dfdr[ii][ii] = -2.0
        dfdr[ii][ii + 1] = +1.0

    # Last point use backward scheme
    dfdr[-1, -1] = +1.0
    dfdr[-1, -2] = -2.0
    dfdr[-1, -3] = +1.0

    return dfdr / (dx * dx)


def NU1_BW2_DM(N, dx):
    """
    Returns first differentiation matrix, second order dx accurate, finite difference backward, nonuniform mesh
    :param N: matrix dimension
    :param dx: vector of grid spacings
    :return: d/dx
    """
    if len(dx)+1 != N:
        print 'Error: dimensions of N and dx mismatch in nonuniform differentiation matrix calculation'
        return

    dfdr = np.zeros([N, N])

    # First point, use forward scheme
    dfdr[0,0] = -(2.0*dx[0] + 1.0*dx[1]) / (dx[0] * (dx[0] + dx[1]))
    dfdr[0,1] = +(dx[0] + dx[1]) / (dx[0] * dx[1])
    dfdr[0,2] = -(dx[0]) / (dx[1] * (dx[0] + dx[1]))

    # Second point, use centered scheme
    dfdr[1,0] = -(dx[1]) / (dx[0] * (dx[1] + dx[0]))
    dfdr[1,1] = (dx[1] - dx[0]) / (dx[0] * dx[1])
    dfdr[1,2] = +(dx[0]) / (dx[1] * (dx[0] + dx[1]))

    for ii in range(2, N):
        dfdr[ii, ii] = (2. * dx[ii-1] + dx[ii - 2]) / (dx[ii-1] * (dx[ii-1] + dx[ii - 2]))
        dfdr[ii, ii - 1] = -(dx[ii - 2] + dx[ii - 1]) / (dx[ii - 2] * dx[ii - 1])
        dfdr[ii, ii - 2] = (dx[ii - 1]) / ((dx[ii - 2] + dx[ii - 1]) * dx[ii - 2])

    return dfdr