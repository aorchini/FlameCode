__author__ = 'Alessandro Orchini'

#!/usr/bin/python
"""Functions that evaluate finite different derivatives or differentiation matrices
"""

import numpy as np

def DfDr(f, r):
    """\brief second order forward finite difference first derivative

    Evaluate the first derivative of f using a second order finite difference forward scheme
    """

    dr = r[1] - r[0]

    DN = np.zeros([len(r), len(r)])

    for ii in range(0, len(r) - 2):
        DN[ii, ii] = -3. / 2.
        DN[ii, ii + 1] = 2.
        DN[ii, ii + 2] = -1. / 2.

    # Vorlast row is order 1
    DN[-2, -3] = -1. / 2.
    DN[-2, -1] = +1. / 2.

    # Last row I use backward FD
    DN[-1, -1] = +3. / 2.
    DN[-1, -2] = -2.
    DN[-1, -3] = +1. / 2.

    DN = DN / dr
    dfdr = np.dot(DN, f)

    return dfdr


def diff_FD(ys, N, dx):
    """
    First derivative - Second order finite difference
    :param ys: current shape
    :param N: number of points
    :param dx: grid size
    :return: first derivative (second order finite difference)
    """
    n = (N + 1) / 2
    Dys = np.zeros(len(ys))
    Dys[0] = (-1.5 * ys[0] + 2. * ys[1] - 0.5 * ys[2]) / dx
    for i in range(1, n - 1):
        Dys[i] = (-0.5 * ys[i - 1] + 0.5 * ys[i + 1]) / dx

    Dys[n - 1] = 0.
    return Dys


def diff2_FD(ys, N, dx):
    """
    Second derivative - Second order finite difference
    :param ys: current shape
    :param N: number of points
    :param dx: grid size
    :return: second derivative (second order finite difference)
    """
    n = (N + 1) / 2
    D2ys = np.zeros(len(ys))
    D2ys[0] = (1.0 * ys[0] - 2.0 * ys[1] + 1.0 * ys[2]) / (dx * dx)
    for i in range(1, n - 1):
        D2ys[i] = (1.0 * ys[i - 1] - 2.0 * ys[i] + 1.0 * ys[i + 1]) / (dx * dx)

    D2ys[n - 1] = (2.0 * ys[n - 2] - 2.0 * ys[n - 1]) / (dx * dx)
    return D2ys


def FinitDiff2Back_D1(y, dx, y0):
    """ Given a vector of length N_x and a function y defined at these values
    subject to the boundary condition y(x(1)) = 0 AND dy/dx(end) = 0 this function evaluates the
    discrete derivative of y dydx using a second order backward finite
    difference scheme, on the points x(2:end).
    The input is the function y on the points y(2), ..., y(end) (NOT y(1)
    which is zero!), and the spacing dx between grid points, which is supposed
    to be zero. If the function is not zero at y(0), the first two points are
    wrong, and need to be corrected by
    dydx(1) = dydx(1) - y(0)/dx;
    dydx(2) = dydx(2) + 1/2*y(0)/dx;
    """
    DN = np.zeros([len(y), len(y)])

    DN[0, 0] = 1.
    DN[1, 1] = 3. / 2.
    DN[1, 0] = -2.

    for ii in range(2, len(y)):
        DN[ii, ii] = 3. / 2.
        DN[ii, ii - 1] = -2.
        DN[ii, ii - 2] = 1. / 2.

    dydx = np.dot(DN, y) / dx

    dydx[0] -= y0 / dx
    dydx[1] += 1. / 2. * y0 / dx

    return dydx


def FinitDiff2Central_D2(y, dx, y0):
    """Given a vector of length N_x and a function y defined at these values
    subject to the boundary condition y(x(1)) = 0 AND dy/dx(end)=0
    this function evaluates the
    discrete derivative of y dydx using a second order central finite
    difference scheme, on the points x(2:end).
    The input is the function y on the points y(2), ..., y(end) (NOT y(1)
    which is zero!), and the spacing dx between grid points, which is supposed
    to be zero. If the function is not zero at y(0), the first two points are
    wrong, and need to be corrected by
    dydx(1) = dydx(1) + y(0)/(dx^2);
    For the last point, I use a backward scheme to avoid the use of the point
    at end which is not in the vector state.
    """
    DN = np.zeros([len(y), len(y)])

    DN[0, 0] = -2.
    DN[0, 1] = +1.

    for ii in range(1, len(y) - 1):
        DN[ii, ii - 1] = 1.
        DN[ii, ii] = -2.
        DN[ii, ii + 1] = 1.

    DN[-1, -1] = +2.
    DN[-1, -2] = -5.
    DN[-1, -3] = +4.
    DN[-1, -4] = -1.

    dydx = np.dot(DN, y) / (dx * dx)

    dydx[0] += y0 / (dx * dx)

    return dydx


def D2_FinitDiff2Central_Coeff(N, dx):
    """
    Returns the coefficient of the derivative matrix as in FinitDiff2Central
    """
    DN = np.zeros([N, N])

    DN[0, 0] = -2.
    DN[0, 1] = 1.

    for ii in range(1, N - 1):
        DN[ii, ii - 1] = 1.
        DN[ii, ii] = -2.
        DN[ii, ii + 1] = 1.

    DN[-1, -1] = 2.
    DN[-1, -2] = -5.
    DN[-1, -3] = 4.
    DN[-1, -4] = -1.

    DN = DN / (dx * dx)

    return DN


def D1NU_FinitDiff2Back_Coeff(N, dx):
    """
    Returns the coefficient of the derivative matrix as in FinitDiff2Back
    The first point is order (1), the other points order (2)
    1 is the value of the SECOND element; the first element, 0, is supposed to be zero.
    If it is not zero then need to correct the first two points: the first
    one is DN(1,0) = -1/dx(1), the second one is D(2,0) = (dx(2))/((dx(1)+dx(2))*dx(1));
    Example of usage:
    lear
    x = []
    x(1) = 0;
    while(x(end)<4)
    x(end+1) = x(end) + (rand(1)*0.05);
    end
    x(end) = 4
    y = sin(x);
    for ii=2:length(y)
    dx(ii-1) = x(ii)-x(ii-1);
    end
    plot(dx,'o')
    y = y';
    DN = D1NU_FinitDiff2Back_Coeff(length(dx),dx);
    dydx = DN*y(2:end);
    plot(x(2:end),dydx,x,cos(x),'ro')
    """

    DN = np.zeros([N, N])

    DN[0, 0] = 1. / dx[0]
    DN[1, 1] = (2. * dx[1] + dx[0]) / (dx[1] * (dx[1] + dx[0]))
    DN[1, 0] = -1.0 * (dx[0] + dx[1]) / (dx[0] * dx[1])

    for ii in range(2, N):
        DN[ii, ii] = (2. * dx[ii] + dx[ii - 1]) / (dx[ii] * (dx[ii] + dx[ii - 1]))
        DN[ii, ii - 1] = -(dx[ii - 1] + dx[ii]) / (dx[ii - 1] * dx[ii])
        DN[ii, ii - 2] = (dx[ii]) / ((dx[ii - 1] + dx[ii]) * dx[ii - 1])

    return DN


def D1_FinitDiff2Back_Coeff(N, dx):
    # Returns the coefficient of the derivative matrix as in FinitDiff2Back

    DN = np.zeros([N, N])

    DN[0, 0] = 1.
    DN[1, 1] = 3. / 2.
    DN[1, 0] = -2.

    for ii in range(2, N):
        DN[ii, ii] = 3. / 2.
        DN[ii, ii - 1] = -2.
        DN[ii, ii - 2] = 1. / 2.

    DN = DN / dx

    return DN

























