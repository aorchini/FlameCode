__author__ = 'Alessandro Orchini'
#!/usr/bin/python

"""Mainfile
\brief: Linear ThermoAcoustic Solver by Alessandro Orchini \n
Started on 22-03-2015 \n
Last Updated on 21-04-2015

Thermoacoustic solver. The program has a few functionalities:
    - Acoustics
        -# Given a geometry, calculate the acoustic eigenvalues (Helmholtz Solver)
        -# Given a geometry, calculate the frequency response of a compact flame (acoustic TF)
        -# Fit a frequency response onto a state-space model
        -# More?
    - Flame
        -# Given a linear flame model, performs time marching simulations of forced input response
        -# Given a linear flame model in the time domain, calculate its Transfer Function
    - Thermoacoustics
        -# Couples flame and acoustics responses into an eigenvalue problem and solves it
        -# Uses the harmonic balance criterion to assess the stability of the system
"""

import numpy as np

import derivs
import derivsnew
import eigProblem
import pylab as pl

from steadyFlame import steady_flame_area_FD3
from loadAcoustics import loadAcoustics
from subMatrices import buildMatrix


def main():
    """Main Function.

    Define variables and call appropriate routines
    """

    varList = {'beta': 6., 'convSpeed': 1.2, 'Mark': 0.02, 'axi': 1, 'acModes': 4, 'Nr': 801, 'Tf': 600., 'xf': 0.15}

    # Solve steady flame.
    [qMean, xMean, yMean] = steady_flame_area_FD3(varList['Mark'], varList['beta'], varList['axi'], varList['Nr'])
    xMean = xMean * varList['beta']

    # Use all points. Remember that the extrems need to be set depending on the BC!
    # I always have the attachment BC at r = 1
    # Plus I need to set dF/dr = 0 at r = 0 if Mark != 0
    r = xMean[:]
    FMean = yMean[:]

    # Calculate mean flame derivatives
    dFMeanDr =  derivsnew.FD1_BW2_D(FMean, r[1] - r[0])
    d2FMeanDr2 = derivsnew.FD2_CT2_D(FMean, r[1] - r[0])

    #Check that BC are satisfied:
    print 'Percentage difference on attachment BC:'
    print abs((FMean[0])/(1e-12))*100

    if(varList['Mark']!=0.0):
        print 'Percentage difference on smooth tip:'
        print abs( (FMean[-1] - (2.0*FMean[-2] - 0.5*FMean[-3])/1.5 ) / ((2.0*FMean[-2] - 0.5*FMean[-3])/1.5))*100


    den = 1 + varList['beta'] * varList['beta'] * dFMeanDr * dFMeanDr

    Nr = varList['Nr'] / 2 - 1
    dR = r[1] - r[0]
    # Set equal to Nr for now.
    # The implementation is more complicated if they differ, and need to interpolate between values.
    Nx = Nr

    # Nonuniform grid spacing along x!
    # I do not need the last point, but I need the first (I use a backward scheme).
    # Nx = length(dx) has to hold.
    dx = np.empty(len(yMean) - 1)
    for ii in range(1, len(yMean)):
        dx[ii - 1] = yMean[ii] - yMean[ii - 1]

    [A, B, C, tau] = loadAcoustics(varList['xf'], varList['Tf'], varList['acModes'], varList['beta'])

    Matrix = buildMatrix(Nr, dR, varList['beta'], den, r, dFMeanDr, d2FMeanDr2, varList['Mark'], varList['acModes'], A,
                         B, C, Nx, dx, tau, qMean, varList['convSpeed'], yMean)

    [d, W, V] = eigProblem.solveEigProb(Matrix)
    [dnew, Wnew, Vnew] = eigProblem.selectUnstable(d, W, V)

    print dnew / (2. * np.pi)


if __name__ == "__main__":
    main()