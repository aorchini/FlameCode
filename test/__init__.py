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
    # BC1: I have the attachment BC at r = 1, always
    # BC2: I need to set dF/dr = 0 at r = 0 iff Mark != 0
    [qMean, r, FMean] = steady_flame_area_FD3(varList['Mark'], varList['beta'], varList['axi'], varList['Nr'])
    r = r * varList['beta']

    # Calculate mean flame derivatives
    dFMeanDr =  derivsnew.FD1_CT2_D(FMean, r[1] - r[0])
    d2FMeanDr2 = derivsnew.FD2_CT2_D(FMean, r[1] - r[0])

    #Apply BC smooth tip:
    if(varList['Mark']!=0.0):
        dFMeanDr[-1] = 0.0

    # Use correct number of points. Remember that the extrems need to be set depending on the BC!
    # The attach BC (first point) is always assumed to be true and removed from the vector list
    if(varList['Mark']==0):
        Nr = varList['Nr'] / 2
        dFMeanDr = dFMeanDr[1:]
        d2FMeanDr2 = d2FMeanDr2[1:]
        r = r[1:]
    # The smooth BC holds only if Mark!=0 (second derivatives appear): remove also the last point
    else:
        Nr = varList['Nr'] / 2 - 1
        dFMeanDr = dFMeanDr[1:-1]
        d2FMeanDr2 = d2FMeanDr2[1:-1]
        r = r[1:-1]

    # Calculate geometric values
    den = 1 + varList['beta'] * varList['beta'] * dFMeanDr * dFMeanDr
    dR = r[1] - r[0]
    # Set Nx equal to Nr for now.
    # The implementation is more complicated if they differ, and need to interpolate between values.
    Nx = Nr

    # Nonuniform grid spacing along x!
    # Nx = length(dx) has to hold.
    dx = np.empty(len(FMean) - 1)
    for ii in range(1, len(FMean)):
        dx[ii - 1] = FMean[ii] - FMean[ii - 1]

    [A, B, C, tau] = loadAcoustics(varList['xf'], varList['Tf'], varList['acModes'], varList['beta'])

    Matrix = buildMatrix(Nr, dR, varList['beta'], den, r, FMean, dFMeanDr, d2FMeanDr2, varList['Mark'], varList['acModes'], A,
                         B, C, Nx, dx, tau, qMean, varList['convSpeed'])

    [d, W, V] = eigProblem.solveEigProb(Matrix)
    [dnew, Wnew, Vnew] = eigProblem.selectUnstable(d, W, V)

    print dnew / (2. * np.pi)


if __name__ == "__main__":
    main()