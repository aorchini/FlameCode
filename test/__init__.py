__author__ = 'Alessandro Orchini'

# ! /usr/bin/python
"""Linear ThermoAcoustic Solver by Alessandro Orchini -
Started on 22-03-2015
Last Updated on 22-03-2015
This programs runs a GUI that enables the user to evaluate a Flame Transfer Function
that then can be coupled with an acoustic network, or to directly build an eigenvalue
problem that calculates the thermoacoutic modes with the largest eigenvalues.
"""

import numpy as np

import derivs
import eigProblem

from steadyFlame import steady_flame_area_FD3
from loadAcoustics import loadAcoustics
from subMatrices import buildMatrix



def main():
    """Main Function. Define variables and call appropriate routines"""

    varList = {'beta': 4., 'convSpeed': 1.1, 'Mark': 0.0, 'axi': 1, 'acModes': 4, 'Nr': 801, 'Tf': 600., 'xf': 0.21}

    # Solve steady flame
    [qMean,xMean,yMean] = steady_flame_area_FD3(varList['Mark'], varList['beta'], varList['axi'], varList['Nr'])
    xMean = xMean * varList['beta']

    # The first and last point are set by the BC. AO: I WANT TO CHANGE THIS TO TREAT THE CASE WITHOUT CURVATURE, TOO!
    r = xMean[1:-1]
    FMean = yMean[1:-1]

    # Calculate mean flame derivatives
    dFMeanDr = 1.0*derivs.FinitDiff2Back_D1(FMean,r[1]-r[0],0) # -1 because dr is going from 1 to 0
    d2FMeanDr2 = derivs.FinitDiff2Central_D2(FMean,r[1]-r[0],0)
    den = 1+varList['beta']*varList['beta']*dFMeanDr*dFMeanDr

    # The following code is Matlab code I used for testing. It is superfluous here.
    # heatRelease = evaluateSteadyHeatRelease(Mark,r,beta,dFMeanDr,d2FMeanDr2,den,yMean); This is unnecessary.
    # disp(['Finite difference: ', num2str(heatRelease)]);
    #
    # %% Check: recompute the right hand side of the Mean G-equation: this has to be equal to 1 at each point. Compare heat release with analytical solution
    # check = den.^(0.5)/sqrt(1+beta*beta) - Mark*beta*beta/sqrt(1+beta*beta)*(d2FMeanDr2./den + dFMeanDr./r);
    # checkhR = sqrt(1+beta*beta)/(2*beta*beta);
    # disp(['Analytical solution: ', num2str(checkhR)]);
    #
    # nS = norm(check-1);
    # nhR = norm(checkhR-heatRelease);
    #
    # fprintf('\n');
    #
    # disp(['The l2 norm of the steady shape is ', num2str(nS)])
    # disp(['The l2 norm of the steady heat release (finite difference is ', num2str(nhR)])
    #
    # fprintf('\n');
    Nr = varList['Nr']/2-1
    dR = r[1]-r[0]
    # Set equal to Nr for now.
    # The implementation is more complicated if they differ, and need to interpolate between values.
    Nx = Nr

    # Nonuniform grid spacing along x!
    # I do not need the last point, but I need the first (I use a backward scheme).
    # Nx = length(dx) has to hold.
    dx = np.empty(len(yMean)-2)
    for ii in range (1,len(yMean)-1):
        dx[ii-1] = yMean[ii]-yMean[ii-1]


    [A, B, C, tau] = loadAcoustics(varList['xf'],varList['Tf'],varList['acModes'],varList['beta'])


    Matrix = buildMatrix(Nr,dR,varList['beta'],den,r,dFMeanDr,d2FMeanDr2,varList['Mark'],varList['acModes'],A,B,C,Nx,dx,tau,qMean,varList['convSpeed'],yMean)

    [d, W, V] = eigProblem.solveEigProb(Matrix)
    [dnew, Wnew, Vnew] = eigProblem.selectUnstable(d, W, V)

    print dnew/(2.*np.pi)


if __name__ == "__main__":
    main()