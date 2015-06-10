#!/usr/bin/python
"""Functions that construct the linearized thermoacoustic system.
The acoustics is not solved here, it is imported from Matlab files generated with my code or LOTAN.
The only requirement is that it is in state space form.
Using methods I can generalize to include for more
- flame models (currently only G-equation)
- flow models (currently convective, incompressible model)
"""

__author__ = 'Alessandro Orchini'


import numpy as np

import derivsnew

def buildMatrix(Nr, dR, beta, den, r, FMean, dFMeanDr, d2FMeanDr2, Mark, M, A, B, C, Nx, dx, tau, heatRelease, K):

    m_ff = M_ff(Nr, dR, beta, den, r, dFMeanDr, d2FMeanDr2, Mark)
    m_fs = M_fs(M, Nr, C, Nx, r, dFMeanDr, dx, beta)
    m_fv = M_fv(Nx, Nr, r, dFMeanDr, dx, beta, Mark)
    m_sf = M_sf(beta, Nr, M, tau, B, r, den, Mark, dFMeanDr, d2FMeanDr2, FMean, heatRelease)
    m_ss = M_ss(M, A, tau)
    m_sv = M_sv(Nx, M)
    m_vf = M_vf(Nx, Nr)
    m_vs = M_vs(M, C, Nx, K, dx)
    m_vv = M_vv(K, Nx, dx, Mark)

    M_f = np.concatenate((m_ff, m_fs, m_fv), axis=1)
    M_s = np.concatenate((m_sf, m_ss, m_sv), axis=1)
    M_v = np.concatenate((m_vf, m_vs, m_vv), axis=1)

    Matrix = np.concatenate((M_f, M_s, M_v), axis=0)

    # Check dimensions
    checkM = np.shape(Matrix)
    if (checkM[0] != Nx + Nr + M or checkM[1] != Nx + Nr + M):
        print 'Dimensions of A are not consistent with M!'
        exit()

    return Matrix


def M_ff(Nr, dR, beta, den, r, dFMeanDr, d2FMeanDr2, Mark):
    """Build M_ff matrix
    Evaluate df1/dr and d2f1/dr2 and build the matrix
    """
    M_ff = np.zeros([Nr, Nr])

    if(Mark!=0):
        D1 = derivsnew.FD1_BW2_DM(Nr+2, dR)
        D2 = derivsnew.FD2_CT2_DM(Nr+2, dR)
    else:
        D1 = derivsnew.FD1_BW2_DM(Nr+1, dR)
        D2 = derivsnew.FD2_CT2_DM(Nr+1, dR)


    D1 = D1[1:,1:] # Apply BC attachment, remove first column and row (ok because the first point is fixed to 0)
    D2 = D2[1:,1:] # Apply BC attachment, remove first column and row (ok because the first point is fixed to 0)
    if(Mark!=0):
        D1 = D1[0:-1,0:-1] # Apply symmetry BC: remove last column and row (ok because it is a backward scheme)
        D2 = D2[0:-1,0:-1] # Apply symmetry BC: remove last column and row (
        # here I need to correct the vorlast point (new last) because I have removed the last column and D2 uses a centered scheme
        # from the smooth tip condition on the backward scheme I have that f(end) = 2.0/1.5 f(end-1) - 0.5/1.5 f(end-2)
        D2[-1,-1] = -2./3. / (dR * dR) # = -2.0 + 2.0/1.5
        D2[-1,-2] = +2./3. / (dR * dR) # = +1.0 - 0.5/1.5

    for ii in range(0, Nr):
        for jj in range(0, Nr):
            if(r[jj]==0):
                M_ff[jj, ii] = D1[jj, ii] * dFMeanDr[jj] / np.sqrt(den[jj]) + \
                            -Mark * (D2[jj, ii] / den[jj] +
                                        -2 * beta * beta * D1[jj, ii] * dFMeanDr[jj] * d2FMeanDr2[jj] / (den[jj] * den[jj]))
            else:
                M_ff[jj, ii] = D1[jj, ii] * dFMeanDr[jj] / np.sqrt(den[jj]) + \
                            -Mark * (D2[jj, ii] / den[jj] +
                                        -2 * beta * beta * D1[jj, ii] * dFMeanDr[jj] * d2FMeanDr2[jj] / (den[jj] * den[jj]) +
                                       +D1[jj, ii] / r[jj])


    M_ff = -beta * beta / np.sqrt(1 + beta * beta) * M_ff

    # Checked, it is ok!
    return M_ff


def M_fs(M, Nr, C, Nx, r, dFMeanDr, dx, beta):
    """BC Correction for v:
    the term reads coeff*1/2*dFMeandr*v(0)  with v(0) = sum C_i s_i
    the first coeff is DN(1,0) = -1/dx(1), the second one is D(2,0) = (dx(2))/((dx(1)+dx(2))*dx(1));
    This works only if Nr = Nx    M_fs = np.zeros([Nr, M])
    """

    M_fs = np.zeros([Nr, M])

    for ii in range(0, M):
        # NB: THE DIVISION BY BETA IS WRONG!! THERE WAS A BUG IN LSGEND2D. IF I WANT TO REPRODUCE
        # THOSE RESULTS I NEED TO REPRODUCE THE BUG (DONE DIVIDING BY BETA!)
        M_fs[0, ii] = -(dx[1]) / (dx[0] * (dx[1] + dx[0])) * (+1. / 2.) * r[0] * dFMeanDr[0] * C[ii] / beta
        M_fs[1, ii] = +(dx[1]) / ((dx[0] + dx[1]) * dx[0]) * (+1. / 2.) * r[1] * dFMeanDr[1] * C[ii] / beta

    # Checked, it is ok!
    return M_fs


def M_fv(Nx, Nr, r, dFMeanDr, dx, beta, Mark):
    """This assumes that Nr = Nx, and that we have f(1) = f(r=R) and v(1) = v(x=0);
    """

    M_fv = np.zeros([Nr, Nx])

    if(Mark!=0):
        DN = derivsnew.NU1_BW2_DM(Nx+2, dx)
    else:
        DN = derivsnew.NU1_BW2_DM(Nx+1, dx)

    DN = DN[1:,1:] # Apply BC attachment, remove first column and row (a correction is needed as is added into M_fs)
    if(Mark!=0):
        DN = DN[0:-1,0:-1] # Apply symmetry BC: remove last column and row (ok because it is a backward scheme)

    for ii in range(0, Nx):
        for jj in range(0, Nr):
            # THE DIVISION BY BETA IS WRONG!! THERE WAS A BUG IN LSGEND2D. IF I WANT TO REPRODUCE
            # THOSE RESULTS I NEED TO REPRODUCE THE BUG (DONE DIVIDING BY BETA!)
            M_fv[jj, ii] = 1. / 2. * r[jj] * DN[jj, ii] * dFMeanDr[jj] / beta

    print 'Note: I am voluntarly imposing a bug in M_fv M_fs with beta to reproduce the results as in Lieuwen paper.'

    for ii in range(0, Nx):
        M_fv[ii, ii] = M_fv[ii, ii] + 1.

    # Checked, it is ok!
    return M_fv


def M_sf(beta, Nr, M, tau, B, r, den, Mark, dFMeanDr, d2FMeanDr2, yMean, heatRelease):
    # Be very careful with this, is the most complicated bit
    # Build finite difference coefficients

    dR = r[1] - r[0]

    if(Mark!=0):
        D1 = derivsnew.FD1_BW2_DM(Nr+2, dR)
        D2 = derivsnew.FD2_CT2_DM(Nr+2, dR)
    else:
        D1 = derivsnew.FD1_BW2_DM(Nr+1, dR)
        D2 = derivsnew.FD2_CT2_DM(Nr+1, dR)

    D1 = D1[1:,1:] # Apply BC attachment, remove first column and row (ok because the first point is fixed to 0)
    D2 = D2[1:,1:] # Apply BC attachment, remove first column and row (ok because the first point is fixed to 0)
    if(Mark!=0):
        D1 = D1[0:-1,0:-1] # Apply symmetry BC: remove last column and row (ok because it is a backward scheme)
        D2 = D2[0:-1,0:-1] # Apply symmetry BC: remove last column and row (
        # here I need to correct the vorlast point (new last) because I have removed the last column and D2 uses a centered scheme
        # from the smooth tip condition on the backward scheme I have that f(end) = 2.0/1.5 f(end-1) - 0.5/1.5 f(end-2)
        D2[-1,-1] = -2./3. / (dR * dR) # = -2.0 + 2.0/1.5
        D2[-1,-2] = +2./3. / (dR * dR) # = +1.0 - 0.5/1.5

    # Integration coefficients
    mu = np.ones(Nr)
    # First point (r = 1) separately later
    mu[0] = 7. / 6.
    mu[1] = 23. / 24.

    # If Mark !=0 the last point is fixed by the BC and removed from the vector list.
    if(Mark!=0):
        mu[-1] = 7. / 6.
        mu[-2] = 23. / 24.
    else:
        mu[-1] = 3. / 8. # This has weight zero if the flame is axisymmetric. This is seen by dq already because r=0 there
        mu[-2] = 7. / 6.
        mu[-3] = 23. / 24.

    mu = mu * abs(dR)


    # Effect of the flame on the acoustics
    # Given by tau*B*int_0^1 (dq')
    # with int_0^1 (dq')= int_0^1 [(1-Mark*k)*sqrt( (\beta*DF/dr)^2 + (dG/dx)^2)*r\beta^2]'*dr
    dq = np.zeros(Nr)

    for jj in range(0, Nr):
        for ii in range(0, Nr):
            if(r[ii]==0):
                dq[jj] += + mu[ii] * r[ii] * (dFMeanDr[ii] * D1[ii, jj] / np.sqrt(den[ii]) +
                                              - Mark * (
                                                  D2[ii, jj] / den[ii] - 2. * beta * beta * d2FMeanDr2[ii] * dFMeanDr[ii] *
                                                  D1[ii, jj] / (den[ii] * den[ii]) ) )
            else:
                dq[jj] += + mu[ii] * r[ii] * (dFMeanDr[ii] * D1[ii, jj] / np.sqrt(den[ii]) +
                                              - Mark * (
                                                  D2[ii, jj] / den[ii] - 2. * beta * beta * d2FMeanDr2[ii] * dFMeanDr[ii] *
                                                  D1[ii, jj] / (den[ii] * den[ii]) + D1[ii, jj] / r[ii] ) )


    # Add contribution of first point
    corrCoeff = 3. / 8. * abs(dR)
    dFdRCorr = (2. * yMean[1] - 1. / 2. * yMean[2]) / dR  # because yMean(1) = 0
    dF2dR2Corr = (-5. * yMean[1] + 4. * yMean[2] - yMean[3]) / (dR * dR)  # because yMean(1) = 0
    denCorr = 1. + beta * beta * dFdRCorr * dFdRCorr

    # Forward first derivative
    D1C1 = (2.0) / dR
    D1C2 = (-1.0 / 2.0) / dR

    # Forward second derivative
    D2C1 = (-5.0) / (dR * dR)
    D2C2 = (+4.0) / (dR * dR)
    D2C3 = (-1.0) / (dR * dR)

    dq[0] += corrCoeff * 1.0 * (dFdRCorr * D1C1 / np.sqrt(denCorr) +
                                -Mark * (D2C1 / denCorr - 2. * beta * beta * dF2dR2Corr * dFdRCorr * D1C1 / (
                                    denCorr * denCorr) + D1C1 / 1.0))

    dq[1] += corrCoeff * 1.0 * (dFdRCorr * D1C2 / np.sqrt(denCorr) +
                                -Mark * (D2C2 / denCorr - 2. * beta * beta * dF2dR2Corr * dFdRCorr * D1C2 / (
                                    denCorr * denCorr) + D1C2 / 1.0))

    dq[2] += corrCoeff * 1.0 * (0.0 - Mark * ( D2C3 / denCorr - 0.0 + 0.0 ) )

    # Need to account also for the last point, but only in the case in which the flame is not axysimmetric
    # If axisymmetirc last point (r = 0) weights zero on the integral.
    M_sf = tau * B * dq
    M_sf = M_sf / heatRelease  # nondimensioanal q'/qMean

    # Checked, it is ok!
    return M_sf


def M_ss(M, A, tau):
    checkS = np.shape(A)
    if (checkS[0] != M or checkS[1] != M):
        print 'Dimensions of A are not consistent with M!'
        exit()

    M_ss = A * tau

    # Checked, it is ok!
    return M_ss


def M_sv(Nx, M):
    """No coupling
    """
    M_sv = np.zeros([M, Nx])

    # Checked it is ok!
    return M_sv


def M_vf(Nx, Nr):
    """ No coupling
    """

    M_vf = np.zeros([Nx, Nr])

    # Checked it is ok!
    return M_vf


def M_vs(M, C, Nx, K, dx):
    """BC Correction for v:
    the term reads coeff*1/2*dFMeandr*v(0)  with v(0) = sum C_i s_i
    the first coeff is DN(1,0) = -1/dx(1), the second one is D(2,0) = (dx(2))/((dx(1)+dx(2))*dx(1))"""

    M_vs = np.zeros([Nx, M])

    # This works only if Nr = Nx
    for ii in range(0, M):
        M_vs[0, ii] = -(dx[1]) / (dx[0] * (dx[1] + dx[0])) * (-1. / K) * C[ii]
        M_vs[1, ii] = +(dx[1]) / ((dx[0] + dx[1]) * dx[0]) * (-1. / K) * C[ii]


    # Checked, it is ok!
    return M_vs


def M_vv(K, Nx, dx, Mark):

    if(Mark!=0):
        DN = derivsnew.NU1_BW2_DM(Nx+2, dx)
    else:
        DN = derivsnew.NU1_BW2_DM(Nx+1, dx)

    DN = DN[1:,1:] # Apply BC attachment, remove first column and row (a correction is needed as is added into M_vs)
    if(Mark!=0):
        DN = DN[0:-1,0:-1] # Apply symmetry BC: remove last column and row (ok because it is a backward scheme)

    M_vv = -1. / K * DN

    # Checked, it is ok
    return M_vv

