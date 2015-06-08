#!/usr/bin/python
"""Time Integrate the linear equations. One can use the time series to evaluate a TF.
"""
__author__ = 'Alessandro Orchini'

import numpy as np

from RungeKutta import RK
from sinefit import sinefit


def timeIntegration(beta, convSpeed, St, NyqRate):
    """Time marching procedure"""


    ## General parameters and plotting
    ## plot(r,1-r+0.1*f,'r-','Linewidth',2)
    ## ylim([0 1.5])
    # plot(r,f,'r-','Linewidth',2)
    # ylim([-0.2 0.2])
    # xlim([0 1])
    # drawnow

    St2 = St * (1. + beta * beta) / (beta * beta)  # Reduced frequency, Preetham notation
    eta = convSpeed * beta * beta / (1. + beta * beta)  # Reduced speed. Preetham notation

    ## Space settings

    lambda1 = 2. * np.pi / St / convSpeed  # Forcing wavelength
    lambda2 = 2. * np.pi * beta * beta / St / (1 + beta * beta)  # Natural wavelength
    lambda3 = NyqRate / 10.
    lambd = min(lambda1, lambda2, lambda3)
    r = np.linspace(0., 1., NyqRate / lambd)
    f = np.zeros(len(r))
    dr = r[1] - r[0]

    ## Time settings
    CFL = 0.6  # CFL number
    nc = 2.  # number of cycles of the steady solution
    dt1 = CFL * dr  # grid dependent time-step
    dt2 = 2. * np.pi / St / NyqRate  # Nyquist condition for q
    dt = min(dt1, dt2)
    t = 0.
    tf = (beta * beta + 1.) / (beta * beta) + 2. * np.pi / St * nc  # = time for the transient to disappear + nc cycles

    ## integration settings
    nucoffs = np.ones(len(r))
    nucoffs[0] = 3. / 8.
    nucoffs[-1] = 3. / 8.
    nucoffs[1] = 7. / 6.
    nucoffs[-2] = 7. / 6.
    nucoffs[2] = 23. / 24.
    nucoffs[-3] = 23. / 24.
    nucoffs = nucoffs * dr

    sf = []
    tstore = []
    q = []
    u = []

    while t < tf:
        f = RK(r, beta, convSpeed, St, t, dt, f)
        t = t + dt

        if (t > (beta * beta + 1.) / (beta * beta)):  # store only steady solution, neglect transient
            sf.append(f)
            tstore.append(t)
            q.append(2. * beta * beta / (1 + beta * beta) * np.dot(nucoffs, f))
            u.append(np.cos(St * t))

            # perc = t/tf*100.;
            # if(perc>=targ)
            # fprintf(['Complete percentage = ', num2str(perc,'%.0f'),'%%\n']);
            #   targ = targ+10;
            # end

            ## plotting
            ## plot(r,1-r+0.1*f,'r-','Linewidth',2)
            ## ylim([0 1.5])
            # plot(r,f,'r-','Linewidth',2)
            # ylim([-0.2 0.2])
            # xlim([0 1])
            # drawnow

    sf = np.asarray(sf)
    tstore = np.asarray(tstore)
    q = np.asarray(q)
    u = np.asarray(u)

    [paramsq, yestq, yresq, rmserrq] = sinefit(q, tstore, St / (2 * np.pi), 0, 0, 0)
    [paramsu, yestu, yresu, rmserru] = sinefit(u, tstore, St / (2 * np.pi), 0, 0, 0)

    gain = paramsq[1] / paramsu[1];
    phase = paramsq[3] - paramsu[3];

    ModATF = 2 * (eta - 1 + np.exp(1j * eta * St2) - eta * np.exp(1j * St2)) / ((1 - eta) * eta * St2 * St2)

    print('St2 = ' + str(St2) + '. Error gain = ' + str(abs(gain - abs(ModATF)) / abs(ModATF) * 100))

    print gain
    print phase


timeIntegration(2., 0.7, 10., 50)