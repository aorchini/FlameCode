__author__ = 'Alessandro Orchini'

# ! /usr/bin/python
"""Function that import acoustics state-space matrices from Matlab files."""

import numpy as np
import scipy.io


def loadAcoustics(xf, Tf, M, beta):
    """For now, import acoustic matrices that have been calculated with an acoustic solver.
    E.g., the diy from me and Simon. Consider write an acoustic solver, too"""

    filename = ('/home/ao352/DIDF/kubiraj_acoustics/Lot_ac_xf_%.8f_Tf_%d.00000000/ss_xf_%.8f_Tf_300.mat') % (xf, Tf, xf)
    ACss = scipy.io.loadmat(filename)

    # NB: mat is a dictionary here!! To see a list of the keys use
    # for key, value in ACss.iteritems() :
    # print key
    #
    # The list of keys is: A,B,C,D,mf,geom,x,normp,normu.
    # mf and geom are also dictionaries!

    index = np.where(ACss['x'][:, 0] > 1e-5)[0]

    if (np.size(index)):
        indxf = index[0] - 1
    else:
        indxf = len(ACss['x'])

    A = ACss['A'][0:M, 0:M]
    B = ACss['B'][0:M]
    C = ACss['C'][indxf + len(ACss['x']), 0:M]
    #D = ACss.D(indxf+length(ACss.x)) %This should always be zero.

    # Scaling between acoustic and convection times
    Ltube = ACss['x'][-1, 0] - ACss['x'][0, 0]
    LM_Lf = ACss['mf']['M'][0, 0][0, 0] * Ltube / (0.005 * beta)  #5mm = radius of tube in the experiment
    tau = 1. / LM_Lf

    return [A, B, C, tau]