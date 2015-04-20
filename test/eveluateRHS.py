__author__ = 'Alessandro Orchini'

import numpy as np
from derivs import DfDr

def evaluateRHS(r,beta,K,St,t,f):
    """Evaluate rhs of df/dt = rhs
    """
    rhs = beta*beta/(1+beta*beta)*DfDr(f,r) + np.cos(St*(t-K*(1-r)))
    rhs[-1] = 0 #BC

    return rhs


