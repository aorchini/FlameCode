__author__ = 'Alessandro Orchini'

from test.eveluateRHS import evaluateRHS


def RK(r, beta, K, St, t, dt, f):
    """Third order TVD scheme
    March in time
    """

    fn = f + dt * evaluateRHS(r, beta, K, St, t, f)

    fn = 3. / 4. * f + 1. / 4. * fn + 1. / 4. * dt * evaluateRHS(r, beta, K, St, t + dt, fn)

    fn = 1. / 3. * f + 2. / 3. * fn + 2. / 3. * dt * evaluateRHS(r, beta, K, St, t + dt / 2., fn)

    return fn

