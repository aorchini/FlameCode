__author__ = 'Alessandro Orchini'

# ! /usr/bin/python
"""Functions that evaluate steady flame area and properties. Based on a code originally written by ICW"""

import numpy as np

import derivs


def steady_flame_area_FD3(Mark, beta, axi, N):
    """
    Function that calculates the conical mean flame shape and
    :param Mark: Markstein number = Markstein length/Flame height
    :param beta: flame aspect ratio
    :param axi: 1 = axisymmetric, 0 = not
    :param N: number of points on the flame. Should be odd (the flame is symmetrical and I solve for 1/2 of the domain)
    :return: list containing flame area, x and y coordinates
    """
    Nvec = N
    n = (N + 1) / 2
    if Mark == 0:
        if axi:
            area = np.sqrt(1.0 + np.power(beta,2.0))/(np.power(beta,2.0) * 2.0)
            xs = np.linspace(1., 0., n) / beta
            ys = beta * (1. / beta - xs)
        else:
            area = np.sqrt(1.0 + np.power(beta,2.0))/(beta)
            xs = np.linspace(1., 0., n) / beta
            ys = beta * (1. / beta - xs)

        return [area, xs, ys]

    # Procedure can use grid refinement if required - doesn't speed up much though
    tolvec = 1e-13
    N = Nvec
    n = (N + 1) / 2
    xs = np.linspace(1, 0, n) / beta
    A = 1 - Mark
    if A < 0.2:
        A = 0.2
    B = -beta * A
    ys = A + B * xs
    x1 = Mark * 2
    if x1 > 1 / beta:
        x1 = 0.9 * 1 / beta
    c = B / (2 * x1)
    a = A + 0.5 * B * x1
    inds = np.where(xs < x1)
    ys[inds] = a + c * np.power(xs[inds], 2.)

    [area, xs, ys] = steady_flame(Mark, beta, axi, N, tolvec, ys)
    return [area, xs, ys]


def steady_flame(Mark, beta, axi, N, tol, ys):
    """
    Iterative Newton scheme
    :param Mark: Markstein number = Markstein length/Flame height
    :param beta: flame aspect ratio
    :param axi: 1 = axisymmetric, 0 = not
    :param N: number of points on the flame. Should be odd (the flame is symmetrical and I solve for 1/2 of the domain)
    :param tol: tolerance for convergence
    :param ys: initial shape guess
    :return: list containing flame area, x and y coordinates
    """
    SL = 1. / np.sqrt(1. + beta * beta)
    n = (N + 1) / 2
    xs = np.linspace(1, 0, n) / beta
    dx = xs[1] - xs[0]
    # xx = np.linspace(1./beta, 0, 500)

    # loop
    i = 1
    itmax = 99
    r = grad(ys, N, dx, axi, SL, Mark, beta)
    n = [np.linalg.norm(r * (-dx))]
    #  sfigure(152);plot(xs,ys,'b.',xx,interp1(xs,ys,xx,'pchip'),-xx,interp1(xs,ys,xx,'pchip'),'b-'),ylim([0 1.2])
    #  sfigure(153);plot(r)
    while i < itmax and np.linalg.norm(r * (-dx)) > tol:
        ys = newton_step(Mark, beta, axi, dx, ys, N)
        r = grad(ys, N, dx, axi, SL, Mark, beta)
        n.append(np.linalg.norm(r * (-dx)))
        #         if (1)
        #            sfigure(152);
        #            plot(xs,ys,'b.',xx,interp1(xs,ys,xx,'pchip'),-xx,interp1(xs,ys,xx,'pchip'),'b-')
        #            ylim([0 1.2])
        #            sfigure(153);plot(r)
        #            drawnow
        #         end
        i = + 1
        #          pause
        #       disp(['The steady flame has converged at iteration ',num2str(i),' out of ',num2str(itmax),'...
        #       with a residual equal to ',num2str(n(end))])
        #       sfigure(152);
        #       plot(xs,ys,'b.',xx,interp1(xs,ys,xx,'pchip'),'b-',-xx,interp1(xs,ys,xx,'pchip'),'b-'),
        #       ylim([0 1.2])
        #       sfigure(154);
        #       plot(log10(n))
        #       drawnow
    # work out length of flame surface
    w = np.ones(len(xs))
    w[0] = 0.
    w[1] = 55. / 24.
    w[2] = -1. / 6.
    w[3] = 11. / 8.
    w[-1] = 0.
    w[-2] = 55. / 24.
    w[-3] = -1. / 6.
    w[-4] = 11. / 8.
    Dys = derivs.diff_FD(ys, N, dx)
    D2ys = derivs.diff2_FD(ys, N, dx)
    temp = np.sqrt(1. + np.power(Dys, 2.))
    k = -D2ys / np.power(temp, 3.)
    if axi:
        temp2 = np.zeros(len(temp))
        temp2[1:-1] = -Dys[1:-1] / (temp[1:-1] * xs[1:-1])
        temp2[-1] = 0   # here we have r =0
        k = k + temp2

    if axi:
        # area = sum(2*pi*0.5*(xs(2:end)+xs(1:end-1)).*sqrt( diff(xs).^2 + diff(ys).^2))%sqrt(1 + (Ds*ys).^2))
        integ = xs * np.sqrt(1 + np.power(Dys, 2.)) * (1 + Mark * k)
        area = -np.dot(w, integ) * (xs[1] - xs[0])
    else:
        # area = sum(sqrt( diff(xs).^2 + diff(ys).^2))
        integ = np.sqrt(1 + np.power(Dys, 2.)) * (1 + Mark * k)
        area = -np.dot(w, integ) * (xs[1] - xs[0])

    print('Steady_flame area routine %f: ' % area)
    return [area, xs, ys]


def newton_step(Mark, beta, axi, dx, ys, N):
    """
    Newton-Raphson step
    :param Mark: Markstein number = Markstein length/Flame height
    :param beta: flame aspect ratio
    :param axi: 1 = axisymmetric, 0 = not
    :param dx: grid space
    :param ys: current shape
    :param N: number of points on the flame.
    :return: updated shape ys
    """
    SL = 1. / np.sqrt(1 + beta * beta)
    r = grad(ys, N, dx, axi, SL, Mark, beta)
    [a, b, c] = make_Jac(ys, N, dx, 1e-8, axi, SL, Mark, beta)
    dys = tridiag(N, a, b, c, r[1:], np.zeros(len(b)))
    ys = ys - np.insert(dys, 0, 0)

    return ys


def tridiag(N, a, b, c, r, u):
    """
    Operations on a tridiagonal matrix
    :param N: number of points on the flame.
    :param a: Upper diagonal values
    :param b: Diagonal values
    :param c: Lower diagonal values
    :param r: residual
    :param u: empty array that will be returned. Useless thing in Python.
    :return: u = J*r = correction to flame shape
    """
    n = (N + 1) / 2 - 1
    gam = np.zeros(n)
    bet = b[0]
    u[0] = r[0] / bet
    for j in range(1, n):
        gam[j] = c[j - 1] / bet
        bet = b[j] - a[j] * gam[j]
        if abs(bet) < 1e-14:
            print 'error in tridag!!'

        u[j] = (r[j] - a[j] * u[j - 1]) / bet
    for j in range(n - 2, -1, -1):
        u[j] = u[j] - gam[j + 1] * u[j + 1]

    return u


def make_Jac(ys, N, dx, delta, axi, SL, Mark, beta):
    """
    Build Jacobian Matrix. With a second order finite difference it is tridiagonal
    :param ys: current flame shape
    :param N: number of points on the flame.
    :param dx: grid space
    :param delta: small number for evaluation of gradients
    :param axi: 1 if axi, 0 if not
    :param SL: flame speed
    :param Mark: Markstein number
    :param beta: flame aspect ratio
    :return: list containing Jacobian, upper, central and lower diagonals
    """
    n = (N + 1) / 2
    dys_dt = grad(ys, N, dx, axi, SL, Mark, beta)
    J = np.empty((n - 1, n - 1))
    a = np.empty(n - 1)
    c = np.empty(n - 1)

    for i in range(0, n - 1):
        d = np.zeros(len(ys) - 1)
        d[i] = delta
        temp = grad(ys + np.insert(d, 0, 0), N, dx, axi, SL, Mark, beta)
        J[:, i] = (temp[1:] - dys_dt[1:]) / delta

    a[0] = 0
    for i in range(0, n - 2):
        a[i + 1] = J[i + 1, i]
        c[i] = J[i, i + 1]

    b = np.diag(J)
    return [a, b, c]


def grad(ys, N, dx, axi, SL, Mark, beta):
    """
    Evaluate residual from the steady, single-valued G-equation
    :param ys: current flame shape
    :param N: number of points on the flame.
    :param dx: grid space
    :param axi: 1 if axi, 0 if not
    :param SL: flame speed
    :param Mark: Markstein number
    :param beta: flame aspect ratio
    :return: res = 1 - SL * (1 - Mark * k)*sqrt(1+Dys^2)
    """
    n = (N + 1) / 2
    Dys = derivs.diff_FD(ys, N, dx)
    D2ys = derivs.diff2_FD(ys, N, dx)
    temp = np.sqrt(1 + np.power(Dys, 2.))
    k = D2ys / np.power(temp, 3.)
    if axi:
        xs = np.linspace(1, 0, n) / beta
        temp2 = np.zeros(len(temp))
        temp2[0:-1] = + Dys[0:-1] / (temp[0:-1] * xs[0:-1])
        temp2[-1] = 0    # here we have r=0
        k += temp2

    dys_dt = (1 - SL * (1 - Mark * k) * temp)  # Markstein not multiplied by beta
    dys_dt[0] = 0

    return dys_dt