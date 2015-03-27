__author__ = 'Alessandro Orchini'

#! /usr/bin/python
"""Functions that evaluate steady flame area and properties. Based on a code originally written by ICW"""

import numpy as np
import derivs

def steady_flame(Mark,beta,axi,N,tol,ys):
    """Calculate Flame Steady Component"""
    SL = 1./np.sqrt(1.+beta*beta)
    n = (N+1)/2
    xs = np.linspace(1, 0, n)/beta
    dx = xs[1]-xs[0]
    #xx = np.linspace(1./beta, 0, 500)

    #loop
    i=1
    itmax = 99
    r = grad(ys,N,dx,axi,SL,Mark,beta)
    n = [np.linalg.norm(r*(-dx))]
    #  sfigure(152);plot(xs,ys,'b.',xx,interp1(xs,ys,xx,'pchip'),-xx,interp1(xs,ys,xx,'pchip'),'b-'),ylim([0 1.2])
    #  sfigure(153);plot(r)
    while i < itmax and np.linalg.norm(r*(-dx)) > tol:
         ys = newton_step(Mark, beta, axi, dx, ys, N)
         r = grad(ys, N, dx, axi, SL, Mark, beta)
         n.append(np.linalg.norm(r*(-dx)))
        #         if (1)
        #            sfigure(152);plot(xs,ys,'b.',xx,interp1(xs,ys,xx,'pchip'),-xx,interp1(xs,ys,xx,'pchip'),'b-'),ylim([0 1.2])
        #            sfigure(153);plot(r)
        #            drawnow
        #         end
         i =+ 1
        #          pause
        #    disp(['The steady flame has converged at iteration ',num2str(i),' out of ',num2str(itmax),' with a residual equal to ',num2str(n(end))])
        #    sfigure(152);plot(xs,ys,'b.',xx,interp1(xs,ys,xx,'pchip'),'b-',-xx,interp1(xs,ys,xx,'pchip'),'b-'),ylim([0 1.2])
        #    sfigure(154);plot(log10(n))
        #    drawnow
    #work out length of flame surface
    w = np.ones(len(xs))
    w[0]=0.
    w[1]=55./24.
    w[2]=-1./6.
    w[3]=11./8.
    w[-1]=0.
    w[-2]=55./24.
    w[-3]=-1./6.
    w[-4]=11./8.
    Dys = derivs.diff_FD(ys, N, dx)
    D2ys = derivs.diff2_FD(ys, N, dx)
    temp = np.sqrt(1.+np.power(Dys, 3.))
    k = -D2ys/np.power(temp, 3.)
    if axi:
        temp2 = -Dys/(temp*xs)
        temp2[-1]=0
        k = k +temp2

    if (axi):
        #area = sum(2*pi*0.5*(xs(2:end)+xs(1:end-1)).*sqrt( diff(xs).^2 + diff(ys).^2))%sqrt(1 + (Ds*ys).^2))
        integ = xs*np.sqrt( 1 + np.power(Dys,2.))*(1+Mark*k)
        area = -np.dot(w*integ)*(xs[1]-xs[0])
    else:
        #area = sum(sqrt( diff(xs).^2 + diff(ys).^2))
        integ = np.sqrt( 1 + np.power(Dys,2.))*(1+Mark*k)
        area = -np.dot(w*integ)*(xs[1]-xs[0])

    print('Steady_flame area routine %f: ' % area)
    return [area, xs, ys]


def newton_step(Mark, beta, axi, dx, ys, N):
    """Newton-Raphson step"""
    SL = 1./np.sqrt(1+beta*beta)
    n = (N+1)/2
    r = grad(ys, N, dx, axi, SL, Mark, beta)
    [J,a,b,c] = make_Jac(ys, N, dx, 1e-8, axi, SL, Mark, beta)
    dys = tridiag(N, a, b, c, r[1:], np.zeros(len(b)))
    ys = ys - [0, dys]

    return ys


def tridiag(N, a, b, c, r, u):
    """Operations on a tri-diagonal Jacobian. Reference Numerical Recipes"""
    n = (N+1)/2 - 1
    gam = np.zeros(n)
    bet = b(1)
    u = np.zeros(n)
    u[0] = r[0]/bet
    for j in range(1, n):
       gam[j] = c[j-1]/bet
       bet = b[j]-a[j]*gam[j]
       if abs(bet)<1e-14:
           print 'error in tridag!!'

       u[j] = (r[j]-a[j]*u[j-1])/bet
    for j in range (n-2, -1, -1):
       u[j] = u[j] - gam[j+1]*u[j+1]

    return u

def make_Jac(ys, N, dx, delta, axi, SL, Mark, beta):
    n = (N+1)/2
    dys_dt = grad(ys, N, dx, axi, SL, Mark, beta)
    for i in range (0,n-1):
        d = zeros(size(ys(2:end)))
        d(i) = delta
        temp = grad(ys + [0;d],N,dx,axi,SL,Mark,beta)
        J(:,i) = ( temp(2:end) - dys_dt(2:end))/delta

    a(1)=0
    for i=1:n-2
        a(i+1) = J(i+1,i)
        c(i) = J(i,i+1)

    b = diag(J)
    return [J,a,b,c]

def grad(ys, N, dx, axi, SL, Mark, beta):
    n = (N+1)/2
    Dys = diff_FD(ys,N,dx)
    D2ys = diff2_FD(ys,N,dx)
    temp = sqrt(1+(Dys).^2)
    k = + (D2ys)./temp.^3
    if (axi)
        xs = linspace(1,0,n)'/beta
        temp2 = + (Dys)./(temp.*xs)
        temp2(end)=0;  %%here we have r=0
        k = k +temp2
    end
    dys_dt = (1 - SL*(1-Mark*k).*temp)        %Markstein not multiplied by beta
    dys_dt(1)=0

    return dys_dt