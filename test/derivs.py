__author__ = 'Alessandro Orchini'

function [Dys] = diff_FD(ys,N,dx)
    n = (N+1)/2
    Dys = zeros(size(ys))
    Dys(1) = (-1.5*ys(1) + 2*ys(2) -0.5*ys(3))/dx
    for i=2:n-1
        Dys(i) = (-0.5*ys(i-1) + 0.5*ys(i+1))/dx
    end
    Dys(n) = 0
end

function [D2ys] = diff2_FD(ys,N,dx)
    n = (N+1)/2
    D2ys = zeros(size(ys))
    D2ys(1) = (1.0*ys(1) - 2.0*ys(2) +1.0*ys(3))/(dx*dx)
    for i=2:n-1
        D2ys(i) = (1.0*ys(i-1) - 2.0*ys(i) + 1.0*ys(i+1))/(dx*dx)
    end
    D2ys(n) = (2.0*ys(n-1) - 2.0*ys(n))/(dx*dx)

end