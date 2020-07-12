############################################################################
#
# Functions to calculate integration points and weights for Gaussian
# quadrature
#
# x, w = gaussxw(N) returns integration points x and integration weights
#                      w s.t. sum_i w[i]*f(x[i]) is the Nth-order Gaussian
#                    approximation to the integral int_{-1}^1 f(x) dx
#
# This code finds the zeros of the nth Legendre polynomial using Newton's 
# method.
#
#
# Taken from code written by Mark Newman <mejn@umich.edu>.
#
############################################################################



from numpy import ones, copy, cos, tan, pi, linspace


def gaussxw(N):


    # Initial approximation to roots of Legendre polynomial (given in Abramowitz a          nd Stegun 22.16.6
    a = linspace(3, 4*N-1, N)/(4*N+2)
    x = cos(pi*a+1/(8*N*N*tan(a)))


    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta>epsilon:
        p0 = ones(N,float)
        p1 = copy(x)
        for k in range (1,N):
            p0,p1 = p1, ((2*k+1)*x*p1-k*p0)/(k+1)
        dp = (N+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(abs(dx))

    # Calculate the weights
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)

    return x,w
