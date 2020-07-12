
# correct solution for DSE using Qin-Chang interation - u/d quark
# Calculates A,B, sigs, sigv, Z2
#!python
#cython: boundscheck=False
# distutils: extra_compile_args = -fopenmp
# distutils: extra_link_args = -fopenmp




#####################################################################
# rainbow-ladder approximation  
#####################################################################

import numpy as np
import matplotlib.pyplot as plt
cimport numpy as np
from libc.math cimport  pow, sqrt, log, exp
from cython.parallel import *
import cython
import time
start = time.clock()
from scipy import interpolate

cdef:
    float E = 2.718281828459045
    float PI = 3.141592653589793



cdef:
    float f = 4.0 #flavors
    float w = 0.5 #Qin param
    float D = 0.82*0.82*0.82/w
    float mu = 2.0 #renorm pt (GeV)
    float Lam = 1000.0      #100.0*mu #UV cutoff
    float reg = 1000000.0*mu # soft cutoff
    float lamQ = 0.234 #lambda QCD (GeV)
    float lam = 1.0/1000.0 #IR cutoff
    float m = 0.0178   #0.0178 # renorm. mass
    float mt = 0.5
    float pr = Lam*lam
    float qu = Lam/lam
    int N = 270
    int M = 180
    float gamma = (12.0)/(33.0-2.0*f)
    float tau = pow(E,2) - 1.0
    int i,j,l




print "==== Loading Gauss Points and Determining Momenta ====  \n"


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]



cdef double[:] x = np.polynomial.legendre.leggauss(N)[0] # points for momentum and weights; log mapping
cdef double[:] weights = np.polynomial.legendre.leggauss(N)[1]
cdef double[:] p = np.zeros(N)      # p values
for i in range(N):
    p[i] = (pr)*pow(qu,x[i])
cdef double[:] k = p                # k values





x0 = np.polynomial.legendre.leggauss(N)[0]
p0 = np.array([(pr)*pow(qu,x0[i]) for i in range(N)])
p1 = [(pr)*pow(qu,x0[i]) for i in range(N)]
cdef:
    int index = p1.index(find_nearest(p0,mu*mu))
print p0[index]


cdef double[:] u = np.polynomial.legendre.leggauss(M)[0] # for angle integration -- need no mapping
cdef double[:] wts = np.polynomial.legendre.leggauss(M)[1]




















# Define kernel (qin-chang interaction. paper: interaction model for gap equation)
@cython.boundscheck(False)
cdef double F(double var1) nogil:
    return (1.0/var1)*(1.0 - pow(E, -var1/(4.0*mt*mt)))
@cython.boundscheck(False)
cdef double alpha(double var2) nogil:
    return (2.0/pow(w,4))*D*pow(E,-var2/(w*w)) +(2.0*gamma*F(var2))/(log(tau+pow(1.0+var2/(lamQ*lamQ),2)))







print "===== Determining Angular Integrals for Each k^2,p^2  ===== \n "



# angle integral for each p,k
cdef double[:,:] zaint = np.zeros((N,N))
cdef double[:,:] zbint = np.zeros((N,N))




cdef inline double q(double p,double k,double y) nogil:
    return p + k - 2.0*sqrt(k)*sqrt(p)*y





with nogil, parallel(num_threads=8):
    for j in prange(N, schedule='dynamic'):      # k 
        for i in range(N):
            for l in range(M): # angle integral
                zbint[i,j] += wts[l]*alpha(q(p[i],k[j],u[l]))*sqrt(1-u[l]*u[l])
                zaint[i,j] += wts[l]*alpha(q(p[i],k[j],u[l]))*sqrt(1-u[l]*u[l])*(sqrt(p[i])*sqrt(k[j])*u[l]+2*(p[i]-sqrt(p[i])*sqrt(k[j])*u[l])*(sqrt(p[i])*sqrt(k[j])*u[l]-k[j])/(q(p[i],k[j],u[l])))














# define initial guess for A,B

B = [1. for i in range(N)]
A = [1.0 for i in range(N)]

G = np.zeros(N)
H = np.zeros(N)




z2 = 1.0
m0 = m

epsilon = 0
delta = 0
count = 0


print "=== Iterate over momenta integrals to find A,B === \n"

# integrations: loop in i is for p-values, loop in j does integration for given p

delta = 0
while epsilon == 0 or delta == 0:
#while epsilon < 1000:
    for i in range(N):           # for given p-value
        s = 0.0
        r = 0.0
        for j in range(N):       # for given k-value
             s += weights[j]*k[j]*k[j]*log(qu)*zaint[i,j]*(A[j]/((k[j])*(A[j])*A[j]+(B[j])*B[j]))*(reg**2/(k[j]**2 + reg**2))
             r += weights[j]*k[j]*k[j]*log(qu)*zbint[i,j]*(B[j]/((k[j])*(A[j])*A[j]+(B[j])*B[j]))*(reg**2/(k[j]**2 + reg**2))
        G[i] = (2.0/(3.0*(PI)))*s*(1.0/p[i])*z2*z2
        H[i] = (2.0/PI)*r*z2*z2
    oldm = m0
    oldz = z2
    zcalc = interpolate.interp1d(k,G)
    m0calc = interpolate.interp1d(k,H)
    z2 = 1.0- zcalc(4.0)
    m0 = (m-m0calc(4.0))/(z2)
    T = [A[i] for i in range(N)]
    R = [B[j] for j in range(N)]  
    # match the boundary conditions for new A,B
    A = [z2+G[i] for i in range(N)]
    B = [z2*m0+H[j] for j in range(N)]
    delta = np.allclose(A,T,1e-11,0)
    epsilon = np.allclose(B,R,1e-11,0)
    count += 1
    print count, z2


print "==== Finished === \n"

print "Function A at the renorm. point = %2.2f. \n" % A[index]
print "Function B at the renorm. point = %2.2f. \n" % B[index]
print 'z2 = ', z2
print 'm0 = ', m0
#z4 = (z2*m0)/(m)
#print 'z4 = ', z4




I = np.zeros(N)
for i in range(N):
    I[i] =B[i]/A[i]


sigv = [A[i]/(p[i]*pow(A[i],2) + pow(B[i],2)) for i in range(N)]
sigs = [B[i]/(p[i]*pow(A[i],2) + pow(B[i],2)) for i in range(N)]


'''
x = np.linspace(lam*lam + lam*lam/1000.0,Lam*Lam - Lam*Lam/1000.0,num = 10000000)
ySS = interpolate.interp1d(p, sigs, kind = 'linear')
yVV = interpolate.interp1d(p, sigv, kind = 'linear')
yS = ySS(x)
yV = yVV(x)

dys = [0.0]*len(x)
dys[0] = (yS[0]-yS[1])/(x[0]-x[1])
for i in range(1,len(yS)-1):
    dys[i] = (yS[i+1]-yS[i-1])/(x[i+1]-x[i-1])
dys[-1] = (yS[-1]-yS[-2])/(x[-1]-x[-2])


dyv = [0.0]*len(x)
dyv[0] = (yV[0]-yV[1])/(x[0]-x[1])
for i in range(1,len(yV)-1):
    dyv[i] = (yV[i+1]-yV[i-1])/(x[i+1]-x[i-1])
dyv[-1] = (yV[-1]-yV[-2])/(x[-1]-x[-2])
'''
 


a_p_Q = open('a_p_Q.dat','w+')
for i in range(N):
    a_p_Q.write("%10.20f     %10.20f \n" % (p[i],A[i]))
a_p_Q.close()
b_p_Q = open('b_p_Q.dat', 'w+')
for i in range(N):
    b_p_Q.write("%10.20f     %10.20f \n" % (p[i],B[i]))
b_p_Q.close()
sig_v_Q = open('sig_v_Q.dat','w+')
for i in range(N):
    sig_v_Q.write("%10.20f     %10.20f \n" % (p[i],sigv[i]))
sig_v_Q.close()
sig_s_Q = open('sig_s_Q.dat', 'w+')
for i in range(N):
    sig_s_Q.write("%10.20f     %10.20f \n" % (p[i],sigs[i]))
sig_s_Q.close()
'''ssp = open('ssp_Q.dat', 'w+')
for i in range(len(dys)):
    ssp.write("%10.20f     %10.20f \n" % (x[i],dys[i]))
ssp.close()
ssv = open('ssv_Q.dat', 'w+')
for i in range(len(dyv)):
    ssv.write("%10.20f     %10.20f \n" % (x[i],dyv[i]))
ssv.close()
'''




end = time.clock()
print(end-start)
