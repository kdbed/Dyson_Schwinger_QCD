
# correct solution for DSE using Qin-Chang interation - u/d quark
# Calculates A,B, sigs, sigv, Z2
#!python
#cython: boundscheck=False
# distutils: extra_compile_args = -fopenmp
# distutils: extra_link_args = -fopenmp




#####################################################################
# rainbow-ladder approximation  
#####################################################################

import matplotlib.pyplot as plt

#from libc.math cimport  pow, sqrt, log
from cython.parallel import *
import cython
import time
start = time.clock()
from scipy import interpolate
import numpy as np

cdef:
    float E = 2.718281828459045
    float PI = 3.141592653589793



cdef:
    float f = 4.0 #flavors
    float w = 0.5 #Qin param
    float D = 0.82*0.82*0.82/w
    float mu = 2.0 #renorm pt (GeV)
    float Lam = 10.0*mu #UV cutoff
    float reg = 10.0*mu # soft cutoff
    float lamQ = 0.234 #lambda QCD (GeV)
    float lam = 1.0/10.0 #IR cutoff
    float m = 0.0   #0.0178 # renorm. mass
    float mt = 0.5
    float pr = Lam*lam
    float qu = Lam/lam
    int N = 22
    int M = 19
    float gamma = (12.0)/(33.0-2.0*f)
    float tau = pow(E,2) - 1.0
    int i,j,l,v





print "==== Loading Gauss Points and Determining Momenta ====  \n"

kk1 = np.loadtxt("sig_s_Q.dat")
kk2 = np.loadtxt("sig_v_Q.dat")

kk3 = [kk1[i][0] for i in range(400)]
kk4 = [kk2[i][0] for i in range(400)]

kk5 = [kk1[i][1] for i in range(400)]
kk6 = [kk2[i][1] for i in range(400)]

sigs = interpolate.interp1d(kk3, kk5)
sigv = interpolate.interp1d(kk4, kk6)


cdef double[:] x = np.polynomial.legendre.leggauss(N)[0] # points for momentum and weights; log mapping
cdef double[:] weights = np.polynomial.legendre.leggauss(N)[1]
cdef double[:] p = np.zeros(N)      # p values
for i in range(N):
    p[i] = (pr)*pow(qu,x[i])
cdef double[:] k = p                # k values





x0 = np.polynomial.legendre.leggauss(N)[0]
p0 = np.array([(pr)*pow(qu,x0[i]) for i in range(N)])
p1 = [(pr)*pow(qu,x0[i]) for i in range(N)]



cdef double[:] u = np.polynomial.legendre.leggauss(M)[0] # for angle integration -- need no mapping
cdef double[:] wts = np.polynomial.legendre.leggauss(M)[1]


























cdef int s = 10
cdef double complex Complex(double x1, double y1):
    return x1 + 1j*y1
cdef double complex Power(double complex x1, double complex y1):
    return pow(x1,y1)
cdef double complex Log(double complex x1):
    return np.log(x1)
cdef double complex Sqrt(double complex x1):
    return np.sqrt(x1)




cdef double complex f11(double r1, double r2, double y1, double z1):
    return (0.8687323181279998*Sqrt(1 - Power(z1,2))*Power((Sqrt(r2)*(z1 - Complex(0,1)*y1*Sqrt(1 - Power(z1,2))))/Sqrt(r1),s)*
     (1393.096714091478*Power(E,-4.*r1 - 4.*r2 + 8.*Sqrt(r1)*Sqrt(r2)*z1) + (37.899280900183136*(1 - Power(E,-1.*r1 - 1.*r2 + 2.*Sqrt(r1)*Sqrt(r2)*z1)))/((r1 + r2 - 2*Sqrt(r1)*Sqrt(r2)*z1)*Log(-1 + Power(E,2) + Power(1 + 18.262838775659286*(r1 + r2 - 2*Sqrt(r1)*Sqrt(r2)*z1),2))))*
     ((1.*r1 + Sqrt(r1)*Sqrt(r2)*(-2.*z1 - Complex(0.,2.)*y1*Sqrt(1 - Power(z1,2))) + r2*(1. + Power(y1,2)*(2. - 2.*Power(z1,2)) + Complex(0.,2.)*y1*z1*Sqrt(1 - Power(z1,2))))*Power(sigs(r2),2) + 
       r2*(Sqrt(r1)*Sqrt(r2)*((-2. + 8.*Power(y1,2))*z1 - 8.*Power(y1,2)*Power(z1,3) - Complex(0.,2.)*y1*Sqrt(1 - Power(z1,2)) + Complex(0.,8.)*y1*Power(z1,2)*Sqrt(1 - Power(z1,2))) + 
          r1*(1. - Complex(0.,2.)*y1*z1*Sqrt(1 - Power(z1,2)) + Power(y1,2)*(-2. + 2.*Power(z1,2))) + r2*(1. - Complex(0.,4.)*y1*z1*Sqrt(1 - Power(z1,2)) + Power(y1,2)*(-4. + 4.*Power(z1,2))))*Power(sigv(r2),2)))/(r1 + r2 - 2.*Sqrt(r1)*Sqrt(r2)*z1)


cdef double complex f12(double r1, double r2, double y1, double z1):
    return (-3.474929272511999*r2*y1*(-0.25*r1 - 0.75*r2 + 1.*Sqrt(r1)*Sqrt(r2)*z1)*(1 - Power(z1,2))*Power((Sqrt(r2)*(z1 - Complex(0,1)*y1*Sqrt(1 - Power(z1,2))))/Sqrt(r1),s)*
     (Complex(0.,1.)*z1 + 1.*y1*Sqrt(1 - Power(z1,2)))*(1393.096714091478*Power(E,-4.*r1 - 4.*r2 + 8.*Sqrt(r1)*Sqrt(r2)*z1) + 
       (37.899280900183136*(1 - Power(E,-1.*r1 - 1.*r2 + 2.*Sqrt(r1)*Sqrt(r2)*z1)))/((r1 + r2 - 2*Sqrt(r1)*Sqrt(r2)*z1)*Log(-1 + Power(E,2) + Power(1 + 18.262838775659286*(r1 + r2 - 2*Sqrt(r1)*Sqrt(r2)*z1),2))))*
     (-1.*Power(sigs(r2),2) + r2*Power(sigv(r2),2)))/(-0.5*r1 - 0.5*r2 + 1.*Sqrt(r1)*Sqrt(r2)*z1)

cdef double complex f13(double r1, double r2, double y1, double z1):
    return (-6.949858545023998*r2*y1*(-0.25*r1 - 0.75*r2 + 1.*Sqrt(r1)*Sqrt(r2)*z1)*(1 - Power(z1,2))*Power((Sqrt(r2)*(z1 - Complex(0,1)*y1*Sqrt(1 - Power(z1,2))))/Sqrt(r1),s)*
     (Complex(0.,1.)*z1 + 1.*y1*Sqrt(1 - Power(z1,2)))*(1393.096714091478*Power(E,-4.*r1 - 4.*r2 + 8.*Sqrt(r1)*Sqrt(r2)*z1) + 
       (37.899280900183136*(1 - Power(E,-1.*r1 - 1.*r2 + 2.*Sqrt(r1)*Sqrt(r2)*z1)))/((r1 + r2 - 2*Sqrt(r1)*Sqrt(r2)*z1)*Log(-1 + Power(E,2) + Power(1 + 18.262838775659286*(r1 + r2 - 2*Sqrt(r1)*Sqrt(r2)*z1),2))))*
     sigs(r2)*sigv(r2))/(-0.5*r1 - 0.5*r2 + 1.*Sqrt(r1)*Sqrt(r2)*z1)

cdef double complex f21(double r1, double r2, double y1, double z1):
    return (0.8687323181279998*Sqrt(1 - Power(z1,2))*Power((Sqrt(r2)*(z1 - Complex(0,1)*y1*Sqrt(1 - Power(z1,2))))/Sqrt(r1),s)*
     (1393.096714091478*Power(E,-4.*r1 - 4.*r2 + 8.*Sqrt(r1)*Sqrt(r2)*z1) + (37.899280900183136*(1 - Power(E,-1.*r1 - 1.*r2 + 2.*Sqrt(r1)*Sqrt(r2)*z1)))/((r1 + r2 - 2*Sqrt(r1)*Sqrt(r2)*z1)*Log(-1 + Power(E,2) + Power(1 + 18.262838775659286*(r1 + r2 - 2*Sqrt(r1)*Sqrt(r2)*z1),2))))*(Power(Sqrt(r1) - Sqrt(r2)*(z1 - Complex(0,1)*y1*Sqrt(1 - Power(z1,2))),2)*Power(sigs(r2),2) + 
       r2*(r1 - 2*r1*z1*(z1 - Complex(0,1)*y1*Sqrt(1 - Power(z1,2))) + (r1 + 2*r2 - 4*Sqrt(r1)*Sqrt(r2)*z1)*Power(Complex(0,1)*z1 + y1*Sqrt(1 - Power(z1,2)),2))*Power(sigv(r2),2)))/(r1*(r1 + r2 - 2*Sqrt(r1)*Sqrt(r2)*z1))

cdef double complex f22(double r1, double r2, double y1, double z1):
    return (0.8687323181279998*Sqrt(r2)*Sqrt(1 - Power(z1,2))*Power((Sqrt(r2)*(z1 - Complex(0,1)*y1*Sqrt(1 - Power(z1,2))))/Sqrt(r1),s)*(Complex(0,1)*z1 + y1*Sqrt(1 - Power(z1,2)))*
     (Complex(0,2)*(-(Sqrt(r1)*r2) + r1*Sqrt(r2)*z1) + Sqrt(r2)*(r1 + 3*r2 - 4*Sqrt(r1)*Sqrt(r2)*z1)*(Complex(0,1)*z1 + y1*Sqrt(1 - Power(z1,2))))*
     (1393.096714091478*Power(E,-4.*r1 - 4.*r2 + 8.*Sqrt(r1)*Sqrt(r2)*z1) + (37.899280900183136*(1 - Power(E,-1.*r1 - 1.*r2 + 2.*Sqrt(r1)*Sqrt(r2)*z1)))/((r1 + r2 - 2*Sqrt(r1)*Sqrt(r2)*z1)*Log(-1 + Power(E,2) + Power(1 + 18.262838775659286*(r1 + r2 - 2*Sqrt(r1)*Sqrt(r2)*z1),2))))*(-Power(sigs(r2),2) + r2*Power(sigv(r2),2)))/(r1*(r1 + r2 - 2*Sqrt(r1)*Sqrt(r2)*z1))

cdef double complex f23(double r1, double r2, double y1, double z1):
    return (1.7374646362559996*Sqrt(r2)*Sqrt(1 - Power(z1,2))*Power((Sqrt(r2)*(z1 - Complex(0,1)*y1*Sqrt(1 - Power(z1,2))))/Sqrt(r1),s)*(Complex(0,1)*z1 + y1*Sqrt(1 - Power(z1,2)))*
     (Complex(0,2)*(-(Sqrt(r1)*r2) + r1*Sqrt(r2)*z1) + Sqrt(r2)*(r1 + 3*r2 - 4*Sqrt(r1)*Sqrt(r2)*z1)*(Complex(0,1)*z1 + y1*Sqrt(1 - Power(z1,2))))*
     (1393.096714091478*Power(E,-4.*r1 - 4.*r2 + 8.*Sqrt(r1)*Sqrt(r2)*z1) + (37.899280900183136*(1 - Power(E,-1.*r1 - 1.*r2 + 2.*Sqrt(r1)*Sqrt(r2)*z1)))/((r1 + r2 - 2*Sqrt(r1)*Sqrt(r2)*z1)*Log(-1 + Power(E,2) + Power(1 + 18.262838775659286*(r1 + r2 - 2*Sqrt(r1)*Sqrt(r2)*z1),2))))*sigs(r2)*sigv(r2))/(r1*(r1 + r2 - 2*Sqrt(r1)*Sqrt(r2)*z1))





cdef double complex f31(double r1, double r2, double y1, double z1):
    return (-2.606196954383999*Sqrt(r2)*Sqrt(1 - Power(z1,2))*Power((Sqrt(r2)*(z1 - Complex(0,1)*y1*Sqrt(1 - Power(z1,2))))/Sqrt(r1),s)*(z1 - Complex(0.,1.)*y1*Sqrt(1 - Power(z1,2)))*
     (1393.096714091478*Power(E,-4.*r1 - 4.*r2 + 8.*Sqrt(r1)*Sqrt(r2)*z1) + (37.899280900183136*(1 - Power(E,-1.*r1 - 1.*r2 + 2.*Sqrt(r1)*Sqrt(r2)*z1)))/((r1 + r2 - 2*Sqrt(r1)*Sqrt(r2)*z1)*Log(-1 + Power(E,2) + Power(1 + 18.262838775659286*(r1 + r2 - 2*Sqrt(r1)*Sqrt(r2)*z1),2))))*sigs(r2)*sigv(r2))/Sqrt(r1)


cdef double complex f32(double r1, double r2, double y1, double z1):
    return (-5.212393908767998*Power(r2,1.5)*Sqrt(1 - Power(z1,2))*Power((Sqrt(r2)*(z1 - Complex(0,1)*y1*Sqrt(1 - Power(z1,2))))/Sqrt(r1),s)*(z1 - Complex(0.,1.)*y1*Sqrt(1 - Power(z1,2)))*
     (1393.096714091478*Power(E,-4.*r1 - 4.*r2 + 8.*Sqrt(r1)*Sqrt(r2)*z1) + (37.899280900183136*(1 - Power(E,-1.*r1 - 1.*r2 + 2.*Sqrt(r1)*Sqrt(r2)*z1)))/((r1 + r2 - 2*Sqrt(r1)*Sqrt(r2)*z1)*Log(-1 + Power(E,2) + Power(1 + 18.262838775659286*(r1 + r2 - 2*Sqrt(r1)*Sqrt(r2)*z1),2))))*sigs(r2)*sigv(r2))/Sqrt(r1)
    
cdef double complex f33(double r1, double r2, double y1, double z1):
    return (2.606196954383999*Sqrt(r2)*Sqrt(1 - Power(z1,2))*Power((Sqrt(r2)*(z1 - Complex(0,1)*y1*Sqrt(1 - Power(z1,2))))/Sqrt(r1),s)*(z1 - Complex(0.,1.)*y1*Sqrt(1 - Power(z1,2)))*
     (1393.096714091478*Power(E,-4.*r1 - 4.*r2 + 8.*Sqrt(r1)*Sqrt(r2)*z1) + (37.899280900183136*(1 - Power(E,-1.*r1 - 1.*r2 + 2.*Sqrt(r1)*Sqrt(r2)*z1)))/((r1 + r2 - 2*Sqrt(r1)*Sqrt(r2)*z1)*Log(-1 + Power(E,2) + Power(1 + 18.262838775659286*(r1 + r2 - 2*Sqrt(r1)*Sqrt(r2)*z1),2))))*(-1.*Power(sigs(r2),2) + r2*Power(sigv(r2),2)))/Sqrt(r1)


print "===== Determining Angular Integrals for Each k^2,p^2  ===== \n "



# angle integral for each p,k
cdef double[:,:] zint11 = np.zeros((N,N))
cdef double[:,:] zint12 = np.zeros((N,N))
cdef double[:,:] zint13 = np.zeros((N,N))
cdef double[:,:] zint21 = np.zeros((N,N))
cdef double[:,:] zint22 = np.zeros((N,N))
cdef double[:,:] zint23 = np.zeros((N,N))
cdef double[:,:] zint31 = np.zeros((N,N))
cdef double[:,:] zint32 = np.zeros((N,N))
cdef double[:,:] zint33 = np.zeros((N,N))





import scipy
'''
cdef int nnn = N*N
cdef int count = 0
#with nogil, parallel(num_threads=8):
for j in range(N):      # k 
    for i in range(N):
        print nnn - count
        for l in range(M): # angle integral
            for v in range(M):
                zint11[i,j] += wts[l]*wts[v]*f11(p[i],k[j],u[l],u[v]).real
                zint12[i,j] += wts[l]*wts[v]*f12(p[i],k[j],u[l],u[v]).real
                zint13[i,j] += wts[l]*wts[v]*f13(p[i],k[j],u[l],u[v]).real
                zint21[i,j] += wts[l]*wts[v]*f21(p[i],k[j],u[l],u[v]).real
                zint22[i,j] += wts[l]*wts[v]*f22(p[i],k[j],u[l],u[v]).real
                zint23[i,j] += wts[l]*wts[v]*f23(p[i],k[j],u[l],u[v]).real
                zint31[i,j] += wts[l]*wts[v]*f31(p[i],k[j],u[l],u[v]).real
                zint32[i,j] += wts[l]*wts[v]*f32(p[i],k[j],u[l],u[v]).real
                zint33[i,j] += wts[l]*wts[v]*f33(p[i],k[j],u[l],u[v]).real
        count += 1
'''



gfun = lambda xt : -1
hfun = lambda yt : 1

cdef int nnn = N*N
cdef int count = 0 
cdef double epsi = 1.45e-5

for j in range(N):
    for i in range(N):
        print nnn - count
        zint11[i,j] = scipy.integrate.dblquad(lambda xa,ya: f11(p[i],k[j],xa,ya).real, -1, 1, gfun, hfun, epsabs=epsi)[0]
        zint12[i,j] = scipy.integrate.dblquad(lambda xa,ya: f12(p[i],k[j],xa,ya).real, -1, 1, gfun, hfun,epsabs=epsi)[0]
        zint13[i,j] = scipy.integrate.dblquad(lambda xa,ya: f13(p[i],k[j],xa,ya).real, -1, 1, gfun, hfun,epsabs=epsi)[0]
        zint21[i,j] = scipy.integrate.dblquad(lambda xa,ya: f21(p[i],k[j],xa,ya).real, -1, 1, gfun, hfun,epsabs=epsi)[0]
        zint22[i,j] = scipy.integrate.dblquad(lambda xa,ya: f22(p[i],k[j],xa,ya).real, -1, 1, gfun, hfun,epsabs=epsi)[0]
        zint23[i,j] = scipy.integrate.dblquad(lambda xa,ya: f23(p[i],k[j],xa,ya).real, -1, 1, gfun, hfun,epsabs=epsi)[0]
        zint31[i,j] = scipy.integrate.dblquad(lambda xa,ya: f31(p[i],k[j],xa,ya).real, -1, 1, gfun, hfun,epsabs=epsi)[0]
        zint32[i,j] = scipy.integrate.dblquad(lambda xa,ya: f32(p[i],k[j],xa,ya).real, -1, 1, gfun, hfun,epsabs=epsi)[0]
        zint33[i,j] = scipy.integrate.dblquad(lambda xa,ya: f33(p[i],k[j],xa,ya).real, -1, 1, gfun, hfun,epsabs=epsi)[0]
        count += 1











# define initial guess for A,B

f1 = np.zeros(N)
f2 = np.zeros(N)
f3 = np.zeros(N)
G = np.zeros(N)
H = np.zeros(N)
W = np.zeros(N)



z2 = 0.80955
m0 = m

epsilon = 0
delta = 0
count = 0





print "=== Iterate over momenta integrals to find A,B === \n"

# integrations: loop in i is for p-values, loop in j does integration for given p

#delta = 0
#while epsilon == 0 or delta == 0:
while epsilon < 100:
    for i in range(N):           # for given p-value
        s1 = 0.0
        s2 = 0.0
        s3 = 0.0
        for j in range(N):       # for given k-value
            s1 += weights[j]*k[j]*k[j]*Log(qu)*(zint11[i,j]*f1[j] + f2[j]*zint12[i,j] + f3[j]*zint13[i,j])
            s2 += weights[j]*k[j]*k[j]*Log(qu)*(zint21[i,j]*f1[j] + f2[j]*zint22[i,j] + f3[j]*zint23[i,j])
            s3 += weights[j]*k[j]*k[j]*Log(qu)*(zint31[i,j]*f1[j] + f2[j]*zint32[i,j] + f3[j]*zint33[i,j])
        G[i] = z2 + (0.5*pow(2.0*PI,-3)*s1).real
        H[i] = (0.5*pow(2.0*PI,-3)*s2).real
        W[i] = (0.5*pow(2.0*PI,-3)*s3).real
    #delta = np.allclose(A,T,1e-15,0)
    #epsilon = np.allclose(B,R,1e-15,0)
    for i in range(N):
        f1[i] = G[i]
        f2[i] = H[i]
        f3[i] = W[i]
    epsilon += 1
    print epsilon


print "==== Finished === \n"

 


a_p_Q = open('f1m10.dat','w+')
for i in range(N):
    a_p_Q.write("%10.20f     %10.20f \n" % (p[i],f1[i]))
a_p_Q.close()
b_p_Q = open('f2m10.dat', 'w+')
for i in range(N):
    b_p_Q.write("%10.20f     %10.20f \n" % (p[i],f2[i]))
b_p_Q.close()
sig_v_Q = open('f3m10.dat','w+')
for i in range(N):
    sig_v_Q.write("%10.20f     %10.20f \n" % (p[i],f3[i]))
sig_v_Q.close()




