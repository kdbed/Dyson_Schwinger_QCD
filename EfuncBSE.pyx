from scipy.interpolate import interp1d
from mpi4py import MPI
import time
from numpy import ndarray
from numpy import linalg
from numpy import dot
from numpy cimport *
import cmath
from cmath import *
from scipy import optimize
import random
from numpy import linspace
from numpy import zeros
from math import factorial
from numpy import mat
from numpy import reshape
Zq = [0.3768+0.7116j,0.3768 -0.7116j,0.1381,0.1381]
mq = [0.7112 +  0.2228j,0.7112 -  0.2228j, -0.7788 +  0.7548j, -0.7788 - 0.7548j]
Zg = [3.03287 +  61.6236j,3.03287 - 61.6236j, -3.04202 - 13.7662j, -3.04202 +13.7662j]
mg = [0.553342 + 0.260293j, 0.553342 - 0.260293j,0.579953 + 0.770662j, 0.579953 - 0.770662j]
cdef float cc = 5.86845



cdef int np = 12

x1 = [0.00921968, 0.0479414, 0.115049, 0.206341, 0.316084, 0.437383, 0.562617, 0.683916, 0.793659, 0.884951, 0.952059, 0.99078]
x2  =  [0.00921968, 0.0479414, 0.115049, 0.206341, 0.316084, 0.437383,0.562617, 0.683916, 0.793659, 0.884951, 0.952059, 0.99078]
x3 =  [0.00921968, 0.0479414, 0.115049, 0.206341, 0.316084, 0.437383,0.562617, 0.683916, 0.793659, 0.884951, 0.952059, 0.99078]
x4 = [0.00921968, 0.0479414, 0.115049, 0.206341, 0.316084, 0.437383,0.562617, 0.683916, 0.793659, 0.884951, 0.952059, 0.99078]
x5 = [0.00921968, 0.0479414, 0.115049, 0.206341, 0.316084, 0.437383,0.562617, 0.683916, 0.793659, 0.884951, 0.952059, 0.99078]
w1  = [0.0235877, 0.0534697, 0.0800392, 0.101584, 0.116746, 0.124574, 0.124574, 0.116746, 0.101584, 0.0800392, 0.0534697, 0.0235877]
w2 = [0.0235877, 0.0534697, 0.0800392, 0.101584, 0.116746, 0.124574, 0.124574, 0.116746, 0.101584, 0.0800392, 0.0534697, 0.0235877]
w3 = [0.0235877, 0.0534697, 0.0800392, 0.101584, 0.116746, 0.124574, 0.124574, 0.116746, 0.101584, 0.0800392, 0.0534697, 0.0235877]
w4 = [0.0235877, 0.0534697, 0.0800392, 0.101584, 0.116746, 0.124574, 0.124574, 0.116746, 0.101584, 0.0800392, 0.0534697, 0.0235877]
w5 = [0.0235877, 0.0534697, 0.0800392, 0.101584, 0.116746, 0.124574, 0.124574, 0.116746, 0.101584, 0.0800392, 0.0534697, 0.0235877]



##############################################################################
cdef d1(int i,double lamg,double lamh,double P,int g, int h):
    return x5[i]*lamg**2+(1.0-x5[i])*lamh**2
cdef d2(int i,double lamg,double lamh,double P,int g, int h):
    return (x5[i])**(g-1)*(1.0-x5[i])**(h-1)*factorial(g+h-1)/(factorial(h-1)*factorial(g-1))
cdef LEgEh(int i,double lamg,double lamh,double P,int g, int h):
    return (w5[i]*d2(i,lamg,lamh,P,g,h))/((4.0*pi)**2*(h+g-1.0)*(h+g-2.0)*(d1(i,lamg,lamh,P,g,h))**(g+h-2))
cdef LTgEh(int i,double lamg,double lamh,double P,int g, int h):
    return (w5[i]*d2(i,lamg,lamh,P,g,h)*P)/(2.0*(4.0*pi)**2*(h+g-1.0)*(h+g-2.0)*(g+h-3.0)*(d1(i,lamg,lamh,P,g,h))**(g+h-3))
cdef LTgTh(int i,double lamg,double lamh,double P,int g, int h):
    return (w5[i]*d2(i,lamg,lamh,P,g,h)*3.0*P**2)/(4.0*(4.0*pi)**2*(h+g-1.0)*(h+g-2.0)*(h+g-3.0)*(g+h-4.0)*(d1(i,lamg,lamh,P,g,h))**(g+h-4))



##############################################################################
##############################################################################



cdef L11(double lam,double lamm,double P):
        L11a = 0.0
        cdef int i=0
        for  i from 0 <= i < np:
               L11a = L11a + LEgEh(i,lam,lam,P,2,1)
        return L11a

cdef L12(double lam,double lamm,double P):
        L13a = 0.0
        cdef int i=0
        for  i from 0 <= i < np:
               L13a = L13a + LTgEh(i,lam,lam,P,2,3)
        return L13a        

cdef L13(double lam,double lamm,double P):
        L15a = 0.0
        cdef int i=0
        for  i from 0 <= i < np:
               L15a = L15a + LEgEh(i,lam,lam,P,2,2)
        return L15a          

cdef L14(double lam,double lamm,double P):
        L17a = 0.0
        cdef int i=0
        for  i from 0 <= i < np:
               L17a = L17a + LTgEh(i,lam,lam,P,2,4)
        return L17a  


cdef L21(double lam,double lamm,double P):
        L31a = 0.0
        cdef int i=0
        for  i from 0 <= i < np:
               L31a = L31a + LTgEh(i,lam,lam,P,4,1)
        return L31a

cdef L22(double lam,double lamm,double P):
        L33a = 0.0
        cdef int i=0
        for  i from 0 <= i < np:
               L33a = L33a + LTgTh(i,lam,lam,P,4,3)
        return L33a

cdef L23(double lam,double lamm,double P):
        L35a = 0.0
        cdef int i=0
        for  i from 0 <= i < np:
               L35a = L35a + LTgEh(i,lam,lam,P,4,2)
        return L35a

cdef L24(double lam,double lamm,double P):
        L37a = 0.0
        cdef int i=0
        for  i from 0 <= i < np:
               L37a = L37a + LTgTh(i,lam,lam,P,4,4)
        return L37a

cdef L31(double lam,double lamm,double P):
        L51a = 0.0
        cdef int i=0
        for  i from 0 <= i < np:
               L51a = L51a + LEgEh(i,lam,lam,P,3,1)
        return L51a

cdef L32(double lam,double lamm,double P):
        L53a = 0.0
        cdef int i=0
        for  i from 0 <= i < np:
               L53a = L53a + LTgEh(i,lam,lam,P,3,3)
        return L53a        

cdef L33(double lam,double lamm,double P):
        L55a = 0.0
        cdef int i=0
        for  i from 0 <= i < np:
               L55a = L55a + LEgEh(i,lam,lam,P,3,2)
        return L55a          

cdef L34(double lam,double lamm,double P):
        L57a = 0.0
        cdef int i=0
        for  i from 0 <= i < np:
               L57a = L57a + LTgEh(i,lam,lam,P,3,4)
        return L57a  

cdef L41(double lam,double lamm,double P):
        L71a = 0.0
        cdef int i=0
        for  i from 0 <= i < np:
               L71a = L71a + LTgEh(i,lam,lam,P,5,1)
        return L71a

cdef L42(double lam,double lamm,double P):
        L73a = 0.0
        cdef int i=0
        for  i from 0 <= i < np:
               L73a = L73a + LTgTh(i,lam,lam,P,5,3)
        return L73a

cdef L43(double lam,double lamm,double P):
        L75a = 0.0
        cdef int i=0
        for  i from 0 <= i < np:
               L75a = L75a + LTgEh(i,lam,lam,P,5,2)
        return L75a

cdef L44(double lam,double lamm,double P):
        L77a = 0.0
        cdef int i=0
        for  i from 0 <= i < np:
               L77a = L77a + LTgTh(i,lam,lam,P,5,4)
        return L77a



##############################################################################
#############################################################################

cdef  LH(double lam,double lamm,double P):
        return [[L11(lam,lamm,P),L12(lam,lamm,P),L13(lam,lamm,P),L14(lam,lamm,P)],[L21(lam,lamm,P),L22(lam,lamm,P),L23(lam,lamm,P),L24(lam,lamm,P)],[L31(lam,lamm,P),L32(lam,lamm,P),L33(lam,lamm,P),L34(lam,lamm,P)],[L41(lam,lamm,P),L42(lam,lamm,P),L43(lam,lamm,P),L44(lam,lamm,P)]]


cdef LHIN(double lam,double lamm,double P):
        return linalg.inv(LH(lam,lamm,P))









cdef c1(int s,int t,int u,int v,int l,int r,int n,double lamg,double lamh,double P,int g, int h):
    return mq[l]*mq[r]-P/4.0
cdef c2(int s,int t,int u,int v,int l,int r,int n,double lamg,double lamh,double P, int g, int h):
    return 4.0*cc*Zq[l]*Zq[r]*Zg[n]
cdef c3(int s,int t,int u,int v,int l,int r,int n,double lamg,double lamh, double P, int g, int h):
    return x4[v]*(1.0-x4[v])
cdef c4(int s,int t,int u,int v,int l,int r,int n,double lamg,double lamh, double P,int g, int h):
    return x4[v]*lamg**2+(1.0-x4[v])*(mg[n])**2
cdef c5(int s,int t,int u,int v,int l,int r,int n,double lamg,double lamh, double P,int g, int h):
    return (c2(s,t,u,v,l,r,n,lamg,lamh,P,g,h)*(x4[v])**(g-1)*(x1[s])**(h-1)*(x2[t])**(g-2)*factorial(h+g))/((4.0*pi)**2*(g-1.0)*(c3(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(g-1)*factorial(h-1)*factorial(g-2))
cdef c6(int s,int t,int u,int v,int l,int r,int n,double lamg,double lamh,double P,int g, int h):
    return (x3[u]*sqrt(P)/2.0)-(1.0-x1[s]-x2[t]-x3[u])*sqrt(P)/2.0
cdef c7(int s,int t,int u,int v,int l,int r,int n,double lamg,double lamh,double P,int g, int h):
    return x3[u]*P/4.0+x3[u]*(mq[l])**2+(1.0-x1[s]-x2[t]-x3[u])*((P/4.0)+(mq[r])**2)+x1[s]*lamh**2+x2[t]*c4(s,t,u,v,l,r,n,lamg,lamh,P,g,h)/c3(s,t,u,v,l,r,n,lamg,lamh,P,g,h)


cdef EgEh(int s,int t,int u,int v,int l,int r,int n,double lamg,double lamh,double P,int g, int h):
    return (c5(s,t,u,v,l,r,n,lamg,lamh,P,g,h)/((4.0*pi)**2))*(((c6(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**2+c1(s,t,u,v,l,r,n,lamg,lamh,P,g,h))/((h+g)*(g+h-1.0)*(c7(s,t,u,v,l,r,n,lamg,lamh,P,g,h)-(c6(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**2)**(g+h-1)))



cdef c8(int s,int t,int u,int v,int l,int r,int n,double lamg,double lamh,double P,int g, int h):
    return c7(s,t,u,v,l,r,n,lamg,lamh, P,g,h)-(c6(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**2
cdef EgTh(int s,int t,int u,int v,int l,int r,int n,double lamg,double lamh,double P,int g, int h):
    return  (c5(s,t,u,v,l,r,n,lamg,lamh,P,g,h)/((4.0*pi)**2))*((3.0*P)/(2.0*(h+g)*(g+h-1)*(g+h-2)*(g+h-3)*(c8(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(h+g-3))+(4.0*(c6(s,t,u,v,l,r,n,lamg,lamh,P,g,h)*sqrt(P))**2)/((h+g)*(h+g-1)*(h+g-2)*(c8(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(h+g-2))+(((c6(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**2+c1(s,t,u,v,l,r,n,lamg,lamh,P,g,h))*P)/(2.0*(g+h)*(g+h-1)*(g+h-2)*(c8(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(g+h-2))+(((c6(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**2+c1(s,t,u,v,l,r,n,lamg,lamh,P,g,h))*(c6(s,t,u,v,l,r,n,lamg,lamh,P,g,h)*sqrt(P))**2)/((h+g)*(g+h-1)*(c8(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(g+h-1)))


cdef c9(int s,int t,int u,int v,int l,int r,int n,double lamg,double lamh,double P,int g, int h):
    return (c2(s,t,u,v,l,r,n,lamg,lamh,P,g,h)*(x4[v])**(g-1))/(4.0*pi)**2


cdef c10(int s,int t,int u,int v,int l,int r,int n,double lamg,double lamh,double P,int g, int h):
    return (c9(s,t,u,v,l,r,n,lamg,lamh,P,g,h)*P)/(2.0*(g-1.0)*(g-2.0)*(c3(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(g-2))

cdef c11(int s,int t,int u,int v,int l,int r,int n,double lamg,double lamh,double P,int g, int h):
    return (c9(s,t,u,v,l,r,n,lamg,lamh,P,g,h)*(1.0-x4[v])**2)/((g-1.0)*(c3(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(g-1))

cdef c12(int s,int t,int u,int v,int l,int r,int n,double lamg,double lamh,double P,int g, int h):
    return (c10(s,t,u,v,l,r,n,lamg,lamh,P,g,h)*(x1[s])**(h-1)*(x2[t])**(g-3)*factorial(g+h-1))/(1.0*factorial(h-1)*factorial(g-3))

cdef c13(int s,int t,int u,int v,int l,int r,int n,double lamg,double lamh,double P,int g, int h):
    return (c11(s,t,u,v,l,r,n,lamg,lamh,P,g,h)*(x1[s])**(h-1)*(x2[t])**(g-2)*factorial(g+h))/(1.0*factorial(h-1)*factorial(g-2))


cdef c14(int s,int t,int u,int v,int l,int r,int n,double lamg,double lamh,double P,int g, int h):
    return (c6(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**2+ c1(s,t,u,v,l,r,n,lamg,lamh,P,g,h)



cdef TgTh(int s,int t,int u,int v,int l,int r,int n,double lamg,double lamh,double P,int g, int h):
    return (c12(s,t,u,v,l,r,n,lamg,lamh,P,g,h)/(4.0*pi)**2)*((3.0*P)/(2.0*(h+g-1.0)*(h+g-2.0)*(h+g-3.0)*(h+g-4.0)*(c8(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(h+g-4))+((sqrt(P)*c6(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**2*4.0)/((h+g-1.0)*(h+g-2.0)*(g+h-3.0)*(c8(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(h+g-3))+(c14(s,t,u,v,l,r,n,lamg,lamh,P,g,h)*P)/(2.0*(h+g-1.0)*(h+g-2.0)*(g+h-3.0)*(c8(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(h+g-3))+(c14(s,t,u,v,l,r,n,lamg,lamh,P,g,h)*(sqrt(P)*c6(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(2))/((g+h-1.0)*(h+g-2.0)*(c8(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(h+g-2)))+((c13(s,t,u,v,l,r,n,lamg,lamh,P,g,h))/(4.0*pi)**2)*((3.0*P**2)/((g+h)*(g+h-1.0)*(g+h-2.0)*(g+h-3.0)*(g+h-4.0)*(c8(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(g+h-4))+(9.0*P*(sqrt(P)*c6(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**2)/((h+g)*(g+h-1.0)*(g+h-2.0)*(g+h-3.0)*(c8(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(h+g-3))+(2.0*(sqrt(P)*c6(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**4)/((g+h)*(g+h-1.0)*(g+h-2.0)*(c8(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(g+h-2))+(6.0*P*(c6(s,t,u,v,l,r,n,lamg,lamh,P,g,h)*sqrt(P))**(2))/((g+h)*(g+h-1.0)*(g+h-2.0)*(g+h-3.0)*(c8(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(g+h-3))+(4.0*(sqrt(P)*c6(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(4))/((h+g)*(g+h-1.0)*(g+h-2.0)*(c8(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(g+h-2.0))+(3.0*P**2*c14(s,t,u,v,l,r,n,lamg,lamh,P,g,h))/(4.0*(g+h)*(g+h-1.0)*(g+h-2.0)*(g+h-3.0)*(c8(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(h+g-3))+(3.0*P*c14(s,t,u,v,l,r,n,lamg,lamh,P,g,h)*(sqrt(P)*c6(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(2))/((h+g)*(g+h-1.0)*(g+h-2.0)*(c8(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(h+g-2))+(c14(s,t,u,v,l,r,n,lamg,lamh,P,g,h)*(sqrt(P)*c6(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(4))/((h+g)*(h+g-1.0)*(c8(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(g+h-1.0)))   


cdef TgEh(int s,int t,int u,int v,int l,int r,int n,double lamg,double lamh,double P,int g, int h):
    return (c12(s,t,u,v,l,r,n,lamg,lamh,P,g,h)/(4.0*pi)**2)*(2.0/((g+h-1.0)*(g+h-2.0)*(g+h-3.0)*(c8(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(g+h-3))+(c14(s,t,u,v,l,r,n,lamg,lamh,P,g,h))/((g+h-1.0)*(g+h-2.0)*(c8(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(g+h-2)))+(c13(s,t,u,v,l,r,n,lamg,lamh,P,g,h)/(4.0*pi)**2)*((3.0*P)/(2.0*(h+g)*(g+h-1.0)*(g+h-2.0)*(g+h-3.0)*(c8(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(g+h-3))+(4.0*(sqrt(P)*c6(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(2))/((g+h)*(g+h-1.0)*(g+h-2.0)*(c8(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(h+g-2))+(c14(s,t,u,v,l,r,n,lamg,lamh,P,g,h)*P)/(2.0*(g+h)*(g+h-1.0)*(g+h-2.0)*(c8(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(g+h-2.0))+(c14(s,t,u,v,l,r,n,lamg,lamh,P,g,h)*(sqrt(P)*c6(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(2))/((g+h)*(g+h-1.0)*(c8(s,t,u,v,l,r,n,lamg,lamh,P,g,h))**(h+g-1)))





##############################################################################
##############################################################################

cdef R11(double lam,double lamm,double P):
        R11a = 0.0
        cdef int s = 0
        cdef int t = 0
        cdef int u = 0
        cdef int v = 0
        cdef int l = 0
        cdef int r = 0
        cdef int n = 0
        for  s from 0 <= s < np:
                for t from 0 <= t < np:
                        for u from 0 <= u <  np:
                                for v from 0 <= v <  np:
                                        for l from 0 <= l < 4:
                                                for r from 0 <= r < 4:
                                                        for  n from 0 <= n < 4:
                                                                R11a = R11a + w1[s]*w2[t]*w3[u]*w4[v]*(1.0 +e**(-2.0*100.0*(1.0 - x1[s] - x2[t] -x3[u])))**(-1)*(EgEh(s,t,u,v,l,r,n,lam,lam,P,2,1))
        return R11a



cdef R12(double lam,double lamm,double P):
        R13a = 0.0
        cdef int s = 0
        cdef int t = 0
        cdef int u = 0
        cdef int v = 0
        cdef int l = 0
        cdef int r = 0
        cdef int n = 0
        for  s from 0 <= s < np:
                for t from 0 <= t < np:
                        for u from 0 <= u <  np:
                                for v from 0 <= v <  np:
                                        for l from 0 <= l < 4:
                                                for r from 0 <= r < 4:
                                                        for  n from 0 <= n < 4:
                                                                R13a = R13a + w1[s]*w2[t]*w3[u]*w4[v]*(1.0 +e**(-2.0*100.0*(1.0 - x1[s] - x2[t] -x3[u])))**(-1)*(EgTh(s,t,u,v,l,r,n,lam,lam,P,2,3))                                             
        return R13a



cdef R13(double lam,double lamm,double P):
        R15a = 0.0
        cdef int s = 0
        cdef int t = 0
        cdef int u = 0
        cdef int v = 0
        cdef int l = 0
        cdef int r = 0
        cdef int n = 0
        for  s from 0 <= s < np:
                for t from 0 <= t < np:
                        for u from 0 <= u <  np:
                                for v from 0 <= v <  np:
                                        for l from 0 <= l < 4:
                                                for r from 0 <= r < 4:
                                                        for  n from 0 <= n < 4:
                                                                R15a = R15a + w1[s]*w2[t]*w3[u]*w4[v]*(1.0 +e**(-2.0*100.0*(1.0 - x1[s] - x2[t] -x3[u])))**(-1)*(EgEh(s,t,u,v,l,r,n,lam,lam,P,2,2))
        return R15a



cdef R14(double lam,double lamm,double P):
        R17a = 0.0
        cdef int s = 0
        cdef int t = 0
        cdef int u = 0
        cdef int v = 0
        cdef int l = 0
        cdef int r = 0
        cdef int n = 0
        for  s from 0 <= s < np:
                for t from 0 <= t < np:
                        for u from 0 <= u <  np:
                                for v from 0 <= v <  np:
                                        for l from 0 <= l < 4:
                                                for r from 0 <= r < 4:
                                                        for  n from 0 <= n < 4:
                                                                R17a = R17a + w1[s]*w2[t]*w3[u]*w4[v]*(1.0 +e**(-2.0*100.0*(1.0 - x1[s] - x2[t] -x3[u])))**(-1)*(EgTh(s,t,u,v,l,r,n,lam,lam,P,2,4))
        return R17a




#############################################################################
############################################################################







                          
cdef R21(double lam,double lamm, double P):
        R31a = 0.0
        cdef int s = 0
        cdef int t = 0
        cdef int u = 0
        cdef int v = 0
        cdef int l = 0
        cdef int r = 0
        cdef int n = 0
        for  s from 0 <= s < np:
                for t from 0 <= t < np:
                        for u from 0 <= u <  np:
                                for v from 0 <= v <  np:
                                        for l from 0 <= l < 4:
                                                for r from 0 <= r < 4:
                                                        for  n from 0 <= n < 4: 
                                                                R31a = R31a + w1[s]*w2[t]*w3[u]*w4[v]*(1.0 +e**(-2.0*100.0*(1.0 - x1[s] - x2[t] -x3[u])))**(-1)*(TgEh(s,t,u,v,l,r,n,lam,lam,P,4,1))                                     
        return R31a





 


cdef R22(double lam,double lamm, double P):
        R33a = 0.0
        cdef int s = 0
        cdef int t = 0
        cdef int u = 0
        cdef int v = 0
        cdef int l = 0
        cdef int r = 0
        cdef int n = 0
        for  s from 0 <= s < np:
                for t from 0 <= t < np:
                        for u from 0 <= u <  np:
                                for v from 0 <= v <  np:
                                        for l from 0 <= l < 4:
                                                for r from 0 <= r < 4:
                                                        for  n from 0 <= n < 4: 
                                                                R33a = R33a + w1[s]*w2[t]*w3[u]*w4[v]*(1.0 +e**(-2.0*100.0*(1.0 - x1[s] - x2[t] -x3[u])))**(-1)*(TgTh(s,t,u,v,l,r,n,lam,lam,P,4,3))                                     
        return R33a








cdef R23(double lam,double lamm, double P):
        R35a = 0.0
        cdef int s = 0
        cdef int t = 0
        cdef int u = 0
        cdef int v = 0
        cdef int l = 0
        cdef int r = 0
        cdef int n = 0
        for  s from 0 <= s < np:
                for t from 0 <= t < np:
                        for u from 0 <= u <  np:
                                for v from 0 <= v <  np:
                                        for l from 0 <= l < 4:
                                                for r from 0 <= r < 4:
                                                        for  n from 0 <= n < 4: 
                                                                R35a = R35a + w1[s]*w2[t]*w3[u]*w4[v]*(1.0 +e**(-2.0*100.0*(1.0 - x1[s] - x2[t] -x3[u])))**(-1)*(TgEh(s,t,u,v,l,r,n,lam,lam,P,4,2))                                     
        return R35a






cdef R24(double lam,double lamm, double P):
        R37a = 0.0
        cdef int s = 0
        cdef int t = 0
        cdef int u = 0
        cdef int v = 0
        cdef int l = 0
        cdef int r = 0
        cdef int n = 0
        for  s from 0 <= s < np:
                for t from 0 <= t < np:
                        for u from 0 <= u <  np:
                                for v from 0 <= v <  np:
                                        for l from 0 <= l < 4:
                                                for r from 0 <= r < 4:
                                                        for  n from 0 <= n < 4: 
                                                                R37a = R37a + w1[s]*w2[t]*w3[u]*w4[v]*(1.0 +e**(-2.0*100.0*(1.0 - x1[s] - x2[t] -x3[u])))**(-1)*(TgTh(s,t,u,v,l,r,n,lam,lam,P,4,4))                                     
        return R37a









##############################################################################
##############################################################################







cdef R31(double lam,double lamm, double P):
        R51a = 0.0
        cdef int s = 0
        cdef int t = 0
        cdef int u = 0
        cdef int v = 0
        cdef int l = 0
        cdef int r = 0
        cdef int n = 0
        for  s from 0 <= s < np:
                for t from 0 <= t < np:
                        for u from 0 <= u <  np:
                                for v from 0 <= v <  np:
                                        for l from 0 <= l < 4:
                                                for r from 0 <= r < 4:
                                                        for  n from 0 <= n < 4: 
                                                                R51a = R51a + w1[s]*w2[t]*w3[u]*w4[v]*(1.0 +e**(-2.0*100.0*(1.0 - x1[s] - x2[t] -x3[u])))**(-1)*(EgEh(s,t,u,v,l,r,n,lam,lam,P,3,1))                                     
        return R51a






cdef R32(double lam,double lamm, double P):
        R53a = 0.0
        cdef int s = 0
        cdef int t = 0
        cdef int u = 0
        cdef int v = 0
        cdef int l = 0
        cdef int r = 0
        cdef int n = 0
        for  s from 0 <= s < np:
                for t from 0 <= t < np:
                        for u from 0 <= u <  np:
                                for v from 0 <= v <  np:
                                        for l from 0 <= l < 4:
                                                for r from 0 <= r < 4:
                                                        for  n from 0 <= n < 4: 
                                                                R53a = R53a + w1[s]*w2[t]*w3[u]*w4[v]*(1.0 +e**(-2.0*100.0*(1.0 - x1[s] - x2[t] -x3[u])))**(-1)*(EgTh(s,t,u,v,l,r,n,lam,lam,P,3,3))                                     
        return R53a









cdef R33(double lam,double lamm, double P):
        R55a = 0.0
        cdef int s = 0
        cdef int t = 0
        cdef int u = 0
        cdef int v = 0
        cdef int l = 0
        cdef int r = 0
        cdef int n = 0
        for  s from 0 <= s < np:
                for t from 0 <= t < np:
                        for u from 0 <= u <  np:
                                for v from 0 <= v <  np:
                                        for l from 0 <= l < 4:
                                                for r from 0 <= r < 4:
                                                        for  n from 0 <= n < 4: 
                                                                R55a = R55a + w1[s]*w2[t]*w3[u]*w4[v]*(1.0 +e**(-2.0*100.0*(1.0 - x1[s] - x2[t] -x3[u])))**(-1)*(EgEh(s,t,u,v,l,r,n,lam,lam,P,3,2))                                     
        return R55a  


 



cdef R34(double lam,double lamm, double P):
        R57a = 0.0
        cdef int s = 0
        cdef int t = 0
        cdef int u = 0
        cdef int v = 0
        cdef int l = 0
        cdef int r = 0
        cdef int n = 0
        for  s from 0 <= s < np:
                for t from 0 <= t < np:
                        for u from 0 <= u <  np:
                                for v from 0 <= v <  np:
                                        for l from 0 <= l < 4:
                                                for r from 0 <= r < 4:
                                                        for  n from 0 <= n < 4: 
                                                                R57a = R57a + w1[s]*w2[t]*w3[u]*w4[v]*(1.0 +e**(-2.0*100.0*(1.0 - x1[s] - x2[t] -x3[u])))**(-1)*(EgTh(s,t,u,v,l,r,n,lam,lam,P,3,4))                                     
        return R57a  







################################################################################
###############################################################################



                          
cdef R41(double lam,double lamm, double P):
        R71a = 0.0
        cdef int s = 0
        cdef int t = 0
        cdef int u = 0
        cdef int v = 0
        cdef int l = 0
        cdef int r = 0
        cdef int n = 0
        for  s from 0 <= s < np:
                for t from 0 <= t < np:
                        for u from 0 <= u <  np:
                                for v from 0 <= v <  np:
                                        for l from 0 <= l < 4:
                                                for r from 0 <= r < 4:
                                                        for  n from 0 <= n < 4: 
                                                                R71a = R71a + w1[s]*w2[t]*w3[u]*w4[v]*(1.0 +e**(-2.0*100.0*(1.0 - x1[s] - x2[t] -x3[u])))**(-1)*(TgEh(s,t,u,v,l,r,n,lam,lam,P,5,1))                                     
        return R71a





cdef R42(double lam,double lamm, double P):
        R73a = 0.0
        cdef int s = 0
        cdef int t = 0
        cdef int u = 0
        cdef int v = 0
        cdef int l = 0
        cdef int r = 0
        cdef int n = 0
        for  s from 0 <= s < np:
                for t from 0 <= t < np:
                        for u from 0 <= u <  np:
                                for v from 0 <= v <  np:
                                        for l from 0 <= l < 4:
                                                for r from 0 <= r < 4:
                                                        for  n from 0 <= n < 4: 
                                                                R73a = R73a + w1[s]*w2[t]*w3[u]*w4[v]*(1.0 +e**(-2.0*100.0*(1.0 - x1[s] - x2[t] -x3[u])))**(-1)*(TgTh(s,t,u,v,l,r,n,lam,lam,P,5,3))                                     
        return R73a










cdef R43(double lam,double lamm, double P):
        R75a = 0.0
        cdef int s = 0
        cdef int t = 0
        cdef int u = 0
        cdef int v = 0
        cdef int l = 0
        cdef int r = 0
        cdef int n = 0
        for  s from 0 <= s < np:
                for t from 0 <= t < np:
                        for u from 0 <= u <  np:
                                for v from 0 <= v <  np:
                                        for l from 0 <= l < 4:
                                                for r from 0 <= r < 4:
                                                        for  n from 0 <= n < 4: 
                                                                R75a = R75a + w1[s]*w2[t]*w3[u]*w4[v]*(1.0 +e**(-2.0*100.0*(1.0 - x1[s] - x2[t] -x3[u])))**(-1)*(TgEh(s,t,u,v,l,r,n,lam,lam,P,5,2))                                     
        return R75a






       


cdef R44(double lam,double lamm, double P):
        R77a = 0.0
        cdef int s = 0
        cdef int t = 0
        cdef int u = 0
        cdef int v = 0
        cdef int l = 0
        cdef int r = 0
        cdef int n = 0
        for  s from 0 <= s < np:
                for t from 0 <= t < np:
                        for u from 0 <= u <  np:
                                for v from 0 <= v <  np:
                                        for l from 0 <= l < 4:
                                                for r from 0 <= r < 4:
                                                        for  n from 0 <= n < 4: 
                                                                R77a = R77a + w1[s]*w2[t]*w3[u]*w4[v]*(1.0 +e**(-2.0*100.0*(1.0 - x1[s] - x2[t] -x3[u])))**(-1)*(TgTh(s,t,u,v,l,r,n,lam,lam,P,5,4))                                     
        return R77a







#################################################################################
###############################################################################


                          

##############################################################################
##############################################################################
###############################################################################
#############################################################################
###############################################################################
###############################################################################









comm = MPI.COMM_WORLD
size = comm.size
rank = comm.rank
name = MPI.Get_processor_name()



cdef RH1(double lam,double lamm,double P):
    return [R11(lam,lamm,P),R12(lam,lamm,P),R13(lam,lamm,P),R14(lam,lamm,P),R21(lam,lamm,P),R22(lam,lamm,P),R23(lam,lamm,P),R24(lam,lamm,P),R31(lam,lamm,P),R32(lam,lamm,P),R33(lam,lamm,P),R34(lam,lamm,P),R41(lam,lamm,P),R42(lam,lamm,P),R43(lam,lamm,P),R44(lam,lamm,P)]
cdef RH2(double lam,double lamm,double P,int rank):
    return RH1(lam,lamm,P)[rank]

 
aa = []
P = -0.035
#P1 = linspace(-0.04,0,num=20)
lam1 = linspace(0.05,2.5,50)
#lam1 = 0.6
lamm1 = 0.65
for i in range(50):
    if rank == 0:
        lam = lam1[i]
        #lam = lam1
        #P = P1[i]
        R1z = RH2(lam,lamm1,P,rank)
        comm.send(R1z,dest = 17)
    elif rank == 1:
        lam = lam1[i]
        #lam = lam1
        #P = P1[i]
        R2z = RH2(lam,lamm1,P,rank)
        comm.send(R2z,dest = 17)
    elif rank == 2:
        lam = lam1[i]
        #lam = lam1
        #P = P1[i]
        R3z = RH2(lam,lamm1,P,rank)
        comm.send(R3z,dest = 17)
    elif rank == 3:
        lam = lam1[i]
        #lam = lam1
        #P = P1[i]
        R4z = RH2(lam,lamm1,P,rank)
        comm.send(R4z,dest = 17)
    elif rank == 4:
        lam = lam1[i]
        #lam = lam1
        #P = P1[i]
        R5z = RH2(lam,lamm1,P,rank)
        comm.send(R5z,dest = 17)
    elif rank == 5:
        lam = lam1[i]
        #lam = lam1
        #P = P1[i]
        R6z = RH2(lam,lamm1,P,rank)
        comm.send(R6z,dest = 17)
    elif rank == 6:
        lam = lam1[i]
        #lam = lam1
        #P = P1[i]
        R7z = RH2(lam,lamm1,P,rank)
        comm.send(R7z,dest = 17)
    elif rank == 7:
        lam = lam1[i]
        #lam = lam1
        #P = P1[i]
        R8z = RH2(lam,lamm1,P,rank)
        comm.send(R8z,dest = 17)
    elif rank == 8:
        lam = lam1[i]
        #lam = lam1
        #P = P1[i]
        R9z = RH2(lam,lamm1,P,rank)
        comm.send(R9z,dest = 17)
    elif rank == 9:
        lam = lam1[i]
        #lam = lam1
        #P = P1[i]
        R10z = RH2(lam,lamm1,P,rank)
        comm.send(R10z,dest = 17)
    elif rank == 10:
        lam = lam1[i]
        #lam = lam1
        #P = P1[i]
        R11z = RH2(lam,lamm1,P,rank)
        comm.send(R11z,dest = 17)
    elif rank == 11:
        lam = lam1[i]
        #lam = lam1
        #P = P1[i]
        R12z = RH2(lam,lamm1,P,rank)
        comm.send(R12z,dest = 17)
    elif rank == 12:
        lam = lam1[i]
        #lam = lam1
        #P = P1[i]
        R13z = RH2(lam,lamm1,P,rank)
        comm.send(R13z,dest = 17)
    elif rank == 13:
        lam = lam1[i]
        #lam = lam1
        #P = P1[i]
        R14z = RH2(lam,lamm1,P,rank)
        comm.send(R14z,dest = 17)
    elif rank == 14:
        lam = lam1[i]
        #lam = lam1
        #P = P1[i]
        R15z = RH2(lam,lamm1,P,rank)
        comm.send(R15z,dest = 17)
    elif rank == 15:
        lam = lam1[i]
        #lam = lam1
        #P = P1[i]
        R16z = RH2(lam,lamm1,P,rank)
        comm.send(R16z,dest = 17)   
    elif rank == 16:
        lam = lam1[i]
        #lam = lam1
        #P = P1[i]
        LIN = LHIN(lam,lamm1,P)
        comm.send(LIN, dest = 17)
    elif rank == 17:
        bbb = [0.0 for k in range(16)]
        for r in range(16):
            bbb[r] = comm.recv(source = r)
        bbbb = mat(bbb)
        RH = bbbb.reshape(4,4)
        LIN = comm.recv(source = 16)
        mat11 = linalg.eigvals(dot(LIN,RH))
        mat1 = [0.0 for u in range(4)]
        for l in range(4):
            mat1[l] = mat11[l].real
        print "{"+repr(lam1[i])+",{"+repr(mat1[0])+","+repr(mat1[1])+","+repr(mat1[2])+","+repr(mat1[3])+"}}"
        #print "{"+repr(P1[i])+",{"+repr(mat1[0])+","+repr(mat1[1])+","+repr(mat1[2])+","+repr(mat1[3])+"}}"
        fin = max(mat1)
        #print "At lambda value " + repr(lam1[i]) + " the maximum eigenvalue is " + repr(fin) + "."
        comm.send(fin,dest = 18)
    elif rank == 18:
        finn = comm.recv(source = 17)
        aa.append(finn)
        if len(aa) == 50:
            bfun = [0.0 for k in range(50)]
            for v in range(50):
                bfun[v] = - aa[v]
            fnct = interp1d(lam1,bfun,bounds_error=False)
            lammax = optimize.fmin(fnct,1.0)
            print "The maximum is reached when lambda = " + repr(lammax)
        else:
            print "Working..." 








#def mat(lam,P):
 #   return linalg.eigvals(dot(LHIN(lam,P),RH(lam,P)))



"""def mat1(lam,P):
    mat2 = [0 for j in range(4)]
    for j in range(4):
        mat2[j] = mat(lam,P)[j].real
    return mat2"""


#def fin(lam,P):
 #   return max(mat1(lam,P))




#def fin1(lam):  
 #   if lam > 0.0:
  #      return -fin(lam,-0.025)
   # else:
    #    return 0

"""test2 = lambda x: fin1(x[0])
Nfeval = 1
guess = [0.4]
maxiter1 = 50
xmin = [0.05]
xmax = [50]
bounds = [(low,high) for low,high in zip(xmin,xmax)]
def callbackF(Xi):
    global Nfeval
    print '{0:4d}   {1: 3.6f}   {2: 3.6f}'.format(Nfeval, Xi[0], test2(Xi))
    Nfeval += 1


#minimizer_kwargs = dict(method='L-BFGS-B', bounds=bounds,callback=callbackF)

def print_fun(x,f,accepted):
        print("at minimum %.4f at %.4f accepted %d" %(f,x,int(accepted)))

class RandomDisplacementBounds(object):
       
        def __init__(self,xmin,xmax,stepsize=0.2):
                self.xmin = xmin
                self.xmax = xmax
                self.stepsize = stepsize
        def __call__(self,x):
                
                while True:
                        xnew = x + random.uniform(-self.stepsize,self.stepsize,shape(x))
                        if all(xnew < xmax) and all(xnew > xmin):
                                break
                return xnew"""

#take_step = RandomDisplacementBounds(xmin,xmax)
#result = optimize.basinhopping(test2,guess,niter = 1000,callback=print_fun,minimizer_kwargs=minimizer_kwargs)#,take_step=take_step)
#print result
#print optimize.fmin(test2,guess,callback=callbackF)
