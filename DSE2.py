# this is the solution at 19gev; the points from mathematica need updated 
# Correct solution to Maris-Tandy model; prints .dat file for A,B,sigs, sigv
# Uses 128 Gauss pts imported from mathematica

#####################################################################
#####################################################################

from numpy import *
import matplotlib.pyplot as plt
#from gausswts import gaussxw




#parameters
# number of flavors
f = 4.0
# parameters for the Maris-Tandy model (Gev)
w = 0.4
#D = 0.781
D = 0.93
# renormalization point (GeV)
mu = 19.0
#  UV cutoff (Gev)
Lam = 100.0*mu
# Lambda QCD (Gev)
lamQ = 0.234
# IR cutoff (Gev)
lam = pow(10.0,-2)
# renormalized mass  (GeV)
m = 0.00374
mt = 0.5
#product 
pr = Lam*lam
# quotient
qu = Lam/lam
# order in Legendre polys
N = 128
M = 128
index = 73
gamma = (12.0)/(33.0-2.0*f)
tau = pow(e,2) - 1.0



# 200-84  , 300 - 126, 128-53/4 



x = loadtxt('gausspts.dat')
weights = loadtxt('gausswts.dat')

#Import gauss pts
#x,weights = gaussxw(N)
# the p-values are determined by the zeros of the Legendre functions
p = [(pr)*(qu)**(x[i]) for i in range(N)]
k = [(pr)*(qu)**(x[i]) for i in range(N)]          # defines integration values















# define intial guess for A,B
B = [m for i in range(N)]
A = [1 for i in range(N)]








# Define kernel (Maris-Tandy model (k^2 in first part of alpha, type var2 before D) or qin-chang (no k^2)
def F(var1):
    return (1.0-e**(-var1/(4.0*pow(mt,2))))/(var1)
def alpha(var2):
    return  (1.0/w**6)*var2*D*e**(-var2/w**2) +(2.0*gamma*F(var2))/(log(tau+(1.0+(var2/lamQ**2))**2))

# for the z-int
#v,wts2 = gaussxw(M)
u = loadtxt('gausspts.dat')
wts = loadtxt('gausswts.dat')
q = [[[0.0 for i in range(N)] for j in range(N)] for l in range(M)]
zaint = [[0.0 for i in range(N)] for j in range(N)]
zbint = [[0.0 for i in range(N)] for j in range(N)]
# do the angular integral
for i in range(M):
    for j in range(N):
        for l in range(N):
            q[i][j][l] = p[j]+k[l]-2*sqrt(k[l])*sqrt(p[j])*u[i]
#
#
#
#
for j in range(N):
    for i in range(N):
        for l in range(M):
            zbint[i][j] += wts[l]*alpha(q[l][i][j])*sqrt(1-u[l]**2)
#
for j in range(N):
    for i in range(N):
        for l in range(M):
            zaint[i][j] += wts[l]*alpha(q[l][i][j])*sqrt(1-u[l]**2)*(sqrt(p[i])*sqrt(k[j])*u[l]+2*(p[i]-sqrt(p[i])*sqrt(k[j])*u[l])*(sqrt(p[i])*sqrt(k[j])*u[l]-k[j])/(q[l][i][j]))
#
#
#
#
G = [0.0 for i in range(N)]
H = [0.0 for i in range(N)]
#
#
#
#
# integrations: loop in i is for p-values, loop in j does integration for given p
epsilon = 0
delta = 0
while epsilon == 0 or delta == 0:
#while epsilon < 1000:
    for i in range(N):           # for given p-value
        s = 0.0
        r = 0.0
        for j in range(N):       # for given k-value
             s += weights[j]*k[j]**2*log(qu)*zaint[i][j]*(A[j]/((k[j])*(A[j])**2+(B[j])**2))
             r += weights[j]*k[j]**2*log(qu)*zbint[i][j]*(B[j]/((k[j])*(A[j])**2+(B[j])**2))
        G[i] = (2.0/(3.0*(pi)))*s*(1.0/(p[i]))
        H[i] = (2.0/(pi))*r
#
    T = [A[i] for i in range(N)]
    R = [B[j] for j in range(N)]  
    # match the boundary conditions for new A,B
    A = [G[i]-G[index]+1.0 for i in range(N)]
    B = [m+H[j]-H[index] for j in range(N)]
    delta = allclose(A,T,1e-15,0)
    epsilon = allclose(B,R,1e-15,0)
   # epsilon += 1

print A[index]
print B[index]

I = []
for i in range(N):
    I.append(i)

for i in range(N):
    I[i] =(B[i])/A[i]
O = []
for i in range(N):
    O.append(i)
for i in range(N):
    O[i] = p[i]

sigv = [A[i]/(O[i]*pow(A[i],2) + pow(B[i],2)) for i in range(N)]
sigs = [B[i]/(O[i]*pow(A[i],2) + pow(B[i],2)) for i in range(N)]


plt.loglog(O,A,'r')
plt.xlabel('p^2 (GeV)')


plt.loglog(O,B,'b')
plt.xlim(10**(-3),10**3)
plt.ylim(10**(-3),3)
plt.show()

plt.loglog(O,I)
plt.xlabel('p (GeV)')
plt.ylabel('M(p)')
plt.xlim(10**(-3),10**3)
plt.ylim(0,10)
plt.show()



m_mt = open('m_mt.dat','w')
for i, j in zip(O,I):
    m_mt.write(str(i) + "   " + str(j) + '\n')
m_mt.close()
a_p = open('a_p.dat','w')
for i,j in zip(O,A):
    a_p.write(str(i) + "   " + str(j) + '\n')
a_p.close()
b_p = open('b_p.dat', 'w')
for i, j in zip(O,B):
    b_p.write(str(i) + "   " + str(j) + '\n')
b_p.close()
sig_v = open('sig_v.dat','w')
for i,j in zip(O,sigv):
    sig_v.write(str(i) + "   " + str(j) + '\n')
sig_v.close()
sig_s = open('sig_s.dat', 'w')
for i,j in zip(O,sigs):
    sig_s.write(str(i) + "    " + str(j) + '\n')
sig_s.close()
