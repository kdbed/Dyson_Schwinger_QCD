
# correct solution for DSE using Qin-Chang interation - u/d quark
# Calculates A,B, sigs, sigv, Z2
# Uses 128 Gauss pts imported from mathematica

#####################################################################
# rainbow-ladder approximation  
#####################################################################

from numpy import *
import matplotlib.pyplot as plt





#parameters
# number of flavors
f = 4.0
# parameters for the Qin-Chang model (Gev)
w = 0.5
#D = 0.781
D = pow(0.82,3)/w
# renormalization point (GeV)
mu = 2.0
#  UV cutoff (Gev)
Lam = 100.0*mu
# Lambda QCD (Gev)
lamQ = 0.234
# IR cutoff (Gev)
lam = pow(10.0,-2)
# renormalized mass  (GeV)
m = 0.0178

mt = 0.5
#product 
pr = Lam*lam
# quotient
qu = Lam/lam
# order in Legendre polys
N = 200
M = 200
index = 104  # index 104 of 200 corresponds to p^2 = 4.01829,  68 of 131 (gausspts2.dat) gives 4.06
gamma = (12.0)/(33.0-2.0*f)
tau = pow(e,2) - 1.0






x = loadtxt('gausspts1.dat')
weights = loadtxt('gausswts1.dat')

#Import gauss pts
#x,weights = gaussxw(N)
# the p-values are determined by the zeros of the Legendre functions
p = [(pr)*(qu)**(x[i]) for i in range(N)]
k = [(pr)*(qu)**(x[i]) for i in range(N)]          # defines integration values















# define initial guess for A,B
B = [m for i in range(N)]
A = [1.0 for i in range(N)]








# Define kernel (qin-chang interaction. paper: interaction model for gap equation)
def F(var1):
    return (1.0-e**(-var1/(4.0*pow(mt,2))))/(var1)
def alpha(var2):
    return  (2.0/w**4)*D*e**(-var2/w**2) +(2.0*gamma*F(var2))/(log(tau+(1.0+(var2/lamQ**2))**2))

# for the z-int
#u,wts = gaussxw(M)
u = loadtxt('gausspts1.dat')
wts = loadtxt('gausswts1.dat')
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

print zbint[0][0]
#
for j in range(N):
    for i in range(N):
        for l in range(M):
            zaint[i][j] += wts[l]*alpha(q[l][i][j])*sqrt(1-u[l]**2)*(sqrt(p[i])*sqrt(k[j])*u[l]+2*(p[i]-sqrt(p[i])*sqrt(k[j])*u[l])*(sqrt(p[i])*sqrt(k[j])*u[l]-k[j])/(q[l][i][j]))
print zaint[0][0]
#
#
#
#
G = [0.0 for i in range(N)]
H = [0.0 for i in range(N)]
#
#
z2 = 1.0
m0 = m



#
#
# integrations: loop in i is for p-values, loop in j does integration for given p
epsilon = 0
delta = 0
eta = 0
while epsilon == 0 or delta == 0:
#while epsilon < 100:
    for i in range(N):           # for given p-value
        s = 0.0
        r = 0.0
        for j in range(N):       # for given k-value
             s += weights[j]*k[j]**2*log(qu)*zaint[i][j]*(A[j]/((k[j])*(A[j])**2+(B[j])**2))
             r += weights[j]*k[j]**2*log(qu)*zbint[i][j]*(B[j]/((k[j])*(A[j])**2+(B[j])**2))
        G[i] = (2.0/(3.0*(pi)))*s*(1.0/(p[i]))*(z2)**2
        H[i] = (2.0/(pi))*r*(z2)**2
    print G[0]
    print H[0]  
    oldm = m0
    oldz = z2
    z2 = 1.0-G[index]
    m0 = (m-H[index])/(z2)
    T = [A[i] for i in range(N)]
    R = [B[j] for j in range(N)]  
    # match the boundary conditions for new A,B
    A = [z2+G[i] for i in range(N)]
    B = [z2*m0+H[j] for j in range(N)]
    print A
    print B
    delta = allclose(A,T,1e-15,0)
    epsilon = allclose(B,R,1e-15,0)
   # epsilon += 1

print A[index]
print B[index]
print 'z2 = ', z2
print 'm0 = ', m0
z4 = (z2*m0)/(m)
print 'z4 = ', z4

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


mhat2 = []
for i in range(N):
    mhat2.append(2)
for i in range(N):
    mhat2[i] = (I[i]-(2.0*pow(pi,2)*gamma*0.012514)/(3.0*O[i]*pow(0.5*log(O[i]/((lamQ)**2)),(1.0-gamma))))*pow(0.5*log(O[i]/((lamQ)**2)),gamma)
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

mhat1 = open('mhat.dat','w')
for i, j in zip(O,mhat2):
    mhat1.write(str(i) + "   " + str(j) + '\n')
mhat1.close()

m_q = open('m_q.dat','w')
for i, j in zip(O,I):
    m_q.write(str(i) + "   " + str(j) + '\n')
m_q.close()


a_p_Q = open('a_p_Q.dat','w')
for i,j in zip(O,A):
    a_p_Q.write(str(i) + "   " + str(j) + '\n')
a_p_Q.close()
b_p_Q = open('b_p_Q.dat', 'w')
for i, j in zip(O,B):
    b_p_Q.write(str(i) + "   " + str(j) + '\n')
b_p_Q.close()
sig_v_Q = open('sig_v_Q.dat','w')
for i,j in zip(O,sigv):
    sig_v_Q.write(str(i) + "   " + str(j) + '\n')
sig_v_Q.close()
sig_s_Q = open('sig_s_Q.dat', 'w')
for i,j in zip(O,sigs):
    sig_s_Q.write(str(i) + "    " + str(j) + '\n')
sig_s_Q.close()
