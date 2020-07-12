from sympy import KroneckerDelta

from numpy import *
from cmath import *
import scipy.optimize as optimize

siv = loadtxt('sig_v_Q2.dat')
sis = loadtxt('sig_s_Q2.dat')

dat1 = [[1,siv[i][0]] for i in range(len(siv))]
dat2 = [[2,sis[i][0]] for i in range(len(sis))]
for i in range(len(siv)):
    dat1.append(dat2[i])
xdat = array(dat1)

ydat1 = [siv[i][1] for i in range(len(siv))]
ydat2 = [sis[i][1] for i in range(len(siv))]
for i in range(len(siv)):
    ydat1.append(ydat2[i])

ydat = array(ydat1)


def func(x,y,z1,z2,mr,m,z3,z4,mr1,m1):
    return KroneckerDelta(x,1)*((complex(z1,z2))/(y + (complex(mr,m))**2) + (complex(z1, - z2))/(y + (complex(mr, - m))**2) + (complex(z3, z4))/(y + (complex(mr1, m1))**2) + (complex(z3, - z4))/(y + (complex(mr1, - m1))**2)) + KroneckerDelta(x,2)*(((complex(z1, z2))*(complex(mr, m)))/(y + (complex(mr, m))**2) + ((complex(z1, - z2))*(complex(mr, - m)))/(y + (complex(mr, - m))**2) + ((complex(z3, z4))*(complex(mr1, m1)))/(y + (complex(mr1,  m1))**2) + ((complex(z3, - z4))*(complex(mr1, - m1)))/(y + (complex(mr1, - m1))**2))
    





def minfunc(z):
    sum = 0.0
    for i in range(2*len(siv)):
        sum += (xdat[i][0])*float((func(xdat[i][0],xdat[i][1],z[0],z[1],z[2],z[3],z[4],z[5],z[6],z[7])-ydat[i])**2/(ydat[i]))
    return sum













p0 = [1.1,-2.1,.1,0.1,0.1,-0.8,1,0.0]
minimizer_kwargs = {"method":"BFGS"}

def print_fun(x, f, accepted):
         print("at  z1 -> %.5f , z2 -> %.5f, mr -> %.5f , m -> %.5f , z3 -> %.5f , z4 -> %.5f , mr1 -> %.5f , m1 -> %.5f  function is %.5f accepted %d" % (x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7], f, int(accepted)))








class MyBounds(object):
     def __init__(self, xmax=[15,15,15,15,15,15,15,15], xmin=[-15,-15,-15,-15,-15,-15,-15,-15] ):
         self.xmax = array(xmax)
         self.xmin = array(xmin)
     def __call__(self, **kwargs):
         x = kwargs["x_new"]
         tmax = bool(all(x <= self.xmax))
         tmin = bool(all(x >= self.xmin))
         return tmax and tmin


mybounds = MyBounds()


ret =  optimize.basinhopping(minfunc,p0,minimizer_kwargs=minimizer_kwargs,callback=print_fun, niter = 100, accept_test=mybounds) ## ,take_step=take_step)
print ret




'''
a = 0.25


rranges = ((-4, 4),(-4, 4), (-4, 4), (-4, 4), (-4, 4), (-4, 4), (-4, 4), (-4, 4))
ret = optimize.brute(minfunc,rranges,full_output = True,finish = optimize.fmin)
print ret[0]
print ret[1]


'''
