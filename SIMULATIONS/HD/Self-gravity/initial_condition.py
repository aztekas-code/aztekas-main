#! /usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy import integrate
from scipy import interpolate

# define functions used in integrate
def f(x, y, z):
    return -z/x**2

def g(x, y, z, n, c):
    return c * x**2 * y**n

# runge-kutta method
def rk(x_0, y_0, z_0, n, c, t_step):
    '''Integrates one time step'''
    k_0 = t_step * f(x_0, y_0, z_0)
    l_0 = t_step * g(x_0, y_0, z_0, n, c)
    k_1 = t_step * f(x_0+1/2*t_step, y_0+1/2*k_0, z_0+1/2*l_0)
    l_1 = t_step * g(x_0+1/2*t_step, y_0+1/2*k_0, z_0+1/2*l_0, n, c)
    k_2 = t_step * f(x_0+1/2*t_step, y_0+1/2*k_1, z_0+1/2*l_1)
    l_2 = t_step * g(x_0+1/2*t_step, y_0+1/2*k_1, z_0+1/2*l_1, n, c)
    k_3 = t_step * f(x_0+t_step, y_0+k_2, z_0+l_2)
    l_3 = t_step * g(x_0+t_step, y_0+k_2, z_0+l_2, n, c)
    x_1 = x_0 + t_step
    y_1 = y_0 + 1/6 * (k_0+2*k_1+2*k_2+k_3)
    z_1 = z_0 + 1/6 * (l_0+2*l_1+2*l_2+l_3)
    return (x_1, y_1, z_1)

def rk_integrate(n, rho0, h0, t_step):
    # define initial conditions
    r_0 = 1.0   # r/R
    h_0 = h0   # h/(GM/R)
    m_0 = 1.0   # m/M
    
    #c = ( rho0/(M/4piR**3) ) *( (GM/R)/h )**n
    c = rho0/(h0**n)

    # do the integration and compile the results into two lists
    r = [r_0]
    h = [h_0]
    m = [m_0]
    while h_0>0:
        r_0, h_0, m_0 = rk(r_0, h_0, m_0, n, c, t_step)
        if type(h_0) == complex: # complex values of h_0 sometimes occur when n is not an integer
            break
        r.append(r_0)
        h.append(h_0)
        m.append(m_0)
    return (r, h, m)

def thermo(rho0,h0,h,gamma):
    n = 1./(gamma-1)
    P0 = h0*rho0*(gamma-1.)/gamma
 
    rho = rho0*(h/h0)**n
    P = P0*(rho/rho0)**gamma
    
    return rho, P

# define integration variables
gamma = 5./3.
n = 1./(gamma-1)
t_step = 0.01

rho0 = 0.1
h0 = 10.0

output='initial.dat'

r, h, m = rk_integrate(n, rho0, h0, t_step)

r = np.array(r)
h = np.array(h)

print ("\n")
print ("cloud mass ",m[-1]-1,"radius = ",r[-1])

interp2 = interpolate.interp1d(r, h, kind = "quadratic")

# aztekas log grid
gc = 3
x1min = 1.0
x1max = 0.9999*r[-1]
Nx1 = 100 + 2*gc

rad = np.zeros(Nx1+1)

for i in range(Nx1+1):
    rad[i] = x1min + np.exp(np.log((x1max - x1min + 1.0))*(i-gc)/(Nx1-2.*gc)) - 1.

#print(rad[:5])
#print(rad[255:])  
print(rad[gc:Nx1-gc+1])  
print(rad[Nx1-gc])

w = interp2(rad[gc:Nx1-gc+1])

print (len(w))
#print (w[-1])

P0 = h0*rho0*(gamma-1.)/gamma


rho_atm = 0.1*rho0*(w[-1]/h0)**n
P_atm = P0*(rho_atm/rho0)**gamma
print( rho_atm, P_atm )

keep = sys.stdout
f = open(output, 'w')
sys.stdout = f


print("###############PARAM###############")
print("%e\n%d"% (0.0,Nx1+1-2*gc ))
print("SPHERICAL")
print("###################################")

for i in range(len(w)):
    rho = rho0*(w[i]/h0)**n
    P = P0*(rho/rho0)**gamma
    print("%e %e %e %e" % ( rad[i+gc], rho, P, 0.0 ))

sys.stdout = keep
f.close()

quit()

########## plotting
# font size for the figure
fs = 12

plt.rcParams.update({'font.size': fs})

fig = plt.figure()
ax1 = fig.add_subplot(111)

#plt.xscale('log') 
#plt.yscale('log') 



plt.plot(r, h,color='blue')
plt.plot(rad[gc:Nx1-gc], fit , '--', color='red')

h0 = 1-1/r[-1]
#plt.plot(r, 1 + (1/r -1)/h0 , '--', color='red')

#axes = plt.gca()
#axes.set_xlim([0,5])
#axes.set_ylim(0.8,1.1)
    
#plt.legend()
plt.xlabel('$r/r_0$')
plt.ylabel('$m/M$')
plt.show()

######################################


