#######################################################
# 
# Calculates numerical values for Bondi accretion model
# using parabolic Halley's method
#
# This program is executed as:  
#        ./bondi gamma rmin rmax N
# where:
#        gamma is the polytropic index
#        rmin is the minimum radius
#        rmax is the maximum radius
#        N is the number of points
# Note that for special values of gamma you can
# enter 1 for gamma == 1
# enter 2 for gamma == 4/3
# enter 3 for gamma == 5/3
# otherwise, provide explicit value of gamma 
# as a floating-point number between 1 and 5/3
#
# The output columns are:
# radius  density  velocity  Mach_number
#
# and units are such that:
# rho_infty = 1 (density at infinity)
# a_infty = 1 (sound's speed at infinity)
# r_B = GM/a^2_infty = 1 (Bondi's radius)
#
# author: Emilio Tejeda
# e-mail: etejeda@astro.unam.mx
# date:   2/jun/2018
#
#######################################################

import sys, math
import numpy as np
import scipy.optimize as sy
import matplotlib.pyplot as plt

if len(sys.argv) != 5 :
    print ("This program calculates numerical values \
for Bondi's accretion model \n\
execute as:\n\
  ./bondi gamma rmin rmax N\nwhere:\n\
  gamma is the polytropic index\n\
  rmin is the minimum radius\n\
  rmax is the maximum radius\n\
  N is the number of points\n\
\nNote that for special values of gamma you can\n\
enter 1 for gamma == 1\n\
enter 2 for gamma == 4/3\n\
enter 3 for gamma == 5/3\n\
otherwise, provide explicit value of gamma \n\
as a floating-point number between 1 and 5/3")
    quit()

poly = sys.argv[1]
rmin = float(sys.argv[2])
rmax = float(sys.argv[3])
npoint = int(sys.argv[4])

# choose gamma and check consistency of input
if poly == "1" : 
    gamma = 1.
    acc_B = 0.25*np.exp(1.5)
elif poly == "2" : 
    gamma = 4./3.
    acc_B = 1./np.sqrt(2.)
elif poly == "3" : 
    gamma = 5./3.
    acc_B = 0.25
else :
    gamma = float(poly)
    if gamma > 5./3. :
        print ("gamma must be less than 5/3")
        quit()
    acc_B = 0.25*(2./(5.-3.*gamma))**(0.5*(5.-3.*gamma)/(gamma-1.))

if rmin<0. or rmax<0. : 
    print ("rmin and rmax have to be positive numbers")
    quit()

# sonic radius
r_s = 0.25*(5.-3.*gamma)


if poly == "1" : 
    def f(x):
      return x**2*(np.log(x) - 1./r) + 0.5*acc_B**2/r**4
    def df(x):
      return 2.*x*(np.log(x) - 1./r + 0.5)
    def ddf(x):
      return 2.*(np.log(x) - 1./r + 1.5)
else :
    def f(x):
      return x**(gamma + 1.) - x**2*(1. + (gamma - 1.)/r) \
             + 0.5*(gamma - 1.)*acc_B**2/r**4
    def df(x):
      return (gamma + 1.)*x**gamma - 2.*x*(1. + (gamma - 1.)/r) 
    def ddf(x):
      return gamma*(gamma + 1.)*x**(gamma - 1.) \
             - 2.*(1. + (gamma - 1.)/r) 
    
def vel(rho,r):
    return acc_B/(rho*r**2)

def Mach(rho,r):
    return acc_B/(rho**(0.5*(gamma+1.))*r**2)
    
rad = np.linspace(rmax, rmin, num = npoint, endpoint = True)

M = 10.
rho0= 1.

for i in range(len(rad)) :
  r = rad[i]
  if r > r_s :
    M = 10.
    # making sure solution is on subsonic branch    
    while M > 1.0:
      rho = sy.newton(f, rho0, fprime=df,maxiter=500,fprime2= ddf)
      vv = vel(rho,r)
      M = Mach(rho,r)
      #rho0 = 2.*rho0
      rho0 = rho
  else :
    M = 0.1
    # making sure solution is on supersonic branch    
    while M < 1.0:
      rho = sy.newton(f, rho0, fprime=df,maxiter=500,fprime2= ddf)
      vv = vel(rho,r)
      M = Mach(rho,r)     
      #rho0 = .5*rho0
      rho0 = rho
  print (r, rho, vv, M)
