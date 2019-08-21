#!/bin/bash

python2.7<<EOF
#######################################################
# 
# Calculates numerical values for Michel accretion model
# using Brent's method
#
# This program is executed as:  
#        ./michel gamma ainf rmin rmax N
# where:
#        gamma is the polytropic index
#        ainf is the asymptotic value of the sound speed
#        rmin is the minimum radius
#        rmax is the maximum radius
#        N is the number of points
# Note that for special values of gamma you can
# enter 1 for gamma == 1
# enter 2 for gamma == 4/3
# enter 3 for gamma == 5/3
# otherwise, provide explicit value of gamma 
# as a floating-point number between 1 and 2
#
# The output columns are:
# radius  density  pressure  velocity  sound_speed
#
# and units are such that:
# rho_infty = 1 (density at infinity)
# G = c = M = 1 
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
import time
import datetime

poly = 3
ainf = 0.01
rmin = 2.5
rmax = 20.0
npoint = 400
gamma = 4.0/3.0

if ainf > np.sqrt(gamma - 1.) : 
    print "ainf should be less than amax =", np.sqrt(gamma - 1.)
    quit()

if rmin<0. or rmax<0. : 
    print "rmin and rmax have to be positive numbers"
    quit()

amax = np.sqrt(gamma - 1.)

# n
n = 1./(gamma - 1.)

# asymptotic value for the enthalpy
winf = 1./(1. - n*ainf**2)

# critical enthalpy
Psi = np.arccos(1.5*(gamma - 2./3.)**(-1.5)*(gamma - 1.)/winf)/3.
wc = 2.*winf*np.sqrt(gamma - 2./3.)*np.sin(Psi + np.pi/6.)

# critical sound speed
ac = np.sqrt((gamma-1.)*(1. - 1./wc))

# critical radial velocity
vc = np.sqrt((1-winf**2/wc**2)/3.)

vc = winf*np.sqrt((wc - 1.)/(n*wc**3))


# critical density
rhoc = ((wc - 1.)/(winf - 1.))**n

# critical radius
rc = 0.5/vc**2

# accretion rate
acc = rc**2*vc*rhoc

# K constant in the polytropic
Kpoly = winf*ainf**2/gamma

rad = np.linspace(rmin, rmax, num = npoint, endpoint = True)

sound = np.zeros(npoint)

ic = int(npoint*(rc-rmin)/(rmax - rmin))

if (ic > 0 and ic < npoint ):
  if (rad[ic]< rc) :
    print "There seems to be a problem with the intervals, "
    print "please try changing the total number of points. "
    quit()

#print "# This file was generated on " , datetime.datetime.now()
print "# Michel's model for spherical accretion "
print "# Model parameters: gamma =",gamma,", a_inf =",ainf
print "# a_max =",amax
print "# Critical values: r_c =",rc,", v_c =",vc 
print "#                  a_c =",ac,", rho_c =",rhoc
print "# Accretion rate: alpha =",acc 
print "# Units: length = M, velocity = c, density = rho_inf"
print "# Columns: radius, density, velocity, sound_speed"


def f(a) :
  return r**3*(np.fabs(n*a**2*(1. - n*ainf**2)))**(2.*n) * (r - 2. - r*((1. - n*a**2)/(1. - n*ainf**2))**2) \
  + acc**2*(np.fabs(n*ainf**2*(1. - n*a**2)))**(2.*n)

# if poly == "1" : 
#     def f(a) :
#       return r**3*(np.fabs(n*a**2*(1. - n*ainf**2)))**(2.*n) * (r - 2. - r*((1. - n*a**2)/(1. - n*ainf**2))**2) \
#       + acc**2*(np.fabs(n*ainf**2*(1. - n*a**2)))**(2.*n)
# else :
#     def f(a) :
#       return r**3*(np.fabs(n*a**2*(1. - n*ainf**2)))**(2.*n) * (r - 2. - r*((1. - n*a**2)/(1. - n*ainf**2))**2) \
#       + acc**2*(np.fabs(n*ainf**2*(1. - n*a**2)))**(2.*n)


""" Finding roots for r < rc """
if (ic > 0) :
  r = rad[0]
  aa = ac
  ab = amax
  s = np.sign(f(aa))
  for i in range(0,min(ic,npoint)):
    r = rad[i]
    while 1 == 1: 
      sign = np.sign(f(ab))
      if (sign != s ) : break
      ab = 0.5*(aa + ab)
    a = sy.brentq(f, aa, ab)
    sound[i] = a
    ab = a

""" Finding roots for r > rc """
if (ic < npoint) :
  ab = ac  
  for i in range(max(0,ic),npoint):
    r = rad[i]
    aa = ainf
    s = np.sign(f(ab))
    while 1 == 1: 
      sign = np.sign(f(aa))
      if (sign != s ) : break
      aa = 0.5*(aa + ab)
    a = sy.brentq(f, aa, ab)
    sound[i] = a
    ab = a    

""" Print results """
for i in range(npoint):
  r = rad[i]
  a = sound[i]
  w = 1./(1. - n*a**2)
  den = ((w - 1.)/(winf - 1.))**n
  press = Kpoly*den**gamma
  vel = acc/(den*r**2) # U^r

  lapse = np.sqrt(1 - 2/r) # Lapse function alpha
  br = 0.0 # Shift vector beta^r
  gtt = - (1 - 2/r) # g_tt
  gTT = - 1/(1 - 2/r) # g^tt
  gtr = 0.0 # g_tr
  gTR = 0.0 # g^tr
  grr = 1/(1 - 2/r) # g_rr
  gRR = 1 - 2/r # g^rr

  Uur = vel # U^r
  Udt = np.sqrt(Uur*Uur - gtt) # U_t
  Udr = grr*Uur # U_r
  Uut = gTT*Udt # U_t
 
  W = lapse*Uut # Lorentz factor

  vdr = (Udr / W) # v_r
  vur = vdr/grr # v^r

  vr = np.sqrt(vur*vdr) # sqrt(v_r * v^r)

  print ("%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e" % ( r, den, press, vel, a, vdr, vur, vr))

#plt.plot(rad,sound,'bo')
#plt.show()
EOF
