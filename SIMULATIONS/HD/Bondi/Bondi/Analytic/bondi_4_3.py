#! /usr/bin/python
#######################################################
#
# Compare Bondi analytic solution for gamma = 4/3
# against 2D simulation data.
# Provide simulation file name as first argument.
# It is assumed a square grid between [0,10]x[0,10],
# with a three-line header
# 
#######################################################

import sys, math
import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy import interpolate
from subprocess import check_output

fname=sys.argv[1]

analytic = 'Bondi_4_3.dat'

fs = 10
##################################
# Analytic solution

Arad = []
Avel = []
Aden = []

with open(analytic) as f:
  data = csv.reader(f, delimiter='\t')
  for line in data:
    r = float(line[0])
    rho = float(line[1])
    vr = float(line[3])
    c = float(line[4])
    Arad.append(r)
    Avel.append(vr)
    Aden.append(rho)
    
##################################
# interpolate analytic solution
avel = interpolate.interp1d(Arad, Avel, kind='cubic')
aden = interpolate.interp1d(Arad, Aden, kind='cubic')

##################################
# Starts reading simulation data
# checks grid size assuming square grid and three-line header
grid_size = int(check_output("wc -l "+fname+" | awk '{print $1}'", shell=True))-3

grid_size = int(np.sqrt(grid_size))

nx = grid_size
nz = grid_size

xx = np.arange(0,nx,1)
zz = np.arange(0,nz,1)

xx = 10.0*xx/nx
zz = 10.0*zz/nz

den = np.zeros((nx,nz))
vel = np.zeros((nx,nz))

vvx = np.zeros((nx,nz))
vvz = np.zeros((nx,nz))

with open(fname) as f:
  data = csv.reader(f, delimiter=' ')
  nline = 0
  # skip three-line header
  next(data, None)  
  next(data, None)  
  next(data, None)    
  for line in data:
    i = nline / nx
    j = nline % nx
    x = float(line[0])
    z = float(line[1])
    rho = float(line[2])
    P = float(line[3])    
    vx = float(line[4])
    vz = float(line[5])
    r = x*x + z*z
    cs = np.sqrt(P/rho)
    if x > 0.0 and z > 0.0 :
      r = np.sqrt(r)
      vr = (x*vx + z*vz)/r
      den[j][i] = rho
      vel[j][i] = -vr
    nline += 1

int_den = interpolate.interp2d(xx, zz, den, kind='cubic')
int_vel = interpolate.interp2d(xx, zz, vel, kind='cubic')

# radial partition taken for plotting numerical data
rad = np.arange(0.1, 10., .1)

# number of rays used for averaging
nth = 200
dth = 0.5/nth
th = np.arange(dth, 0.5 - dth, dth)
th = th*np.pi

nden =[]
nvel =[]
verror = []
derror = []

for i in range(len(rad)) :
  R0 = rad[i]
  aux_den = 0.
  aux_vel = 0.  
  for j in range(len(th)) :
    x = R0*np.sin(th[j])
    z = R0*np.cos(th[j])
    aux_den += int_den(x,z)[0]
    aux_vel += int_vel(x,z)[0]
  aux_den = aux_den/len(th)
  aux_vel = aux_vel/len(th)
  nden.append(aux_den)  
  nvel.append(aux_vel)  
  derror.append(aux_den/aden(R0)-1.)
  verror.append(aux_vel/avel(R0)-1.)

########################
#### plot velocity #####
########################

fig1 = plt.figure()
ax1 = fig1.add_subplot(211)

ax1.semilogy(Arad[2:],Avel[2:],'b-',lw=2,label='Bondi')
ax1.plot(rad,nvel, 'rx', label='simulation')

y = [0.001,0.01,0.1,1.0,10.0]
plt.yticks(y, y)

legend = plt.legend(loc='best', prop={'size': fs},handlelength=3, frameon=False)
plt.setp(legend.get_title(),fontsize=fs)

plt.ylabel('$|v_r/c_o|$',fontsize=fs)

# hide xticks
plt.setp(ax1.get_xticklabels(), visible=False) 

axes = plt.gca()
axes.set_xlim([0.0,10.0])
axes.set_ylim([0.005,10.0])

plt.setp(axes.get_xticklabels(),fontsize=fs)
plt.setp(axes.get_yticklabels(),fontsize=fs)


#  error plot  #
ax2 = fig1.add_subplot(212) 
ax2.plot(rad,verror,"r.")
plt.ylabel('Relative error', fontsize=fs)
plt.xlabel('$r/r_B$',fontsize=fs)

axes = plt.gca()
axes.set_xlim([0.0,10.0])

plt.setp(axes.get_xticklabels(),fontsize=fs)
plt.setp(axes.get_yticklabels(),fontsize=fs)

plt.subplots_adjust(hspace=0.1)

plt.savefig('bondi_4_3_velocity.png', bbox_inches='tight', pad_inches=0.0,dpi=200)

########################
##### plot density #####
########################

fig2 = plt.figure()
ax1 = fig2.add_subplot(211)

ax1.semilogy(Arad[2:],Aden[2:],'b-',lw=2,label='Bondi')
ax1.plot(rad,nden, 'rx', label='simulation')

y = [1.0,2.0,4.0,10.0,20.0,40.0]
plt.yticks(y, y)

legend = plt.legend(loc='best', prop={'size': fs},handlelength=3, frameon=False)
plt.setp(legend.get_title(),fontsize=fs)


#plt.xlabel('$r/r_B$',fontsize=fs)
plt.ylabel('$\\rho/\\rho_o$',fontsize=fs)

# hide xticks
plt.setp(ax1.get_xticklabels(), visible=False) 

axes = plt.gca()
axes.set_xlim([0.0,10.0])
axes.set_ylim([1.0,60.0])

plt.setp(axes.get_xticklabels(),fontsize=fs)
plt.setp(axes.get_yticklabels(),fontsize=fs)

#  error plot  #
ax2 = fig2.add_subplot(212)

ax2.plot(rad,derror,'r.')
plt.ylabel('Relative error', fontsize=fs)
plt.xlabel('$r/r_B$',fontsize=fs)

axes = plt.gca()
axes.set_xlim([0.0,10.0])

plt.setp(axes.get_xticklabels(),fontsize=fs)
plt.setp(axes.get_yticklabels(),fontsize=fs)

plt.subplots_adjust(hspace=0.1)

plt.savefig('bondi_4_3_density.png', bbox_inches='tight', pad_inches=0.0,dpi=200)

#plt.show()
