#! /usr/bin/python
#######################################################
#
# Compute accretion rate accross a spherical shell
# as a function of radius from 2D simulation data.
# Provide polytropic index (gamma) as first argument.
# Provide simulation file name as second argument.
# Provide output name as last argument
# It is assumed an evenly spaced square grid 
# between [0,10]x[0,10], with a three-line header
# 
#######################################################

import sys, math
import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt
from subprocess import check_output
import csv
from scipy import interpolate

gamma= float(sys.argv[1])

fname= sys.argv[2]

out= sys.argv[3]

# font size for the figure
fs = 10

# Bondi accretion rate:
dotM_B = 0.25*(2./(5. - 3.*gamma))**((5. - 3.*gamma)/(2.*(gamma - 1.)))

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

flow = np.zeros((nx,nz))
inflow = np.zeros((nx,nz))
outflow = np.zeros((nx,nz))

with open(fname) as f:
  data = csv.reader(f, delimiter=' ')
  nline = 0
  # skip three-line header
  next(data, None)  
  time = next(data)  
  next(data, None)    
  for line in data:
    i = nline / nx
    j = nline % nx
    x = float(line[0])
    z = float(line[1])
    rho = float(line[2]) 
    vx = float(line[4])
    vz = float(line[5])
    r = x*x + z*z
    if r > 0.0 :
        r = np.sqrt(r)
        vr = (x*vx + z*vz)/r
        flow[j][i] = rho*vr
        if vr < 0.0 :
            inflow[j][i] = rho*vr
        if vr > 0.0 :
            outflow[j][i] = rho*vr
    nline += 1

time = "%.1f" % float(time[0])    

flux = interpolate.interp2d(xx, zz, flow, kind='cubic')
influx = interpolate.interp2d(xx, zz, inflow, kind='cubic')
outflux = interpolate.interp2d(xx, zz, outflow, kind='cubic')

# radial partition taken for compiting accretion rate
rad = np.arange(.2, 10.2, .2)

# number of rays used for averaging
nth = 200
dth = 0.5/nth
th = np.arange(0.0, 0.5 + dth, dth)
th = th*np.pi

FTOT = []
FIN = []
FOUT = []
acc = []

for j in range(len(rad)) :

    R0 = rad[j]
    norm = R0
    IN = []
    OUT = []
    TOT = []

    for i in range(len(th)) :
      x = R0*np.sin(th[i])
      z = R0*np.cos(th[i])

      IN.append(norm*x*influx(x,z)[0])
      OUT.append(norm*x*outflux(x,z)[0])  
      TOT.append(norm*x*flux(x,z)[0])
      
    #plt.plot(th,TOT,label=R0)
      
    numIN = interpolate.interp1d(th,IN, kind='cubic')    
    resIN = integrate.quad(lambda x: numIN(x), 0.0,0.5*np.pi)[0]
    numOUT = interpolate.interp1d(th,OUT, kind='cubic')
    resOUT = integrate.quad(lambda x: numOUT(x), 0.0,0.5*np.pi)[0]
    numTOT = interpolate.interp1d(th,TOT, kind='cubic')
    resTOT = integrate.quad(lambda x: numTOT(x), 0.0,0.5*np.pi)[0]

    FIN.append(resIN)
    FOUT.append(resOUT)
    FTOT.append(-resTOT/dotM_B)    
    acc.append(-(resIN + resOUT)/dotM_B)
        
    #print rad, resIN, resOUT, resIN + resOUT, resOUT/resIN

#plt.show()

acc_min = min(acc)
acc_max = max(acc)
acc_avg = np.average(acc)

print np.average(FTOT)
print acc_min,acc_max,acc_avg
print acc_avg - acc_min , acc_max - acc_avg

figure = plt.figure(figsize=(6,4))
fig = figure.add_subplot(111)
fig.set_aspect('auto')


#plt.plot(rad,FIN, 'b-', lw=1, label='inflow')
#plt.plot(rad,FOUT, 'r-', lw=1, label='outflow')
plt.plot(rad,acc, 'b-x', lw=1) #, label='accretion')

plt.xlabel('$r/r_B$',fontsize=fs)
plt.ylabel('$\dot{M}/\dot{M}_B$',fontsize=fs)

#legend = plt.legend(loc='best', prop={'size': fs},handlelength=3)
#plt.setp(legend.get_title(),fontsize=fs)

axes = plt.gca()
#axes.set_ylim([0,z0max])

plt.setp(axes.get_xticklabels(),fontsize=fs)
plt.setp(axes.get_yticklabels(),fontsize=fs)

#plt.show()

plt.savefig(out+'.png', bbox_inches='tight', pad_inches=0.1, dpi=200)


quit()

# save data file
f = open(out+'.dat','w')

for i in range(len(rad)):
  f.write(str(rad[i]) + " " + str(FIN[i]) + " " + str(FOUT[i]) + " " + str(acc[i]) + " " + str(-FOUT[i]/FIN[i]) +"\n")

