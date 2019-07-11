#! /usr/bin/python
#######################################################
#
# Plot density from simulation data
# Provide simulation file name as first argument.
# Provide output name as second argument
# It is assumed an evenly spaced square grid 
# between [0,10]x[0,10], with a three-line header
# 
#######################################################

import sys, math
import numpy as np
import matplotlib.pyplot as plt
from subprocess import check_output
import csv

fname= sys.argv[1]

out= sys.argv[2]

# font size for the figure
fs = 12

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

vvx = np.zeros((nx,nz))
vvz = np.zeros((nx,nz))
den = np.ones((nx,nz))

# set cutoff for maximum density
rho_max = 1.

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
    rho = float(line[2])
    vx = float(line[4])
    vz = float(line[5])    
    vvx[j][i] = vx
    vvz[j][i] = vz
    rho = np.log10(rho)
    if rho < rho_max :
      den[j][i] = rho
    else :
      den[j][i] = rho_max      
    nline += 1

time = "%.1f" % float(time[0])

fig = plt.figure()

cm = plt.cm.get_cmap('YlOrRd')

plt.xlabel('$R/r_B$',fontsize=fs)
plt.ylabel('$z/r_B$',fontsize=fs)

plt.contour(xx, zz, den, 20, linewidths=0.5, colors='k')
CS = plt.contourf(xx, zz, den, 100, cmap=cm, extend='max')

#vmin=0.1,vmax=2.0,
#CS.cmap.set_under('white')

cbar=plt.colorbar()
cbar.set_label('log($\\rho/\\rho_0$)', labelpad=-40, y=1.06, rotation=0,fontsize=fs)

cbar.ax.tick_params(labelsize=fs)

#plt.clim(0.1,2.0)

plt.streamplot(xx,zz,vvx,vvz,density=1,color='b',linewidth=1.5)

axes = plt.gca()
axes.set_xlim([0,10])
axes.set_ylim([0,10])

x = [0,2,4,6,8,10]
y = [0,2,4,6,8,10]
plt.xticks(x, x)
plt.yticks(y, y)

#plt.text(8, 9, "t = "+time,fontsize=fs)
  
plt.setp(axes.get_xticklabels(),fontsize=fs)
plt.setp(axes.get_yticklabels(),fontsize=fs)

#plt.show()

#plt.savefig('density_'+out+'.png', bbox_inches='tight', pad_inches=0.0,dpi=200)
#plt.savefig('density_'+out+'.eps', bbox_inches='tight', pad_inches=0.0)
plt.savefig('density_'+out+'.pdf', bbox_inches='tight', pad_inches=0.1)
