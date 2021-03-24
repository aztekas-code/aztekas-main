#! /usr/bin/python
#######################################################
# 
#
#######################################################

import sys, math
import matplotlib
import linecache
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as color
import scipy.constants as const
from scipy import optimize

import csv

rho0 = 1.0
r0 = 10.0

def time(r):
   dum=np.sqrt(1.-r/r0) 
   dum1 = dum*np.sqrt(r/r0)  + np.arcsin(dum)
   C = np.sqrt(8.*np.pi*rho0/3.)
   return dum1/C
   
def f(r,T) :
   return T - time(r)

def dust(T) :
   r = optimize.brentq(f, 0.1, r0, args=(T))
   rho = rho0*(r0/r)**3
   return r, rho

print (dust(0.5))   

matplotlib.rcParams['text.usetex'] = True
cmap = color.get_cmap('plasma')

#fname = sys.argv[1] # Filename
path = 'DATA/test_'
flist= [1,2,3,4,5]

screen = True

if len(sys.argv) > 1 :
  screen = False
  out = sys.argv[1]

# font size for the figure
fs = 12

plt.rcParams.update({'font.size': fs})

fig = plt.figure()
ax1 = fig.add_subplot(111)
plt.rcParams.update({'font.size': fs})

colors = np.linspace(0.9,0.1,len(flist))

colors = cmap(colors)

ar = []
aden = []

for j in range(len(flist)) :
  
  fname = path+str(flist[j])+".dat"

#Time reading
  tt = float(linecache.getline(fname,2))

  rCW, rhoCW = dust(tt)
  ar.append(rCW)
  aden.append(rhoCW)  

for j in range(len(flist)) :
  
  fname = path+str(flist[j])+".dat"

#Time reading
  tt = float(linecache.getline(fname,2))

  Nr = int(linecache.getline(fname,3))

  rad= np.zeros(Nr)
  den= np.zeros(Nr)
  CW= np.zeros(Nr)
  
  #print (rCW, rhoCW)
  
  with open(fname) as f:
    data = csv.reader(f, delimiter=' ')
  # skip five-line header
    next(data, None)
    next(data, None)
    next(data, None)    
    next(data, None)    
    next(data, None)                
    i = 0
    for line in data:
      rad[i] = float(line[0])
      den[i] = float(line[1])
      if rad[i] < ar[j] :
          CW[i] = aden[j]
      i += 1

  plt.plot(rad, CW, '-', color=colors[j], label="$t={:.1f}$".format(tt),lw=2)
  plt.plot(rad, den,'--', color=colors[j],lw=2)



axes = plt.gca()
#axes.set_ylim([-0.5,5.0])
#axes.set_ylim(-6,4)

plt.xlabel('$r$',fontsize=fs)
plt.ylabel('$\\rho$',fontsize=fs)


legend = plt.legend(loc='best',fontsize=fs, frameon=False)

if screen : 
  plt.show()
else :
  plt.savefig(out, bbox_inches='tight', pad_inches=0.1)
#  plt.savefig(out+'.pdf', bbox_inches='tight', pad_inches=0.1)
