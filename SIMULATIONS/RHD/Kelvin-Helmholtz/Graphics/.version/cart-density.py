#!/usr/bin/env python 

import sys, math
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib import colors, ticker, cm
from matplotlib.ticker import FormatStrFormatter
import linecache

fname = sys.argv[1]
oname = sys.argv[2]

#Time reading
time = float(linecache.getline(fname,2))

#Mesh reading
Nx1 = int(linecache.getline(fname,3))
Nx2 = int(linecache.getline(fname,4))

#Data reading
xx, yy, n, p, u, v = np.loadtxt(fname,skiprows=5,unpack=True)

n = n.reshape(Nx1,Nx2)
n = n.T
 
p = p.reshape(Nx1,Nx2)
p = p.T
 
u = u.reshape(Nx1,Nx2)
u = u.T
 
v = v.reshape(Nx1,Nx2)
v = v.T
 
nx = np.linspace(np.min(xx), np.max(xx), Nx1)
ny = np.linspace(np.min(yy), np.max(yy), Nx2)
X,Y = np.meshgrid(nx,ny)

plt.rc('font', family='serif')
plt.rc('text', usetex=True)

plt.figure(figsize=(18,12),dpi=300)

fig, ax = plt.subplots(1, 1)

cbar_min = np.min(n)
cbar_max = np.max(n)
cbarlabels = np.linspace(cbar_min,cbar_max,num=5, endpoint=True)

levels = np.linspace(cbar_min,cbar_max,400)
cmap = plt.get_cmap('RdBu')
norm = BoundaryNorm(levels, ncolors=cmap.N)
# RdYlGn
# YlOrBr
# RdBu
# gist_heat
# afmhot
# jet
# RdGy

xlabels = np.linspace((np.min(X)), (np.max(X)), num=5, endpoint=True)
ylabels = np.linspace((np.min(Y)), (np.max(Y)), num=5, endpoint=True)

cn = ax.contourf(X,Y,n,cmap=cmap,levels=levels,norm=norm)
cn = ax.contourf(X,Y,n,cmap=cmap,levels=levels,norm=norm)
cn = ax.contourf(X,Y,n,cmap=cmap,levels=levels,norm=norm)
cn = ax.contourf(X,Y,n,cmap=cmap,levels=levels,norm=norm)

#Vector plotting
#qx = 5
#qy = 5

#uut= u/np.sqrt(u**2 + v**2)
#vvt= v/np.sqrt(u**2 + v**2)
#uu = np.nan_to_num(u)
#vv = np.nan_to_num(v)
#plt.quiver(Rt_x[::qx,::qy],Rt_y[::qx,::qy],uu[::qx,::qy],vv[::qx,::qy],
#       units='width',scale=205.5,cmap='Reds')

fontsize = 6

# Vertical colorbar
cbarn = plt.colorbar(cn,orientation='vertical',aspect=30,fraction=0.05)
# Horizontal colorbar
#cbarn = plt.colorbar(cn,orientation='horizontal',fraction=0.1)

cbarn.set_label(r'Density $\rho$',rotation=270,fontsize=fontsize,labelpad=20)
cbarn.set_ticks(cbarlabels)
# Vertical colorbar
cbarn.ax.set_yticklabels(['{:.2f}'.format(x) for x in cbarlabels],fontsize=fontsize)
# Horizontal colorbar
#cbarn.ax.set_xticklabels(['{:.2f}'.format(x) for x in cbarlabels],fontsize=fontsize)

plt.xticks(xlabels)
ax.set_xticklabels(['{:.1f}'.format(x) for x in xlabels],fontsize=fontsize)
plt.yticks(ylabels)
ax.set_yticklabels(['{:.1f}'.format(x) for x in ylabels],fontsize=fontsize)
plt.xlabel(r'$x$',fontsize=fontsize)
plt.ylabel(r'$y$',fontsize=fontsize)
ax.set_aspect('equal')

#Graph name
plt.savefig('N'+oname+'_py.png',dpi=300,bbox_inches="tight")
#plt.savefig('N'+oname+'_py.eps',dpi=300,bbox_inches="tight")
#plt.savefig('N'+oname+'_py.pdf',bbox_inches="tight")
plt.close()
