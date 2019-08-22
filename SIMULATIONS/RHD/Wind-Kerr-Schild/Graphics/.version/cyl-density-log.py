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

# Input/output arguments
fname = sys.argv[1] # Filename
oname = sys.argv[2] # Number output

#Time reading
time = float(linecache.getline(fname,2))

#Mesh reading
Nx1 = int(linecache.getline(fname,3))
Nx2 = int(linecache.getline(fname,4))
 
# Data reading
rr, zz, n, p, u, v = np.loadtxt(fname,skiprows=5,unpack=True)

# Data reshape into a Nx1XNx2 matrix
n = n.reshape(Nx1,Nx2)
n = n.T
 
p = p.reshape(Nx1,Nx2)
p = p.T
 
u = u.reshape(Nx1,Nx2)
u = u.T
 
v = v.reshape(Nx1,Nx2)
v = v.T

# Fontsize and orientation
fontsize = 6
orientation = 'v'

# Mesh grid  
nr = np.linspace(np.min(rr), np.max(rr), Nx1)
nz = np.linspace(np.min(zz), np.max(zz), Nx2)
R,Z = np.meshgrid(nr,nz)

# LaTeX text
plt.rc('font', family='serif')
plt.rc('text', usetex=True)

# Figure size
plt.figure(figsize=(18,18),dpi=300)
fig, ax = plt.subplots(1, 1)

# Colorbar limits
cbar_min = np.min(np.log10(n)) # Min
cbar_max = np.max(np.log10(n)) # Max
cbarlabels = np.linspace(cbar_min,cbar_max,num=11, endpoint=True) # num of ticks in the colorbar
levels = np.linspace(cbar_min,cbar_max,400) # Levels
cmap = plt.get_cmap('YlOrBr') # Colormap
norm = BoundaryNorm(levels, ncolors=cmap.N) # Normalization
# RdYlGn
# YlOrBr
# RdBu
# gist_heat
# afmhot
# jet
# RdGy

# R and Z labels 
rlabels = np.linspace((np.min(-R)), (np.max(R)), num=5, endpoint=True) # num of ticks in R axis
zlabels = np.linspace((np.min(Z)), (np.max(Z)), num=6, endpoint=True) # num of ticks in Z axis

# Contour
if orientation == 'v':
    cn = ax.contourf(R,Z,np.log10(n),cmap=cmap,levels=levels,norm=norm)
    cn = ax.contourf(R,Z,np.log10(n),cmap=cmap,levels=levels,norm=norm)
    cn = ax.contourf(R,Z,np.log10(n),cmap=cmap,levels=levels,norm=norm)

    cn = ax.contourf(-R,Z,np.log10(n),cmap=cmap,levels=levels,norm=norm)
    cn = ax.contourf(-R,Z,np.log10(n),cmap=cmap,levels=levels,norm=norm)
    cn = ax.contourf(-R,Z,np.log10(n),cmap=cmap,levels=levels,norm=norm)

if orientation == 'h':
    cn = ax.contourf(Z,R,np.log10(n),cmap=cmap,levels=levels,norm=norm)
    cn = ax.contourf(Z,R,np.log10(n),cmap=cmap,levels=levels,norm=norm)
    cn = ax.contourf(Z,R,np.log10(n),cmap=cmap,levels=levels,norm=norm)

    cn = ax.contourf(Z,-R,np.log10(n),cmap=cmap,levels=levels,norm=norm)
    cn = ax.contourf(Z,-R,np.log10(n),cmap=cmap,levels=levels,norm=norm)
    cn = ax.contourf(Z,-R,np.log10(n),cmap=cmap,levels=levels,norm=norm)


# Vector plotting
vector = 0

if vector == 1:
    qx = 10
    qy = 10

    uut= u/np.sqrt(u**2 + v**2)
    vvt= v/np.sqrt(u**2 + v**2)
    uu = np.nan_to_num(u)
    vv = np.nan_to_num(v)
    plt.quiver(R[::qx,::qy],Z[::qx,::qy],uu[::qx,::qy],vv[::qx,::qy],
                units='width',scale=10,cmap=cmap)
    plt.quiver(-R[::qx,::qy],Z[::qx,::qy],-uu[::qx,::qy],vv[::qx,::qy],
                units='width',scale=10,cmap=cmap)


# Colorbar orientation
if orientation == 'v':
    # Vertical colorbar
    cbarn = plt.colorbar(cn, orientation='vertical',aspect=40,fraction=0.05)
    cbarn.set_label(r'Density $\log (\rho)$',rotation=270,fontsize=fontsize,labelpad=20)
if orientation == 'h':
    # Vertical colorbar
    cbarn = plt.colorbar(cn, orientation='horizontal',aspect=40,fraction=0.05)
    cbarn.set_label(r'Density $\log (\rho)$',rotation=0,fontsize=fontsize,labelpad=20)

# Colorbar ticks
cbarn.set_ticks(cbarlabels)
if orientation == 'v':
    cbarn.ax.set_yticklabels(['{:.2f}'.format(x) for x in cbarlabels],fontsize=fontsize)
if orientation == 'h':
    cbarn.ax.set_xticklabels(['{:.2f}'.format(x) for x in cbarlabels],fontsize=fontsize)

# R and Z ticks
if orientation == 'v':
    plt.xticks(rlabels)
    ax.set_xticklabels(['{:.1f}'.format(x) for x in rlabels],fontsize=fontsize)
    plt.yticks(zlabels)
    ax.set_yticklabels(['{:.1f}'.format(x) for x in zlabels],fontsize=fontsize)
    plt.xlabel(r'$r$',fontsize=fontsize)
    plt.ylabel(r'$z$',fontsize=fontsize)

if orientation == 'h':
    plt.xticks(zlabels)
    ax.set_xticklabels(['{:.1f}'.format(x) for x in zlabels],fontsize=fontsize)
    plt.yticks(rlabels)
    ax.set_yticklabels(['{:.1f}'.format(x) for x in rlabels],fontsize=fontsize)
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$r$',fontsize=fontsize)
ax.set_aspect('equal')

#Graph name
plt.savefig('N'+oname+'_py.png',dpi=300,bbox_inches="tight")
plt.close()
