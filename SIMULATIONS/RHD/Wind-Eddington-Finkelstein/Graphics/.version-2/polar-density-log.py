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
from mpl_toolkits.axes_grid1 import make_axes_locatable
import linecache

matplotlib.rcParams['axes.linewidth'] = 0.3
matplotlib.rcParams['xtick.major.width'] = 0.3
matplotlib.rcParams['ytick.major.width'] = 0.3
matplotlib.rcParams['xtick.major.size'] = 2.0
matplotlib.rcParams['ytick.major.size'] = 2.0
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
matplotlib.rcParams['xtick.top'] = False
matplotlib.rcParams['ytick.right'] = False

# Input/output arguments
fname = sys.argv[1] # Filename
oname = sys.argv[2] # Number output

#Time reading
time = float(linecache.getline(fname,2))

#Mesh reading
Nx1 = int(linecache.getline(fname,3))
Nx2 = int(linecache.getline(fname,4))
 
# Data reading
rr, tt, n, p, u, v = np.loadtxt(fname,skiprows=5,unpack=True)

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
nt = np.linspace(np.min(tt), np.max(tt), Nx2)
R,T = np.meshgrid(nr,nt)

# Turn radial grid points into (x,y)
X = R * np.cos(T)
Y = R * np.sin(T)

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

# X and Y labels 
xlabels = np.linspace((np.min(X)), (np.max(X)), num=5, endpoint=True) # num of ticks in X axis
ylabels = np.linspace((np.min(Y)), (np.max(Y)), num=6, endpoint=True) # num of ticks in Y axis

# Contour
cn = ax.contourf(X,Y,np.log10(n),cmap=cmap,levels=levels,norm=norm)
cn = ax.contourf(X,Y,np.log10(n),cmap=cmap,levels=levels,norm=norm)
cn = ax.contourf(X,Y,np.log10(n),cmap=cmap,levels=levels,norm=norm)

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

# X and Y ticks
plt.xticks(xlabels)
ax.set_xticklabels(['{:.1f}'.format(x) for x in xlabels],fontsize=fontsize)
plt.yticks(ylabels)
ax.set_yticklabels(['{:.1f}'.format(x) for x in ylabels],fontsize=fontsize)
plt.xlabel(r'$x$',fontsize=fontsize)
plt.ylabel(r'$y$',fontsize=fontsize)

# Colorbar positon (right,left,top,bottom)
cbpos = "right"
if (cbpos == "right") or (cbpos == "left"):
    cbor = 'vertical'
    rotation = 270
if (cbpos == "top") or (cbpos == "bottom"):
    cbor = 'horizontal'
    rotation = 0

divider = make_axes_locatable(ax)
cax = divider.append_axes(cbpos,size="4%",pad=0.05)
cbarn = fig.colorbar(cn,orientation=cbor,cax=cax)
cbarn.set_label(r'Density $\log (\rho)$',rotation=rotation,fontsize=fontsize,labelpad=20)
if (cbpos == "right") or (cbpos == "left"):
   cax.yaxis.set_ticks_position(cbpos) 
   cax.yaxis.set_label_position(cbpos) 
if (cbpos == "top") or (cbpos == "bottom"):
   cax.xaxis.set_ticks_position(cbpos) 
   cax.xaxis.set_label_position(cbpos) 

# Colorbar ticks
cbarn.set_ticks(cbarlabels)
if (cbpos == "right") or (cbpos == "left"):
    cbarn.ax.set_yticklabels(['{:.2f}'.format(x) for x in cbarlabels],fontsize=fontsize)
if (cbpos == "top") or (cbpos == "bottom"):
    cbarn.ax.set_xticklabels(['{:.2f}'.format(x) for x in cbarlabels],fontsize=fontsize)

ax.set_aspect('equal')

#Graph name
plt.savefig(oname,dpi=300,bbox_inches="tight")
plt.close()
