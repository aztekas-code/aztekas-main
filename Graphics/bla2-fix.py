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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

lwidth = 0.5
tsize  = 2.0

matplotlib.rcParams['axes.linewidth'] = lwidth
matplotlib.rcParams['xtick.major.width'] = lwidth
matplotlib.rcParams['ytick.major.width'] = lwidth
matplotlib.rcParams['xtick.major.size'] = tsize
matplotlib.rcParams['ytick.major.size'] = tsize
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
matplotlib.rcParams['xtick.top'] = False
matplotlib.rcParams['ytick.right'] = False
matplotlib.rcParams['axes.facecolor'] = 'black'

# Fontsize and orientation
fontsize = 12
orientation = 'h'

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

# Mesh grid  
nr = np.linspace(np.min(rr), np.max(rr), Nx1)
nt = np.linspace(np.min(tt), np.max(tt), Nx2)
r,T = np.meshgrid(nr,nt)

# Turn radial grid points into (r,z)
normalization = 50/13
r = r*normalization
Z = r * np.cos(T)
R = r * np.sin(T)

# LaTeX text
plt.rc('font', family='serif')
plt.rc('text', usetex=True)

# Graph limits
rmin = -15 # R.min()
rmax =  15 # R.max()

zmin = -30 # Z.min()
zmax =  20 # Z.max()

# R and Z labels 
rlabels = np.linspace(rmin,rmax, num=7, endpoint=True) # num of ticks in R axis
zlabels = np.linspace(rmin,zmax, num=6, endpoint=True) # num of ticks in Z axis

# Figure size
plt.figure(figsize=(12,6),dpi=300)
fig, ax = plt.subplots(1, 1)

# R and Z ticks
if orientation == 'v':
    plt.xticks(rlabels)
    ax.set_xticklabels(['{:.1f}'.format(x) for x in rlabels],fontsize=fontsize)
    plt.yticks(zlabels)
    ax.set_yticklabels(['{:.1f}'.format(x) for x in zlabels],fontsize=fontsize)
    plt.xlabel(r'$r/M$',fontsize=fontsize)
    plt.ylabel(r'$z/M$',fontsize=fontsize)

if orientation == 'h':
    plt.xticks(zlabels)
    ax.set_xticklabels(['{:.0f}'.format(x) for x in zlabels],fontsize=fontsize)
    plt.yticks(rlabels)
    ax.set_yticklabels(['{:.0f}'.format(x) for x in rlabels],fontsize=fontsize)
    plt.xlabel(r'$z/M$',fontsize=fontsize)
    plt.ylabel(r'$r/M$',fontsize=fontsize)

# Colorbar limits
cbar_min = -10.2/10# Min
cbar_max = -8.0/10# Max
cbarlabels = np.linspace(cbar_min,cbar_max,num=5, endpoint=True) # num of ticks in the colorbar
levels = np.linspace(cbar_min,cbar_max,400) # Levels
cmap = plt.get_cmap('afmhot_r') # Colormap
norm = BoundaryNorm(levels, ncolors=cmap.N) # Normalization

cn = ax.contourf(-Z,R,np.log10(n)/10,cmap=cmap,levels=levels,norm=norm)
cn = ax.contourf(-Z,R,np.log10(n)/10,cmap=cmap,levels=levels,norm=norm)
cn = ax.contourf(-Z,R,np.log10(n)/10,cmap=cmap,levels=levels,norm=norm)
plt.xlim(-35,20)
plt.ylim(-15,15)

# Colorbar positon (right,left,top,bottom)
cbpos = "right"
if (cbpos == "right") or (cbpos == "left"):
    cbor = 'vertical'
    rotation = 270
if (cbpos == "top") or (cbpos == "bottom"):
    cbor = 'horizontal'
    rotation = 0

#divider = make_axes_locatable(ax)
#cax = divider.append_axes(cbpos,size="5%",pad=0.05)
cax = inset_axes(ax,width="5%",height="45%",loc='upper right',bbox_to_anchor=(0.06,-0.01,1,1),bbox_transform=ax.transAxes,borderpad=0)
cbarn = fig.colorbar(cn,orientation=cbor,cax=cax)
cbarn.set_label(r'$\log(\rho/\rho_\infty)$',rotation=rotation,fontsize=fontsize,labelpad=15)
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

# Colorbar limits
cbar_min = 0.0# Min
cbar_max = 1.0# Max
cbarlabels = np.linspace(cbar_min,cbar_max,num=5, endpoint=True) # num of ticks in the colorbar
levels = np.linspace(cbar_min,cbar_max,400) # Levels
cmap = plt.get_cmap('RdGy') # Colormap
norm = BoundaryNorm(levels, ncolors=cmap.N) # Normalization

cn = ax.contourf(-Z,-R,np.sqrt(u*u/(1 + 2/r) + v*v/(r*r)),cmap=cmap,levels=levels,norm=norm)
cn = ax.contourf(-Z,-R,np.sqrt(u*u/(1 + 2/r) + v*v/(r*r)),cmap=cmap,levels=levels,norm=norm)
cn = ax.contourf(-Z,-R,np.sqrt(u*u/(1 + 2/r) + v*v/(r*r)),cmap=cmap,levels=levels,norm=norm)

# Colorbar positon (right,left,top,bottom)
cbpos = "right"
if (cbpos == "right") or (cbpos == "left"):
    cbor = 'vertical'
    rotation = 270
if (cbpos == "top") or (cbpos == "bottom"):
    cbor = 'horizontal'
    rotation = 0

cax = inset_axes(ax,width="5%",height="45%",loc='lower right',bbox_to_anchor=(0.06,0.01,1,1),bbox_transform=ax.transAxes,borderpad=0)
cbarn = fig.colorbar(cn,orientation=cbor,cax=cax)
cbarn.set_label(r'$v$',rotation=rotation,fontsize=fontsize,labelpad=15)
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

##Graph name
plt.savefig(oname,dpi=300,bbox_inches="tight")
plt.close()
