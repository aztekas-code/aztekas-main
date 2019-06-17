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
from scipy.interpolate import griddata

matplotlib.rcParams['axes.linewidth'] = 0.3
matplotlib.rcParams['xtick.major.width'] = 0.3
matplotlib.rcParams['ytick.major.width'] = 0.3
matplotlib.rcParams['xtick.major.size'] = 2.0
matplotlib.rcParams['ytick.major.size'] = 2.0
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
matplotlib.rcParams['xtick.top'] = False
matplotlib.rcParams['ytick.right'] = False
matplotlib.rcParams['axes.facecolor'] = 'black'

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
fontsize = 12
orientation = 'h'

# Mesh grid  
nr = np.linspace(np.min(rr), np.max(rr), Nx1)
nt = np.linspace(np.min(tt), np.max(tt), Nx2)
R,T = np.meshgrid(nr,nt)

# Turn radial grid points into (r,z)
r = R*50/13
R = r * np.sin(T)
Z = r * np.cos(T)

# LaTeX text
plt.rc('font', family='serif')
plt.rc('text', usetex=True)

# Figure size
plt.figure(figsize=(12,6),dpi=300)
fig, ax = plt.subplots(1, 1)

# Colorbar limits
cbar_min = -10.2/10#np.min(np.log10(n)) # Min
cbar_max = -8.0/10#np.max(np.log10(n)) # Max
cbarlabels = np.linspace(cbar_min,cbar_max,num=7, endpoint=True) # num of ticks in the colorbar
levels = np.linspace(cbar_min,cbar_max,400) # Levels
cmap = plt.get_cmap('afmhot_r') # Colormap
norm = BoundaryNorm(levels, ncolors=cmap.N) # Normalization
# RdYlGn
# YlOrBr
# RdBu
# gist_heat
# afmhot
# jet
# RdGy
# inferno

# R and Z labels 
#rlabels = np.linspace((np.min(R)), (np.max(R)), num=5, endpoint=True) # num of ticks in R axis
#zlabels = np.linspace((np.min(Z)), (np.max(Z)), num=6, endpoint=True) # num of ticks in Z axis
rlabels = np.linspace((np.min(-15)), (np.max(15)), num=7, endpoint=True) # num of ticks in R axis
zlabels = np.linspace((np.min(-30)), (np.max(20)), num=6, endpoint=True) # num of ticks in Z axis

# Contour
if orientation == 'v':
    cn = ax.contourf(R,Z,np.log10(n)/10,cmap=cmap,levels=levels,norm=norm)
    cn = ax.contourf(R,Z,np.log10(n)/10,cmap=cmap,levels=levels,norm=norm)
    cn = ax.contourf(R,Z,np.log10(n)/10,cmap=cmap,levels=levels,norm=norm)

    cn = ax.contourf(-R,Z,np.log10(n)/10,cmap=cmap,levels=levels,norm=norm)
    cn = ax.contourf(-R,Z,np.log10(n)/10,cmap=cmap,levels=levels,norm=norm)
    cn = ax.contourf(-R,Z,np.log10(n)/10,cmap=cmap,levels=levels,norm=norm)

if orientation == 'h':
    cn = ax.contourf(-Z,R,np.log10(n)/10,cmap=cmap,levels=levels,norm=norm)
    cn = ax.contourf(-Z,R,np.log10(n)/10,cmap=cmap,levels=levels,norm=norm)
    cn = ax.contourf(-Z,R,np.log10(n)/10,cmap=cmap,levels=levels,norm=norm)

    cn = ax.contourf(-Z,-R,np.log10(n)/10,cmap=cmap,levels=levels,norm=norm)
    cn = ax.contourf(-Z,-R,np.log10(n)/10,cmap=cmap,levels=levels,norm=norm)
    cn = ax.contourf(-Z,-R,np.log10(n)/10,cmap=cmap,levels=levels,norm=norm)


# Vector plotting
vector = 0

if vector == 1:
    qx = 10
    qy = 10

    uut= u/np.sqrt(u**2 + v**2)
    vvt= v/np.sqrt(u**2 + v**2)
    uu = np.nan_to_num(u)
    vv = np.nan_to_num(v)

    alfa = 1.0/np.sqrt(1 + 2/r)
    beta = (2.0/r)*(1/(1 + 2/r))

    uu = (alfa*uu/(1 + 2/r)) - beta
    vv = (alfa*vv/(r*r))
    ur = uu*R/r + Z*vv
    uz = uu*Z/r - R*vv

    if orientation == 'v':
        plt.quiver(R[::qx,::qy],Z[::qx,::qy],ur[::qx,::qy],uz[::qx,::qy],
                    units='width',scale=10,cmap=cmap)
        plt.quiver(-R[::qx,::qy],Z[::qx,::qy],-ur[::qx,::qy],uz[::qx,::qy],
                    units='width',scale=10,cmap=cmap)
    if orientation == 'h':
        plt.quiver(Z[::qx,::qy],R[::qx,::qy],uz[::qx,::qy],ur[::qx,::qy],
                    units='width',scale=10,cmap=cmap)
        plt.quiver(Z[::qx,::qy],-R[::qx,::qy],uz[::qx,::qy],-ur[::qx,::qy],
                    units='width',scale=10,cmap=cmap)

# Stream plotting
stream = 1

if stream == 1:
    alfa = 1.0/np.sqrt(1 + 2/r)
    beta = (2.0/r)*(1/(1 + 2/r))

    uu = (alfa*u/(1 + 2/r)) - beta # dr/dt
    vv = (alfa*v/(r*r))            # dtheta/dt

    vr = uu*R/r + Z*vv             # dR/dt
    vz = uu*Z/r - R*vv             # dz/dt 

    nr = np.linspace(np.min(R),np.max(R),Nx1)
    nz = np.linspace(np.min(Z),np.max(Z),Nx2)
    RR, ZZ = np.meshgrid(nr,nz)

    pr = R.flatten()
    pz = Z.flatten()
    pu = vr.flatten()
    pv = vz.flatten()

    gu = griddata((pr,pz),pu,(RR,ZZ))
    gv = griddata((pr,pz),pv,(RR,ZZ))

    if orientation == 'v':
        nr = np.linspace(np.min(R),np.max(R),Nx1)
        nz = np.linspace(np.min(Z),np.max(Z),Nx2)
        RR, ZZ = np.meshgrid(nr,nz)

        pr = R.flatten()
        pz = Z.flatten()
        pu = vr.flatten()
        pv = vz.flatten()

        gu = griddata((pr,pz),pu,(RR,ZZ))
        gv = griddata((pr,pz),pv,(RR,ZZ))

        ax.streamplot(nr,nz,gu,gv,density=1,cmap=cmap,linewidth=0.5,color='k')
        ax.streamplot(-nr,nz,-gu,gv,density=1,cmap=cmap,linewidth=0.5,color='k')

    if orientation == 'h':
        nr = np.linspace(np.min(R),np.max(R),Nx1)
        nz = np.linspace(np.min(Z),np.max(Z),Nx2)
        ZZ, RR = np.meshgrid(nz,nr)

        pr = R.flatten()
        pz = Z.flatten()
        pu = vr.flatten()
        pv = vz.flatten()

        gu = griddata((pz,pr),pu,(ZZ,RR))
        gv = griddata((pz,pr),pv,(ZZ,RR))

#        ax.streamplot(-nz, nr,-gv, gu,density=5,cmap=plt.cm.binary,linewidth=0.5,color=np.sqrt(gv*gv + gu*gu))
#        ax.streamplot(-nz,-nr,-gv,-gu,density=5,cmap=plt.cm.binary,linewidth=0.5,color=np.sqrt(gv*gv + gu*gu))
        ax.streamplot(-nz, nr,-gv, gu,density=[2,1],color='k',linewidth=np.sqrt(gv*gv + gu*gu))
        ax.streamplot(-nz,-nr,-gv,-gu,density=[2,1],color='k',linewidth=np.sqrt(gv*gv + gu*gu))

# Set zoom
if orientation == 'v':
    plt.xlim(-15,15)
    plt.ylim(-35,20)
if orientation == 'h':
    plt.xlim(-35,20)
    plt.ylim(-15,15)

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

# Colorbar positon (right,left,top,bottom)
cbpos = "right"
if (cbpos == "right") or (cbpos == "left"):
    cbor = 'vertical'
    rotation = 270
if (cbpos == "top") or (cbpos == "bottom"):
    cbor = 'horizontal'
    rotation = 0

divider = make_axes_locatable(ax)
cax = divider.append_axes(cbpos,size="5%",pad=0.05)
cbarn = fig.colorbar(cn,orientation=cbor,cax=cax)
cbarn.set_label(r'Density $\log (\rho / \rho_\infty)$',rotation=rotation,fontsize=fontsize,labelpad=15)
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
