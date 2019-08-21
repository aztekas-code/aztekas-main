#!/usr/bin/env python 

import sys, math
import numpy as np
import matplotlib
import linecache
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib import colors, ticker, cm
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.interpolate import griddata

#########################################
# Matplotlib params and other variables #
#########################################

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
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True

# Fontsize and orientation
fontsize = 12
orientation = 'c'

# LaTeX text
plt.rc('font', family='serif')
plt.rc('text', usetex=True)

##########################
# Input/output arguments #
##########################

fname = sys.argv[1] # Filename
oname = sys.argv[2] # Number output

######################
# Reading from fname #
######################

#Time reading
time = float(linecache.getline(fname,2))

#Mesh reading
Nx1 = int(linecache.getline(fname,3))
Nx2 = int(linecache.getline(fname,4))
 
# Data reading
x1, x2, n, p, u, v = np.loadtxt(fname,skiprows=5,unpack=True)

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
nx1 = np.linspace(x1.min(), x1.max(), Nx1)
nx2 = np.linspace(x2.min(), x2.max(), Nx2)
X1,X2 = np.meshgrid(nx1,nx2)

###############################
## Turn (r,theta) into (R,z) ##
###############################

r = X1;
t = X2;

X1 = r * np.sin(t);
X2 = r * np.cos(t);

###############################################
############## Contour graphic ################
###############################################

# Figure size
plt.figure(figsize=(18,18),dpi=300)
fig, ax = plt.subplots(1, 1)

# Graph limits
x1min = -5.5#-X1.max()
x1max =  5.5#X1.max()

x2min = -5.5#X2.min()
x2max =  5.5#X2.max()

# Colorbar limits
cbar_min = n.min() # Min
cbar_max = n.max() # Max

# Contour levels
levels = np.linspace(cbar_min,cbar_max,400) # Levels

# Contour colormap 
cmap = plt.get_cmap('YlOrBr') # Colormap
# RdYlGn
# YlOrBr
# RdBu
# gist_heat
# afmhot
# jet
# RdGy

# Contour normalization
norm = BoundaryNorm(levels, ncolors=cmap.N) # Normalization

################
# Contour plot #
################
if orientation == 'v':
    cn = ax.contourf(X1,X2,n,cmap=cmap,levels=levels,norm=norm)
    cn = ax.contourf(-X1,X2,n,cmap=cmap,levels=levels,norm=norm)

if orientation == 'h':
    cn = ax.contourf(X2,X1,n,cmap=cmap,levels=levels,norm=norm)
    cn = ax.contourf(X2,-X1,n,cmap=cmap,levels=levels,norm=norm)

if orientation == 'c':
    cn = ax.contourf(X1,X2,n,cmap=cmap,levels=levels,norm=norm)
    cn = ax.contourf(-X1,X2,n,cmap=cmap,levels=levels,norm=norm)
    cn = ax.contourf(X1,-X2,n,cmap=cmap,levels=levels,norm=norm)
    cn = ax.contourf(-X1,-X2,n,cmap=cmap,levels=levels,norm=norm)

###############
# Vector plot #
###############
vector = 0

if vector == 1:
    qx = 10
    qy = 10

    uu = np.nan_to_num(u)
    vv = np.nan_to_num(v)
    if orientation == 'v':
        q = plt.quiver(X1[::qx,::qy],X2[::qx,::qy],uu[::qx,::qy],vv[::qx,::qy],
                    units='width')
        ax.quiverkey(q,0.9,0.9,1,r'',labelpos='E',coordinates='figure')
        q = plt.quiver(-X1[::qx,::qy],X2[::qx,::qy],-uu[::qx,::qy],vv[::qx,::qy],
                    units='width')
        ax.quiverkey(q,0.9,0.9,1,r'',labelpos='E',coordinates='figure')
    if orientation == 'h':
        q = plt.quiver(X2[::qx,::qy],X1[::qx,::qy],vv[::qx,::qy],uu[::qx,::qy],
                    units='width')
        ax.quiverkey(q,0.9,0.9,1,r'',labelpos='E',coordinates='figure')
        q = plt.quiver(X2[::qx,::qy],-X1[::qx,::qy],vv[::qx,::qy],-uu[::qx,::qy],
                    units='width')
        ax.quiverkey(q,0.9,0.9,1,r'',labelpos='E',coordinates='figure')

###############
# Stream plot #
###############
stream = 0

if stream == 1:
    if orientation == 'v':
        sx1 = np.linspace(X1.min(),X1.max(),Nx1)
        sx2 = np.linspace(X2.min(),X2.max(),Nx2)
        sX1, sX2 = np.meshgrid(sx1,sx2)

        px1 = X1.flatten()
        px2 = X2.flatten()
        pu = u.flatten()
        pv = v.flatten()

        gu = griddata((px1,px2),pu,(sX1,sX2))
        gv = griddata((px1,px2),pv,(sX1,sX2))

        ax.streamplot( sx1,sx2, gu,gv,density=[1,2],color='k',linewidth=np.sqrt(gu*gu + gv*gv))
        ax.streamplot(-sx1,sx2,-gu,gv,density=[1,2],color='k',linewidth=np.sqrt(gu*gu + gv*gv))

    if orientation == 'h':
        sx1 = np.linspace(X1.min(),X1.max(),Nx1)
        sx2 = np.linspace(X2.min(),X2.max(),Nx2)
        sX2, sX1 = np.meshgrid(sx2,sx1)

        px1 = X1.flatten()
        px2 = X2.flatten()
        pu = u.flatten()
        pv = v.flatten()

        gu = griddata((px2,px1),pu,(sX2,sX1))
        gv = griddata((px2,px1),pv,(sX2,sX1))

        ax.streamplot(sx2, sx1,gv, gu,density=[2,1],color='k',linewidth=np.sqrt(gu*gu + gv*gv))
        ax.streamplot(sx2,-sx1,gv,-gu,density=[2,1],color='k',linewidth=np.sqrt(gu*gu + gv*gv))

############
# Set ZOOM #
############
if orientation == 'v':
    plt.xlim(x1min,x1max)
    plt.ylim(x2min,x2max)
if orientation == 'h':
    plt.xlim(x2min,x2max)
    plt.ylim(x1min,x1max)
if orientation == 'c':
    plt.xlim(x1min,x1max)
    x2min = -x2max
    plt.ylim(x2min,x2max)

#################
# X1 and X2 ticks #
#################
x1labels = np.linspace(x1min,x1max, num=5, endpoint=True) # num of ticks in X1 axis
x2labels = np.linspace(x2min,x2max, num=5, endpoint=True) # num of ticks in X2 axis

if orientation == 'v':
    plt.xticks(x1labels)
    ax.set_xticklabels(['{:.1f}'.format(x) for x in x1labels],fontsize=fontsize)
    plt.yticks(x2labels)
    ax.set_yticklabels(['{:.1f}'.format(x) for x in x2labels],fontsize=fontsize)
    plt.xlabel(r'$r$',fontsize=fontsize)
    plt.ylabel(r'$z$',fontsize=fontsize)

if orientation == 'h':
    plt.xticks(x2labels)
    ax.set_xticklabels(['{:.1f}'.format(x) for x in x2labels],fontsize=fontsize)
    plt.yticks(x1labels)
    ax.set_yticklabels(['{:.1f}'.format(x) for x in x1labels],fontsize=fontsize)
    plt.xlabel(r'$z$',fontsize=fontsize)
    plt.ylabel(r'$r$',fontsize=fontsize)

if orientation == 'c':
    plt.xticks(x1labels)
    ax.set_xticklabels(['{:.1f}'.format(x) for x in x1labels],fontsize=fontsize)
    plt.yticks(x2labels)
    ax.set_yticklabels(['{:.1f}'.format(x) for x in x2labels],fontsize=fontsize)
    plt.xlabel(r'$r$',fontsize=fontsize)
    plt.ylabel(r'$z$',fontsize=fontsize)

############################################
# Colorbar positon (right,left,top,bottom) #
############################################
cbpos = "right"
if (cbpos == "right") or (cbpos == "left"):
    cbor = 'vertical'
    rotation = 270
if (cbpos == "top") or (cbpos == "bottom"):
    cbor = 'horizontal'
    rotation = 0

cax = inset_axes(ax,width='5%',height="100%",loc = 'lower right',bbox_to_anchor = (0.1,0.0,1,1),bbox_transform = ax.transAxes,borderpad = 0)
cbarn = fig.colorbar(cn,orientation=cbor,cax=cax)
cbarn.set_label(r'$\rho$',rotation=rotation,fontsize=fontsize,labelpad=20)
if (cbpos == "right") or (cbpos == "left"):
   cax.yaxis.set_ticks_position(cbpos) 
   cax.yaxis.set_label_position(cbpos) 
if (cbpos == "top") or (cbpos == "bottom"):
   cax.xaxis.set_ticks_position(cbpos) 
   cax.xaxis.set_label_position(cbpos) 

##################
# Colorbar ticks #
##################
cbax1labels = np.linspace(cbar_min,cbar_max,num=6, endpoint=True) # num of ticks in the colorbar
cbarn.set_ticks(cbax1labels)
if (cbpos == "right") or (cbpos == "left"):
    cbarn.ax.set_yticklabels(['{:.2f}'.format(x) for x in cbax1labels],fontsize=fontsize)
if (cbpos == "top") or (cbpos == "bottom"):
    cbarn.ax.set_xticklabels(['{:.2f}'.format(x) for x in cbax1labels],fontsize=fontsize)

ax.set_aspect('equal')

#Graph name
plt.savefig(oname,dpi=300,bbox_inches="tight")
plt.close()
