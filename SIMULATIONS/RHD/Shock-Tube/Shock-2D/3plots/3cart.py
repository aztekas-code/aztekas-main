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
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['xtick.top'] = True
matplotlib.rcParams['ytick.right'] = True
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.unicode'] = True

# Fontsize and orientation
fontsize = 8
orientation = 'v'

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
x1, x2, n1, p, u, v = np.loadtxt("./cart-h/shock_20.dat",skiprows=5,unpack=True)
x1, x2, n2, p, u, v = np.loadtxt("./cart-v/shock_20.dat",skiprows=5,unpack=True)
x1, x2, n3, p, u, v = np.loadtxt("./cart-d/shock_20.dat",skiprows=5,unpack=True)

# Data reshape into a Nx1XNx2 matrix
n1 = n1.reshape(Nx1,Nx2)
n1 = n1.T
 
n2 = n2.reshape(Nx1,Nx2)
n2 = n2.T
 
n3 = n3.reshape(Nx1,Nx2)
n3 = n3.T

# Mesh grid  
nx1 = np.linspace(x1.min(), x1.max(), Nx1)
nx2 = np.linspace(x2.min(), x2.max(), Nx2)
X1,X2 = np.meshgrid(nx1,nx2)

###############################################
############## Contour graphic ################
###############################################

# Figure size
plt.figure(figsize=(18,18),dpi=300)
fig, ax = plt.subplots(1, 3)

# Graph limits
x1min = X1.min()
x1max = X1.max()

x2min = X2.min()
x2max = X2.max()

# Colorbar limits
cbar_min = 0.0#np.log10(n.min()) # Min
cbar_max = 1.0#np.log10(n.max()) # Max

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

####################
## X1 and X2 ticks #
####################
################
# Contour plot #
################
if orientation == 'v':
    x2labels = np.linspace(0.2,.8, num=4, endpoint=True) # num of ticks in X1 axis
    x1labels = np.linspace(0.2,.8, num=4, endpoint=True) # num of ticks in X2 axis
    cn = ax[0].contourf(X1,X2,n1,cmap=cmap,levels=levels,norm=norm)
    plt.sca(ax[0])
    plt.xticks(x1labels)
    ax[0].set_xticklabels(['{:.1f}'.format(x) for x in x1labels],fontsize=8)
    plt.yticks(x2labels)
    ax[0].set_yticklabels(['{:.1f}'.format(x) for x in x2labels],fontsize=8)

    cn = ax[1].contourf(X1,X2,n2,cmap=cmap,levels=levels,norm=norm)
    plt.sca(ax[1])
    plt.xticks(x1labels)
    ax[1].set_xticklabels(['{:.1f}'.format(x) for x in x1labels],fontsize=8)
    plt.yticks(x2labels)
#    ax[1].set_yticklabels(['{:.1f}'.format(x) for x in x2labels],fontsize=8)
    ax[1].set_yticklabels([])

    cn = ax[2].contourf(X1,X2,n3,cmap=cmap,levels=levels,norm=norm)
    plt.sca(ax[2])
    plt.xticks(x1labels)
    ax[2].set_xticklabels(['{:.1f}'.format(x) for x in x1labels],fontsize=8)
    plt.yticks(x2labels)
    #ax[2].set_yticklabels(['{:.1f}'.format(x) for x in x2labels],fontsize=8)
    ax[2].set_yticklabels([])

    plt.xlim(x1min,x1max)
    plt.ylim(x2min,x2max)

plt.subplots_adjust(wspace=0.0)

###############
# Vector plot #
###############
#vector = 0
#
#if vector == 1:
#    qx = 10
#    qy = 10
#
#    uu = np.nan_to_num(u)
#    vv = np.nan_to_num(v)
#    if orientation == 'v':
#        q = plt.quiver(X1[::qx,::qy],X2[::qx,::qy],uu[::qx,::qy],vv[::qx,::qy],
#                    units='width')
#        ax.quiverkey(q,0.9,0.9,1,r'',labelpos='E',coordinates='figure')
#    if orientation == 'h':
#        q = plt.quiver(X2[::qx,::qy],X1[::qx,::qy],vv[::qx,::qy],uu[::qx,::qy],
#                    units='width')
#        ax.quiverkey(q,0.9,0.9,1,r'',labelpos='E',coordinates='figure')
#
################
## Stream plot #
################
#stream = 0
#
#if stream == 1:
#    if orientation == 'v':
#        sx1 = np.linspace(X1.min(),X1.max(),Nx1)
#        sx2 = np.linspace(X2.min(),X2.max(),Nx2)
#        sX1, sX2 = np.meshgrid(sx1,sx2)
#
#        px1 = X1.flatten()
#        px2 = X2.flatten()
#        pu = u.flatten()
#        pv = v.flatten()
#
#        gu = griddata((px1,px2),pu,(sX1,sX2))
#        gv = griddata((px1,px2),pv,(sX1,sX2))
#
#        ax.streamplot(sx1,sx2,gu,gv,density=[1,2],color='k',linewidth=np.sqrt(gu*gu + gv*gv))
#
#    if orientation == 'h':
#        sx1 = np.linspace(X1.min(),X1.max(),Nx1)
#        sx2 = np.linspace(X2.min(),X2.max(),Nx2)
#        sX2, sX1 = np.meshgrid(sx2,sx1)
#
#        px1 = X1.flatten()
#        px2 = X2.flatten()
#        pu = u.flatten()
#        pv = v.flatten()
#
#        gu = griddata((px2,px1),pu,(sX2,sX1))
#        gv = griddata((px2,px1),pv,(sX2,sX1))
#
#        ax.streamplot(sx2, sx1,gv,gu,density=[2,1],color='k',linewidth=np.sqrt(gu*gu + gv*gv))
#
##################
## X1 and X2 ticks #
##################
#x1labels = np.linspace(x1min,x1max, num=5, endpoint=True) # num of ticks in X1 axis
#x2labels = np.linspace(x2min,x2max, num=5, endpoint=True) # num of ticks in X2 axis
#
#if orientation == 'v':
#    plt.xticks(x1labels)
#    ax.set_xticklabels(['{:.2f}'.format(x) for x in x1labels],fontsize=fontsize)
#    plt.yticks(x2labels)
#    ax.set_yticklabels(['{:.2f}'.format(x) for x in x2labels],fontsize=fontsize)
#    plt.xlabel(r'$x$',fontsize=fontsize)
#    plt.ylabel(r'$y$',fontsize=fontsize)
#
#if orientation == 'h':
#    plt.xticks(x2labels)
#    ax.set_xticklabels(['{:.2f}'.format(x) for x in x2labels],fontsize=fontsize)
#    plt.yticks(x1labels)
#    ax.set_yticklabels(['{:.2f}'.format(x) for x in x1labels],fontsize=fontsize)
#    plt.xlabel(r'$y$',fontsize=fontsize)
#    plt.ylabel(r'$x$',fontsize=fontsize)
#
#############################################
## Colorbar positon (right,left,top,bottom) #
#############################################
cbpos = "right"
if (cbpos == "right") or (cbpos == "left"):
    cbor = 'vertical'
    rotation = 270
if (cbpos == "top") or (cbpos == "bottom"):
    cbor = 'horizontal'
    rotation = 0

#cax = inset_axes(ax,width='5%',height="100%",loc = 'lower right',bbox_to_anchor = (0.2,0.0,1,1),borderpad = 0)
#cbarn = fig.colorbar(cn,orientation=cbor,cax=cax)
#if (cbpos == "right") or (cbpos == "left"):
#   cax.yaxis.set_ticks_position(cbpos) 
#   cax.yaxis.set_label_position(cbpos) 
#if (cbpos == "top") or (cbpos == "bottom"):
#   cax.xaxis.set_ticks_position(cbpos) 
#   cax.xaxis.set_label_position(cbpos) 

cbarn = fig.colorbar(cn,ax=list(ax[:]),orientation=cbor,label=r'Density',shrink=0.37,pad=0.02)
cbarn.set_label(r'$\rho/\rho_0$',rotation=rotation,fontsize=fontsize,labelpad=20)

##################
# Colorbar ticks #
##################
cbarlabels = np.linspace(cbar_min,cbar_max,num=5, endpoint=True) # num of ticks in the colorbar
cbarn.set_ticks(cbarlabels)
if (cbpos == "right") or (cbpos == "left"):
    cbarn.ax.set_yticklabels(['{:.2f}'.format(x) for x in cbarlabels],fontsize=fontsize)
if (cbpos == "top") or (cbpos == "bottom"):
    cbarn.ax.set_xticklabels(['{:.2f}'.format(x) for x in cbarlabels],fontsize=fontsize)

ax[0].set_aspect('equal')
ax[1].set_aspect('equal')
ax[2].set_aspect('equal')

##Graph name
plt.savefig(oname,dpi=300,bbox_inches="tight")
plt.close()
