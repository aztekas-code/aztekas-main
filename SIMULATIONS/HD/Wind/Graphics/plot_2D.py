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

if len(sys.argv) < 4:
    print(" ")
    print("Python-Matplotlib for 2D aztekas contour-plots")
    print(" ")
    print("Execute this file in your terminal as")
    print("   $ ./plot_2D.py file.dat output.ext scale")
    print("where")
    print("   file.dat: output file from aztekas")
    print("   output.ext: name of the output, with")
    print("               'ext' the extension (png,eps,pdf)")
    print("   scale: linear or log (lin,log)")
    print(" ")
    sys.exit()

#########################################
# Matplotlib params and other variables #
#########################################

lwidth = 0.3
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
matplotlib.rcParams['axes.facecolor'] = 'black'

# Fontsize and orientation
fontsize = 12
orientation = 'h'

# LaTeX text
plt.rc('font', family='serif')
plt.rc('text', usetex=True)

##########################
# Input/output arguments #
##########################

file_name = sys.argv[1] # Filename
ouput_name = sys.argv[2] # Number output
scale = sys.argv[3] # Scale

##########################
# Reading from file_name #
##########################

#Time reading
time = float(linecache.getline(file_name,2))

#Mesh reading
Nx1   = int(linecache.getline(file_name,3))
Nx2   = int(linecache.getline(file_name,4))
COORD = str(linecache.getline(file_name,5)).rstrip('\n')

# Data reading
U = np.loadtxt(file_name,skiprows=6,unpack=True)

# Data reshape into a Nx1XNx2 matrix
x1 = U[0].reshape(Nx1,Nx2)
x1 = x1.T

x2 = U[1].reshape(Nx1,Nx2)
x2 = x2.T
 
rho = U[2].reshape(Nx1,Nx2)
rho = rho.T
 
pre = U[3].reshape(Nx1,Nx2)
pre = pre.T
 
vx1 = U[4].reshape(Nx1,Nx2)
vx1 = vx1.T
 
vx2 = U[5].reshape(Nx1,Nx2)
vx2 = vx2.T
 
if len(U) == 7:
    vx3 = U[6].reshape(Nx1,Nx2)
    vx3 = vx3.T
 
# Mesh grid  
nx1 = x1[::Nx2,:]
nx2 = x2[:,::Nx1]
X1,X2 = np.meshgrid(nx1,nx2)

if COORD == 'SPHERICAL' or COORD == 'POLAR':
    r = X1
    t = X2

    X1 = r * np.sin(t)
    X2 = r * np.cos(t)

###############################################
############## Contour graphic ################
###############################################

# Figure size
plt.figure(figsize=(200,100),dpi=300)
fig, ax = plt.subplots(1, 1)

# Graph limits
x1min = X1.min()
if (COORD == 'CYLINDRICAL') or (COORD == 'SPHERICAL'):
    x1min = -X1.max()
x1max = X1.max()

x2min = X2.min()
x2max = X2.max()

#x2min = -5
#x2max = 15
x2min = -7
x2max = 4
x1min = -3
x1max = 3

if scale == 'lin':
    plot = rho
if scale == 'log':
    plot = np.log10(rho)

# Colorbar limits
cbar_min = 0.0#plot.min() # Min
cbar_max = 1.5#plot.max() # Max

# Contour levels
levels = np.linspace(cbar_min,cbar_max,400) # Levels

# Contour colormap 
cmap = plt.get_cmap('CMRmap_r') # Colormap
# RdYlGn
# YlOrBr
# RdBu
# gist_heat
# afmhot
# jet
# RdGy

mask = plot < 0.0

plot[mask] = 0.0

################
# Contour plot #
################
if orientation == 'v':
    cn = ax.contourf(X1,X2,plot,cmap=cmap,levels=levels)
    if (COORD == 'CYLINDRICAL') or (COORD == 'SPHERICAL'):
        cn = ax.contourf(-X1,X2,plot,cmap=cmap,levels=levels)

if orientation == 'h':
    cn = ax.contourf(-X2,X1,plot,cmap=cmap,levels=levels)
    if (COORD == 'CYLINDRICAL') or (COORD == 'SPHERICAL'):
        cn = ax.contourf(-X2,-X1,plot,cmap=cmap,levels=levels)

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
    if orientation == 'h':
        q = plt.quiver(X2[::qx,::qy],X1[::qx,::qy],vv[::qx,::qy],uu[::qx,::qy],
                    units='width')
        ax.quiverkey(q,0.9,0.9,1,r'',labelpos='E',coordinates='figure')

###############
# Stream plot #
###############
stream = 1

if stream == 1:
    u = vx1
    v = vx2/r

    ur = u*X1/r + X2*v
    uz = u*X2/r - X1*v

    if orientation == 'v':
        sx1 = np.linspace(x1min,x1max,Nx1)
        sx2 = np.linspace(x2min,x2max,Nx2)
        sX1, sX2 = np.meshgrid(sx1,sx2)

        px1 = X1.flatten()
        px2 = X2.flatten()
        pu = ur.flatten()
        pv = uz.flatten()

        gu = griddata((px1,px2),pu,(sX1,sX2))
        gv = griddata((px1,px2),pv,(sX1,sX2))

        ax.streamplot( sx1,sx2, gu,gv,density=[1,1],color='k',linewidth=0.1)
        ax.streamplot(-sx1,sx2,-gu,gv,density=[1,1],color='k',linewidth=0.1)

    if orientation == 'h':
        sx1 = np.linspace(x1min,x1max,Nx1)
        sx2 = np.linspace(x2min,x2max,Nx2)
        sX2, sX1 = np.meshgrid(sx2,sx1)

        px1 = X1.flatten()
        px2 = -X2.flatten()
        pu = ur.flatten()
        pv = -uz.flatten()

        gu = griddata((px2,px1),pu,(sX2,sX1))
        gv = griddata((px2,px1),pv,(sX2,sX1))

        ax.streamplot(sx2, sx1,gv, gu,density=[2,1],color='k',linewidth=0.2*np.sqrt(gu*gu + gv*gv))
        ax.streamplot(sx2,-sx1,gv, -gu,density=[2,1],color='k',linewidth=0.2*np.sqrt(gu*gu + gv*gv))

############
# Set ZOOM #
############
if orientation == 'v':
    plt.xlim(x1min,x1max)
    plt.ylim(x2min,x2max)
if orientation == 'h':
    plt.xlim(x2min,x2max)
    plt.ylim(x1min,x1max)

###################
# X1 and X2 ticks #
###################

dens = (x2max-x2min)/(x1max-x1min)
num_tick = 7
num_x1 = int(num_tick)
num_x2 = 6#int(num_tick*dens)

x1labels = np.linspace(x1min,x1max, num=num_x1, endpoint=True) # num of ticks in X1 axis
x2labels = np.linspace(-6,4, num=num_x2, endpoint=True) # num of ticks in X2 axis

if orientation == 'v':
    plt.xticks(x1labels)
    ax.set_xticklabels(['{:.1f}'.format(x) for x in x1labels],fontsize=fontsize)
    plt.yticks(x2labels)
    ax.set_yticklabels(['{:.1f}'.format(x) for x in x2labels],fontsize=fontsize)
    if COORD == 'CARTESIAN':
        plt.xlabel(r'$x$',fontsize=fontsize)
        plt.ylabel(r'$y$',fontsize=fontsize)
    if (COORD == 'CYLINDRICAL') or (COORD == 'SPHERICAL'):
        plt.xlabel(r'$R/\zeta_\mathrm{HL}$',fontsize=fontsize)
        plt.ylabel(r'$z/\zeta_{\mathrm{HL}}$',fontsize=fontsize)

if orientation == 'h':
    plt.xticks(x2labels)
    ax.set_xticklabels(['{:.0f}'.format(x) for x in x2labels],fontsize=fontsize)
    plt.yticks(x1labels)
    ax.set_yticklabels(['{:.0f}'.format(x) for x in x1labels],fontsize=fontsize)
    if COORD == 'CARTESIAN':
        plt.xlabel(r'$x$',fontsize=fontsize)
        plt.ylabel(r'$y$',fontsize=fontsize)
    if (COORD == 'CYLINDRICAL') or (COORD == 'SPHERICAL'):
        plt.ylabel(r'$R/\zeta_\mathrm{HL}$',fontsize=fontsize)
        plt.xlabel(r'$z/\zeta_{\mathrm{HL}}$',fontsize=fontsize)

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

#cax = inset_axes(ax,width='3%',height="100%",loc = 'lower right',bbox_to_anchor = (0.05,0.0,1,1),bbox_transform = ax.transAxes,borderpad = 0)
#cbarn = fig.colorbar(cn,orientation=cbor,cax=cax)
divider = make_axes_locatable(ax)
cax = divider.append_axes(cbpos,size="3%",pad=0.05)
cbarn = fig.colorbar(cn,orientation=cbor,cax=cax)

if scale == 'lin':
    cbarn.set_label(r'$\rho/\rho_0$',rotation=rotation,fontsize=fontsize,labelpad=20)
if scale == 'log':
    cbarn.set_label(r'$\log(\rho/\rho_\infty)$',rotation=rotation,fontsize=fontsize,labelpad=15)

if (cbpos == "right") or (cbpos == "left"):
   cax.yaxis.set_ticks_position(cbpos) 
   cax.yaxis.set_label_position(cbpos) 
if (cbpos == "top") or (cbpos == "bottom"):
   cax.xaxis.set_ticks_position(cbpos) 
   cax.xaxis.set_label_position(cbpos) 

##################
# Colorbar ticks #
##################
cbax1labels = np.linspace(cbar_min,cbar_max,num=5, endpoint=True) # num of ticks in the colorbar
cbarn.set_ticks(cbax1labels)
if (cbpos == "right") or (cbpos == "left"):
    cbarn.ax.set_yticklabels(['{:.2f}'.format(x) for x in cbax1labels],fontsize=fontsize)
if (cbpos == "top") or (cbpos == "bottom"):
    cbarn.ax.set_xticklabels(['{:.2f}'.format(x) for x in cbax1labels],fontsize=fontsize)

ax.set_aspect('equal')

#Graph name
plt.savefig(ouput_name,dpi=300,bbox_inches="tight")
plt.close()
