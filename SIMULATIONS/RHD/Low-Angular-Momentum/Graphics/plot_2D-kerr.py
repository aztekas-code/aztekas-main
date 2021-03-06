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
from scipy.interpolate import griddata, interp1d
from scipy.optimize import fsolve

if len(sys.argv) < 9:
    print(" ")
    print("Python-Matplotlib for 2D aztekas contour-plots")
    print(" ")
    print("Execute this file in your terminal as")
    print("   $ ./plot_2D.py file.dat output.ext scale a plot plot_max gamma Rmax")
    print("where")
    print("   file.dat: output file from aztekas")
    print("   output.ext: name of the output, with")
    print("               'ext' the extension (png,eps,pdf)")
    print("   scale: linear or log (lin,log)")
    print("   a: spin parameter 0..1") 
    print("   plot: rho, pre, vel")
    print("   plot_max: max value of the graph")
    print("   gamma: adiabatic index")
    print(" ")
    sys.exit()

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

# Fontsize and orientation
fontsize = 11
orientation = 'v'

# LaTeX text
plt.rc('font', family='serif')
plt.rc('text', usetex=True)

##########################
# Input/output arguments #
##########################

file_name  = sys.argv[1] # Filename
ouput_name = sys.argv[2] # Number output
scale      = sys.argv[3] # Scale
a          = float(sys.argv[4]) # Spin
plot_var   = sys.argv[5] # Variable to plot
plot_max   = float(sys.argv[6]) # Variable plot max
gamma      = float(sys.argv[7]) # Variable plot max
Rmax       = float(sys.argv[8]) # Variable plot max

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

h   = 1.0 + gamma * pre / (rho * (gamma - 1.0))
c_s = np.sqrt(gamma * pre / (h * rho))
 
if len(U) == 7:
    vx3 = U[6].reshape(Nx1,Nx2)
    vx3 = vx3.T
 
# Mesh grid  
nx1 = x1[::Nx2,:]
nx2 = x2[:,::Nx1]
X1,X2 = np.meshgrid(nx1,nx2)

if COORD == 'SPHERICAL' or COORD == 'POLAR':
    r  = X1
    th = X2

    X1 = np.sqrt(r**2 + a**2) * np.sin(th)
    X2 = r * np.cos(th)

    # If Kerr-Schild
    rho2  = r**2. + (a**2.)*(np.cos(th)**2.)
    Delta = r**2. - 2.*r + a**2.
    Sigma = (r**2. + a**2.)**2. - (a**2.)*Delta*(np.sin(th)**2.)

    alpha   = 1./np.sqrt(1. + 2.*r/rho2)
    
    beta_r   = 2.*r/rho2
    beta_th  = 0.
    beta_phi = 0.

    betar    = (2.*r/rho2)*(1./(1. + 2.*r/rho2))
    betath   = 0.
    betaphi  = 0.

    gam_rr     = 1. + 2.*r/rho2
    gam_rth    = 0
    gam_rphi   = -a*(1. + 2.*r/rho2)*np.sin(th)**2
    gam_thth   = rho2
    gam_thphi  = 0
    gam_phiphi = (np.sin(th)**2.)*((a**2.)*(1. + 2.*r/rho2)*(np.sin(th)**2.) + rho2)

    gamrr     = ((a**2)*(rho2 + 2.*r)*(np.sin(th)**2.) + rho2**2.)/(rho2*(rho2 + 2.*r))
    gamrth    = 0.
    gamrphi   = a/rho2
    gamthth   = 1./rho2
    gamthphi  = 0
    gamphiphi = 1./(rho2*np.sin(th+0.000001)**2.)

    # Define v_i
    v_r  = vx1
    v_th = vx2 
    if len(U) == 7:
        v_phi = vx3

    # Compute v^i = gam^{ij} v_j
    vr  = gamrr*v_r  + gamrth*v_th
    vth = gamrth*v_r + gamthth*v_th
    if len(U) == 7:
        vr  += gamrphi*v_phi
        vth += gamthphi*v_phi
        vphi = gamrphi*v_r + gamthphi*v_th + gamphiphi*v_phi

    # Compute Horizon Penetrating (HP) vel magnitud square VV_HP = v^i v_i
    VV_HP = vr*v_r + vth*v_th
    if len(U) == 7:
        VV_HP += vphi*v_phi

    # Compute HP Lorentz factor W_HP = 1/sqrt(1 - VV_HP)
    W_HP  = 1./np.sqrt(1. - VV_HP)
    VV_HP = np.sqrt(VV_HP)

    # Compute the four velocity vector U^\mu
    UT   = W_HP/alpha
    Ur   = W_HP*(vr  - betar/alpha)
    Uth  = W_HP*(vth - betath/alpha)
    if len(U) == 7:
        Uphi = W_HP*(vphi - betaphi/alpha)

    UhatT     = UT/np.sqrt(1.0 + 2.0*r/rho2)
    Uhatr     = np.sqrt(1.0 + 2.0*r/rho2)*(Ur + (2.0*r/(rho2 + 2.0*r))*UT - Uphi*a*np.sin(th)**2.0)
    Uhattheta = np.sqrt(rho2)*Uth;
    Uhatphi   = np.sqrt(rho2)*np.sin(th)*Uphi

    vX1 = Ur*X1*r/(r**2. + a**2.) + Uth*np.sqrt(r**2 + a**2)*X2/r
    vX2 = Ur*X2/r - Uth*r*X1/np.sqrt(r**2. + a**.2)

    if plot_var == 'rho':
        rmin = 1 + np.sqrt(1 - a**2) + 0.01
    if plot_var == 'vel':
        rmin = 1 + np.sqrt(1 - a**2) + 0.01
    if plot_var == 'ur':
        rmin = 1 + np.sqrt(1 - a**2) + 0.01
    if plot_var == 'utheta':
        rmin = 1 + np.sqrt(1 - a**2) + 0.01
    if plot_var == 'uphi':
        rmin = 1 + np.sqrt(1 - a**2) + 0.01

    # Convert to non-horizon penetrating
    alpha_NHP  = np.sqrt(rho2*Delta/Sigma)
    beta       = -2.*a*r/Sigma
    gam_rr     = rho2/Delta
    gam_thth   = rho2
    gam_phiphi = Sigma*(np.sin(th)**2.)/rho2

    UT_NHP  = UT - 2.*r*Ur/Delta
    Ur_NHP  = Ur
    Uth_NHP = Uth
    if len(U) == 7:
        Uphi_NHP = Uphi - a*Ur/Delta

    W_NHP  = alpha_NHP*UT_NHP

    vr_NHP  = Ur_NHP/W_NHP
    vth_NHP = Uth_NHP/W_NHP
    vphi_NHP = Uphi_NHP/W_NHP + beta/alpha_NHP

    VV_NHP = np.sqrt(1. - 1./W_NHP**2.)
    vX1 = Ur_NHP*X1*r/(r**2. + a**2.) + Uth_NHP*np.sqrt(r**2 + a**2)*X2/r
    vX2 = Ur_NHP*X2/r - Uth_NHP*r*X1/np.sqrt(r**2. + a**.2)

    sonic = (1.0 - c_s**2.0)*W_NHP**2.0

# User instructions
vel_axis = abs(Ur[0][:])
#vel = interp1d(r[0,:],Ur[0][:],kind='cubic')
stag = np.where(vel_axis == np.amin(vel_axis))
print("Stagnation point S = ",10.0)
print("Velocity of ejection V_ej = ",np.sqrt(VV_NHP[0,Nx1-1]))
print("Lorentz factor ejection W_ej = ",W_NHP[0,Nx1-1])
print("Velocity of injection V_0 = ",np.sqrt(VV_NHP[Nx2-1,Nx1-1]))
print("Lorentz factor injection W_0 = ",W_NHP[Nx2-1,Nx1-1])
print("Ratio Vej/V_0 = ",np.sqrt(VV_NHP[0,Nx1-1]/VV_NHP[Nx2-1,Nx1-1]))
print("Ratio Wej/W_0 = ",np.sqrt(W_NHP[0,Nx1-1]/W_NHP[Nx2-1,Nx1-1]))

###############################################
############## Contour graphic ################
###############################################

# Figure size
plt.figure(figsize=(18,18),dpi=300)
fig, ax = plt.subplots(1, 1)
#plt.text(-9.1666666666,9.16666666666666,"$a = "+str(a)+"$",fontsize=fontsize)

# Black background
#back      = plt.Rectangle((-X1.max() + 1.,0),25,15,color='k',fill=True,lw=0.0,zorder=0)
#ax.add_patch(back)

# Graph limits
x1min = X1.min()
if (COORD == 'CYLINDRICAL') or (COORD == 'SPHERICAL'):
    x1min = -X1.max()
x1max = X1.max()

x2min = X2.min()
if orientation == 'c':
    x2min = -X2.max()
x2max = X2.max()

x1min = -Rmax
x2min = -Rmax/2.0
x1max = Rmax
x2max = Rmax/2.0

# Variable to plot (rho,pre,vx1,vx2,vel)
if plot_var == 'rho':
    if scale == 'lin':
        plot   = rho
    if scale == 'log':
        plot = np.log10(rho)
if plot_var == 'temp':
    if scale == 'lin':
        plot = (pre/rho)/(pre[Nx2-1,Nx1-1]/rho[Nx2-1,Nx1-1])
    if scale == 'log':
        plot = np.log10((pre/rho)/(pre[Nx2-1,Nx1-1]/rho[Nx2-1,Nx1-1]**(5./3.)))
if plot_var == 'vel':
    if scale == 'lin':
        plot = VV_HP
    if scale == 'log':
        plot = np.log10(VV_NHP)
if plot_var == 'ur':
    if scale == 'lin':
        plot   = np.abs(Uhatr)
    if scale == 'log':
        plot = np.log10(np.abs(Uhatr)+0.00000001)
if plot_var == 'utheta':
    if scale == 'lin':
        plot   = np.abs(Uhattheta)
    if scale == 'log':
        plot = np.log10(np.abs(Uhattheta)+0.00000001)
if plot_var == 'uphi':
    if scale == 'lin':
        plot   = np.abs(Uhatphi)
    if scale == 'log':
        plot = np.log10(np.abs(Uhatphi)+0.00000001)


# Colorbar limits
cbar_ticks = 5
if plot_var == 'rho':
    cbar_min = plot.min() # Min
    cbar_max = plot.max()#cbar_ticks*cbar_min # Max
if plot_var == 'temp':
    cbar_min = -0.0#plot.min()# Min
    cbar_max = 0.5#plot.max()#cbar_ticks*cbar_min # Max
if plot_var == 'vel':
    cbar_min = -2.0 # Min
    cbar_max = 0.0 # Max
if plot_var == 'ur':
    cbar_min = -1.2#plot.min() # Min
    cbar_max = -0.2#plot.max() # Max
if plot_var == 'utheta':
    cbar_min = -5.2#plot.min() # Min
    cbar_max = -1.8#plot.max() # Max
if plot_var == 'uphi':
    cbar_min = plot.min() # Min
    cbar_max = plot.max() # Max


# Contour levels
levels = np.linspace(cbar_min,cbar_max,400) # Levels

# Contour colormap 
if plot_var == 'rho':
    cmap = plt.get_cmap('jet')#jet') # Colormap
if plot_var == 'temp':
    cmap = plt.get_cmap('YlOrRd') # Colormap
if plot_var == 'vel':
    cmap = plt.get_cmap('RdYlBu_r') # Colormap
if plot_var == 'ur':
    cmap = plt.get_cmap('magma') # Colormap
if plot_var == 'utheta':
    cmap = plt.get_cmap('viridis') # Colormap
if plot_var == 'uphi':
    cmap = plt.get_cmap('cividis') # Colormap
# RdYlGn
# YlOrBr
# RdBu
# gist_heat
# afmhot
# jet
# RdGy

#ax.set_xscale("log")
#ax.set_yscale("log")

# Contour normalization
norm = BoundaryNorm(levels, ncolors=cmap.N) # Normalization
if plot_var == 'rho':
    lev = 20
    lev = np.linspace(cbar_min,cbar_max,20)
if plot_var == 'temp':
    lev = 20
    lev = np.linspace(cbar_min,cbar_max,20)
if plot_var == 'vel':
    lev = [0.05,0.1,0.15,0.255555,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.9999]
    lev = 20
if plot_var == 'ur':
    lev = 20
    lev = np.linspace(cbar_min,cbar_max,20)
if plot_var == 'utheta':
    lev = 20
    lev = np.linspace(cbar_min,cbar_max,20)
if plot_var == 'uphi':
    lev = 20
    lev = np.linspace(cbar_min,cbar_max,20)

################
# Contour plot #
################
if orientation == 'v':
    #cn = ax.contour(X1,X2,plot,levels=lev,linewidths=0.3,colors='k',linestyles='dashed')
    if plot_var == 'rho':
        cn = ax.contourf( X1, X2,plot,cmap=cmap,levels=levels,norm=norm,extend='both')
        cn = ax.contourf( X1,-X2,plot,cmap=cmap,levels=levels,norm=norm,extend='both')
        cn = ax.contourf(-X1, X2,plot,cmap=cmap,levels=levels,norm=norm,extend='both')
        cn = ax.contourf(-X1,-X2,plot,cmap=cmap,levels=levels,norm=norm,extend='both')
    if plot_var == 'temp':
        cn = ax.contourf(X1,X2,plot,cmap=cmap,levels=levels,norm=norm,extend='both')
    if plot_var == 'vel':
        cn = ax.contourf(X1,X2,plot,cmap=cmap,levels=levels,norm=norm,extend='min')
    if plot_var == 'ur':
        cn = ax.contourf(X1,X2,plot,cmap=cmap,levels=levels,norm=norm,extend='both')
    if plot_var == 'utheta':
        cn = ax.contourf(X1,X2,plot,cmap=cmap,levels=levels,norm=norm,extend='both')
    if plot_var == 'uphi':
        cn = ax.contourf(X1,X2,plot,cmap=cmap,levels=levels,norm=norm,extend='both')

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
    sx1 = np.linspace(x1min,x1max,Nx1)
    sx2 = np.linspace(x2min,x2max,Nx2)
    sX1, sX2 = np.meshgrid(sx1,sx2)

    px1  = X1.flatten()
    px2  = X2.flatten()
    pvX1 = vX1.flatten()
    pvX2 = vX2.flatten()

    gvX1 = griddata((px1,px2),pvX1,(sX1,sX2))
    gvX2 = griddata((px1,px2),pvX2,(sX1,sX2))

    if orientation == 'v':
        ax.streamplot( sx1,sx2, gvX1,gvX2,density=2,color='white',
                                linewidth=0.5)#9.5*np.sqrt(gvX1*gvX1 + gvX2*gvX2))
        ax.streamplot(-sx1,sx2,-gvX1,gvX2,density=2,color='white',
                                linewidth=0.5)#9.5*np.sqrt(gvX1*gvX1 + gvX2*gvX2))
        ax.streamplot( sx1,-sx2, gvX1,-gvX2,density=2,color='white',
                                linewidth=0.5)#9.5*np.sqrt(gvX1*gvX1 + gvX2*gvX2))
        ax.streamplot(-sx1,-sx2,-gvX1,-gvX2,density=2,color='white',
                                linewidth=0.5)#9.5*np.sqrt(gvX1*gvX1 + gvX2*gvX2))

    if orientation == 'h':
        ax.streamplot(sx2, sx1,gvX2, gvX1,density=[2,1],color='k',linewidth=np.sqrt(gvX1*gvX1 + gvX2*gvX2))
        ax.streamplot(sx2,-sx1,gvX2, gvX1,density=[2,1],color='k',linewidth=np.sqrt(gvX1*gvX1 + gvX2*gvX2))

    if orientation == 'c':
        ax.streamplot( sx1,sx2, gvX1,gvX2,density=1,color='k',
                                linewidth=2.0*np.sqrt(gvX1*gvX1 + gvX2*gvX2))
        ax.streamplot(-sx1,sx2,-gvX1,gvX2,density=1,color='k',
                                linewidth=2.0*np.sqrt(gvX1*gvX1 + gvX2*gvX2))
        ax.streamplot( sx1,-sx2, gvX1,-gvX2,density=1,color='k',
                                linewidth=2.0*np.sqrt(gvX1*gvX1 + gvX2*gvX2))
        ax.streamplot(-sx1,-sx2,-gvX1,-gvX2,density=1,color='k',
                                linewidth=2.0*np.sqrt(gvX1*gvX1 + gvX2*gvX2))

rplus = 1. + np.sqrt(1. - a**2)

ax.contour( X1,X2,r,levels=[rplus],linewidths=0.4,colors='w')
ax.contour(-X1,X2,r,levels=[rplus],linewidths=0.4,colors='w')
ax.contour( X1,-X2,r,levels=[rplus],linewidths=0.4,colors='w')
ax.contour(-X1,-X2,r,levels=[rplus],linewidths=0.4,colors='w')
if orientation == 'c':
    ax.contour( X1,-X2,r,levels=[rplus],linewidths=0.4,colors='w')
    ax.contour(-X1,-X2,r,levels=[rplus],linewidths=0.4,colors='w')

#ax.contour( X1,X2,r,levels=[RR],linewidths=1.0,colors='k')
#ax.contour(-X1,X2,r,levels=[RR],linewidths=1.0,colors='k')
if orientation == 'c':
    ax.contour( X1,-X2,r,levels=[RR],linewidths=1.0,colors='k')
    ax.contour(-X1,-X2,r,levels=[RR],linewidths=1.0,colors='k')

if plot_var == 'rho':
    rr = 1.1+0.0001
if plot_var == 'temp':
    rr = 1.1+0.0001
if plot_var == 'ur':
    rr = 1.1+0.0001
if plot_var == 'utheta':
    rr = 1.1+0.0001
if plot_var == 'uphi':
    rr = 1.1+0.0001
if plot_var == 'vel':
    rr = rplus;

ax.contour( X1,X2,r,levels=[rr],linewidths=1.0,colors='k')
ax.contour(-X1,X2,r,levels=[rr],linewidths=1.0,colors='k')
ax.contour( X1,-X2,r,levels=[rr],linewidths=1.0,colors='k')
ax.contour(-X1,-X2,r,levels=[rr],linewidths=1.0,colors='k')
if orientation == 'c':
    ax.contour( X1,-X2,r,levels=[rr],linewidths=1.0,colors='k')
    ax.contour(-X1,-X2,r,levels=[rr],linewidths=1.0,colors='k')

ax.contour( X1, X2,sonic,levels=[1.0],linewidths=1.0,colors='r')
ax.contour( X1,-X2,sonic,levels=[1.0],linewidths=1.0,colors='r')
ax.contour(-X1, X2,sonic,levels=[1.0],linewidths=1.0,colors='r')
ax.contour(-X1,-X2,sonic,levels=[1.0],linewidths=1.0,colors='r')

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
    plt.ylim(x2min,x2max)

###################
# X1 and X2 ticks #
###################

dens = (x2max-x2min)/(x1max-x1min)
num_tick = 6
num_x1 = 5#int(num_tick)
num_x2 = 5#int(num_tick*dens)
if orientation == 'c':
    num_x2 = num_x1

x1labels = np.linspace(x1min,x1max, num=7, endpoint=True) # num of ticks in X1 axis
x2labels = np.linspace(x2min,x2max, num=7, endpoint=True) # num of ticks in X2 axis

if (orientation == 'v') or (orientation == 'c'):
    plt.xticks(x1labels)
    ax.set_xticklabels(['{:.0f}'.format(x) for x in x1labels],fontsize=fontsize)
    #ax.set_xticklabels([-10,-8,-6,-4,-2,0,2,4,6,8,10],fontsize=fontsize)
    plt.yticks(x2labels)
    ax.set_yticklabels(['{:.0f}'.format(x) for x in x2labels],fontsize=fontsize)
    #ax.set_yticklabels([0,2,4,6,8,10],fontsize=fontsize)
    if COORD == 'CARTESIAN':
        plt.xlabel(r'$x$',fontsize=fontsize)
        plt.ylabel(r'$y$',fontsize=fontsize)
    if (COORD == 'CYLINDRICAL') or (COORD == 'SPHERICAL'):
        plt.xlabel(r'$R/M$',fontsize=fontsize)
        plt.ylabel(r'$z/M$',fontsize=fontsize)

if orientation == 'h':
    plt.xticks(x2labels)
    ax.set_xticklabels(['{:.1f}'.format(x) for x in x2labels],fontsize=fontsize)
    plt.yticks(x1labels)
    ax.set_yticklabels(['{:.1f}'.format(x) for x in x1labels],fontsize=fontsize)
    if COORD == 'CARTESIAN':
        plt.xlabel(r'$x$',fontsize=fontsize)
        plt.ylabel(r'$y$',fontsize=fontsize)
    if (COORD == 'CYLINDRICAL') or (COORD == 'SPHERICAL'):
        plt.xlabel(r'$z$',fontsize=fontsize)
        plt.ylabel(r'$R$',fontsize=fontsize)

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

cax = inset_axes(ax,width='3%',height="100%",loc = 'lower right',bbox_to_anchor = (0.07,0.0,1,1),bbox_transform = ax.transAxes,borderpad = 0)
cbarn = fig.colorbar(cn,orientation=cbor,cax=cax)

if scale == 'lin':
    if plot_var == 'rho':
        cbarn.set_label(r'$\rho/\rho_0$',rotation=rotation,fontsize=fontsize,labelpad=-25,y=-0.05)
    if plot_var == 'pre':
        cbarn.set_label(r'$P/P_0$',rotation=rotation,fontsize=fontsize,labelpad=20)
    if plot_var == 'vel':
        cbarn.set_label(r'$V$',rotation=rotation,fontsize=fontsize,labelpad=20)
    if plot_var == 'ur':
        cbarn.set_label(r'$U^{\hat{\r}}$',rotation=rotation,fontsize=fontsize,labelpad=20)
    if plot_var == 'utheta':
        cbarn.set_label(r'$U^{\hat{\theta}}$',rotation=rotation,fontsize=fontsize,labelpad=20)
    if plot_var == 'uphi':
        cbarn.set_label(r'$U^{\hat{\phi}}$',rotation=rotation,fontsize=fontsize,labelpad=20)
if scale == 'log':
    if plot_var == 'rho':
        cbarn.set_label(r'$\log \left(\rho/\rho_\infty \right)$',rotation=rotation,fontsize=fontsize,labelpad=25)
    if plot_var == 'temp':
        cbarn.set_label(r'$\log \left(\Theta/\Theta_\infty \right)$',rotation=rotation,fontsize=fontsize,labelpad=25)
    if plot_var == 'vel':
        cbarn.set_label(r'$\log \left(V/\mathrm{c} \right)$',rotation=rotation,fontsize=fontsize,labelpad=25)
    if plot_var == 'ur':
        cbarn.set_label(r'$\log \left| U^{\hat{r}}/\mathrm{c} \right|$',rotation=rotation,fontsize=fontsize,labelpad=25)
    if plot_var == 'utheta':
        cbarn.set_label(r'$\log \left| U^{\hat{\theta}}/\mathrm{c} \right|$',rotation=rotation,fontsize=fontsize,labelpad=25)
    if plot_var == 'uphi':
        cbarn.set_label(r'$\log \left| U^{\hat{\phi}}/\mathrm{c} \right|$',rotation=rotation,fontsize=fontsize,labelpad=25)

if (cbpos == "right") or (cbpos == "left"):
   cax.yaxis.set_ticks_position(cbpos) 
   cax.yaxis.set_label_position(cbpos) 
if (cbpos == "top") or (cbpos == "bottom"):
   cax.xaxis.set_ticks_position(cbpos) 
   cax.xaxis.set_label_position(cbpos) 

##################
# Colorbar ticks #
##################
cbax1labels = np.linspace(cbar_min,cbar_max,num=cbar_ticks, endpoint=True) # num of ticks in the colorbar
cbarn.set_ticks(cbax1labels)

if plot_var == 'rho':
    if (cbpos == "right") or (cbpos == "left"):
        cbarn.ax.set_yticklabels(['{:.1f}'.format(x) for x in cbax1labels],fontsize=fontsize)
    if (cbpos == "top") or (cbpos == "bottom"):
        cbarn.ax.set_xticklabels(['{:.1f}'.format(x) for x in cbax1labels],fontsize=fontsize)
if plot_var == 'temp':
    if (cbpos == "right") or (cbpos == "left"):
        cbarn.ax.set_yticklabels(['{:.1f}'.format(x) for x in cbax1labels],fontsize=fontsize)
    if (cbpos == "top") or (cbpos == "bottom"):
        cbarn.ax.set_xticklabels(['{:.1f}'.format(x) for x in cbax1labels],fontsize=fontsize)
if plot_var == 'vel':
    if (cbpos == "right") or (cbpos == "left"):
        cbarn.ax.set_yticklabels(['{:.1f}'.format(x) for x in cbax1labels],fontsize=fontsize)
    if (cbpos == "top") or (cbpos == "bottom"):
        cbarn.ax.set_xticklabels(['{:.1f}'.format(x) for x in cbax1labels],fontsize=fontsize)
if plot_var == 'ur':
    if (cbpos == "right") or (cbpos == "left"):
        cbarn.ax.set_yticklabels(['{:.2f}'.format(x) for x in cbax1labels],fontsize=fontsize)
    if (cbpos == "top") or (cbpos == "bottom"):
        cbarn.ax.set_xticklabels(['{:.2f}'.format(x) for x in cbax1labels],fontsize=fontsize)
if plot_var == 'utheta':
    if (cbpos == "right") or (cbpos == "left"):
        cbarn.ax.set_yticklabels(['{:.2f}'.format(x) for x in cbax1labels],fontsize=fontsize)
    if (cbpos == "top") or (cbpos == "bottom"):
        cbarn.ax.set_xticklabels(['{:.2f}'.format(x) for x in cbax1labels],fontsize=fontsize)
if plot_var == 'uphi':
    if (cbpos == "right") or (cbpos == "left"):
        cbarn.ax.set_yticklabels(['{:.2f}'.format(x) for x in cbax1labels],fontsize=fontsize)
    if (cbpos == "top") or (cbpos == "bottom"):
        cbarn.ax.set_xticklabels(['{:.2f}'.format(x) for x in cbax1labels],fontsize=fontsize)
ax.set_aspect('equal')

ax.set_xlim([x1min, x1max])
ax.set_ylim([x2min, x2max])
if orientation == 'c':
    ax.set_xlim([-x1max - 0.2, x1max + 0.2])
    ax.set_ylim([-x2max - 0.2, x2max + 0.2])

#blackhole = plt.Circle((0,0),2.9,color='k',fill = True,ls='-',lw=0.5,zorder=0)
#ax.add_patch(blackhole)

#Graph name
plt.savefig(ouput_name,dpi=300,bbox_inches="tight")
plt.close()
