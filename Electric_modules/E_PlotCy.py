#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from numpy import zeros, ones, arange, dot, array, linspace, meshgrid, cos, sin, pi

def Cells(q):
    A = []
    for i in range(q):
        A.append([])
    
    return A

def Pol2Cart(rho, phi, z):
    x = rho * cos(phi)
    y = rho * sin(phi)
    return( x, y, z)

#%% Filaments

def draw_line(ax, center, operator, l, I, need_curr_dir):
    
    fil_1  = dot( operator, array([0,0, l/2]).reshape([3,1], order='F') ) + center.reshape([3,1], order='F')
    fil_2  = dot( operator, array([0,0,-l/2]).reshape([3,1], order='F') ) + center.reshape([3,1], order='F')

    x = array([float( fil_1[0]), float( fil_2[0]) ])
    y = array([float( fil_1[1]), float( fil_2[1]) ])
    z = array([float( fil_1[2]), float( fil_2[2]) ])
    
    ax.plot( x, y, z, color = 'fuchsia', lw = 3)
    
    if need_curr_dir == 1:
        if I > 0:
            ax.quiver(center[0], center[1], center[2], float( fil_1[0] - center[0] ), float( fil_1[1] - center[1] ), float( fil_1[2] - center[2] ) , lw = 3, Color = 'yellow', normalize = True)
        elif I < 0:
            ax.quiver(center[0], center[1], center[2], float( fil_2[0] - center[0] ), float( fil_2[1] - center[1] ), float( fil_2[2] - center[2] ) , lw = 3, Color = 'yellow', normalize = True)
   
#%% Torus
    
def Torus_plot(ax, Ro, r, center, Euler):
  
    tht = linspace(0, 2* pi, 25)
    phi = linspace(0, 2* pi, 25)
    
    tht, phi = meshgrid(tht, phi)
    
    x = (Ro + r* cos(tht))* cos(phi)
    y = (Ro + r* cos(tht))* sin(phi)
    z = r* sin(tht)
    
    # Euler rotation
    rot = dot( Euler, array([ x.ravel(), y.ravel(), z.ravel()]) )
    
    x_rot = rot[0,:].reshape( x.shape)
    y_rot = rot[1,:].reshape( y.shape)
    z_rot = rot[2,:].reshape( z.shape)

    # transformation
    x = x_rot + center[0]
    y = y_rot + center[1]
    z = z_rot + center[2]
    
    ax.plot_surface(x, y, z, alpha = 0.2, color = (0.1, 0.5, 1), edgecolors='k')