#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from numpy import zeros, ones, arange, dot, array, linspace, meshgrid, cos, sin

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

#%% Slabs

def draw_slab(ax, mp3d, center, operator, ver_org, J, need_curr_dir):
    
    # transformation of slab
    vertices_trnf = Cells(8)
    for j in range(0,8):
        vertices_trnf[j] = dot(operator, ver_org[j].reshape(3,1) ) + center.reshape([3,1])

    s_verts1 = [( float(vertices_trnf[0][0]), float(vertices_trnf[0][1]), float(vertices_trnf[0][2])),
                ( float(vertices_trnf[1][0]), float(vertices_trnf[1][1]), float(vertices_trnf[1][2])),
                ( float(vertices_trnf[2][0]), float(vertices_trnf[2][1]), float(vertices_trnf[2][2])),
                ( float(vertices_trnf[3][0]), float(vertices_trnf[3][1]), float(vertices_trnf[3][2])) ]

    s_verts2 = [( float(vertices_trnf[0][0]), float(vertices_trnf[0][1]), float(vertices_trnf[0][2])),
                ( float(vertices_trnf[4][0]), float(vertices_trnf[4][1]), float(vertices_trnf[4][2])),
                ( float(vertices_trnf[5][0]), float(vertices_trnf[5][1]), float(vertices_trnf[5][2])),
                ( float(vertices_trnf[1][0]), float(vertices_trnf[1][1]), float(vertices_trnf[1][2])) ]

    s_verts3 = [( float(vertices_trnf[0][0]), float(vertices_trnf[0][1]), float(vertices_trnf[0][2])),
                ( float(vertices_trnf[3][0]), float(vertices_trnf[3][1]), float(vertices_trnf[3][2])),
                ( float(vertices_trnf[7][0]), float(vertices_trnf[7][1]), float(vertices_trnf[7][2])),
                ( float(vertices_trnf[4][0]), float(vertices_trnf[4][1]), float(vertices_trnf[4][2])) ]

    s_verts4 = [( float(vertices_trnf[3][0]), float(vertices_trnf[3][1]), float(vertices_trnf[3][2])),
                ( float(vertices_trnf[7][0]), float(vertices_trnf[7][1]), float(vertices_trnf[7][2])),
                ( float(vertices_trnf[6][0]), float(vertices_trnf[6][1]), float(vertices_trnf[6][2])),
                ( float(vertices_trnf[2][0]), float(vertices_trnf[2][1]), float(vertices_trnf[2][2])) ]

    s_verts5 = [( float(vertices_trnf[1][0]), float(vertices_trnf[1][1]), float(vertices_trnf[1][2])),
                ( float(vertices_trnf[2][0]), float(vertices_trnf[2][1]), float(vertices_trnf[2][2])),
                ( float(vertices_trnf[6][0]), float(vertices_trnf[6][1]), float(vertices_trnf[6][2])),
                ( float(vertices_trnf[5][0]), float(vertices_trnf[5][1]), float(vertices_trnf[5][2])) ]

    s_verts6 = [( float(vertices_trnf[4][0]), float(vertices_trnf[4][1]), float(vertices_trnf[4][2])),
                ( float(vertices_trnf[7][0]), float(vertices_trnf[7][1]), float(vertices_trnf[7][2])),
                ( float(vertices_trnf[6][0]), float(vertices_trnf[6][1]), float(vertices_trnf[6][2])),
                ( float(vertices_trnf[5][0]), float(vertices_trnf[5][1]), float(vertices_trnf[5][2])) ]
    ###
    f1 = mp3d.art3d.Poly3DCollection([ s_verts1 ], alpha=0.5, linewidth=0.7, edgecolors='k' )
    f2 = mp3d.art3d.Poly3DCollection([ s_verts2 ], alpha=0.5, linewidth=0.7, edgecolors='k' )
    f3 = mp3d.art3d.Poly3DCollection([ s_verts3 ], alpha=0.5, linewidth=0.7, edgecolors='k' )
    f4 = mp3d.art3d.Poly3DCollection([ s_verts4 ], alpha=0.5, linewidth=0.7, edgecolors='k' )
    f5 = mp3d.art3d.Poly3DCollection([ s_verts5 ], alpha=0.5, linewidth=0.7, edgecolors='k' )
    f6 = mp3d.art3d.Poly3DCollection([ s_verts6 ], alpha=0.5, linewidth=0.7, edgecolors='k' )
    #
    f1.set_facecolor((0.1, 0.5, 1, 0.5)) # Last argument is alpha = 0.5
    f2.set_facecolor((0.1, 0.5, 1, 0.5))
    f3.set_facecolor((0.1, 0.5, 1, 0.5))
    f4.set_facecolor((0.1, 0.5, 1, 0.5))
    f5.set_facecolor((0.1, 0.5, 1, 0.5))
    f6.set_facecolor((0.1, 0.5, 1, 0.5))
    #
    ax.add_collection3d(f1)
    ax.add_collection3d(f2)
    ax.add_collection3d(f3)
    ax.add_collection3d(f4)
    ax.add_collection3d(f5)
    ax.add_collection3d(f6)
    #
    if need_curr_dir == 1:
        if J > 0:
            ax.quiver( vertices_trnf[0][0], vertices_trnf[0][1], vertices_trnf[0][2], vertices_trnf[4][0] - vertices_trnf[0][0], vertices_trnf[4][1] - vertices_trnf[0][1], vertices_trnf[4][2] - vertices_trnf[0][2], normalize = True, color = 'yellow' )
        elif J < 0:
            ax.quiver( vertices_trnf[4][0], vertices_trnf[4][1], vertices_trnf[4][2], vertices_trnf[0][0] - vertices_trnf[4][0], vertices_trnf[0][1] - vertices_trnf[4][1], vertices_trnf[0][2] - vertices_trnf[4][2], normalize = True, color = 'yellow' )

#%% Arcs

def circle(ax, r, theta_S, theta_E, center, Euler):
    ''' r       ---> radius of the circle
        theta_S ---> starting angle of the circle
        theta_E ---> ending angle of the circle
        center  ---> center of the circle (shape 3*1)
        Euler   ---> Euler rotation matrix of the circle '''
    # Circle at origin
    t = arange( theta_S, theta_E + 0.05, 0.05)
    x, y, z = Pol2Cart( r, t, 0)
    z = zeros( len(x) )
    C = array([x, y, z])

    # Euler rotation of the circle
    cir = dot( Euler, C)
    
    # center shifting
    circle = cir + center
    X = circle[0,:]
    Y = circle[1,:]
    Z = circle[2,:]
    
    ax.plot(X, Y, Z)

#%%
        
# roll to draw the cylindre
def roll(ax, R, zi, zf, ang_i, ang_f, center, Euler):
    
    t = linspace( ang_i, ang_f, 50)
    z = array([zi, zf])
    t, z = meshgrid(t, z)
    p, q = t.shape
    r = R* ones([p,q], float)
    
    # cylindrical coordinates to Cartesian coordinate
    x, y, z = Pol2Cart(r,t,z)

    # Euler rotation
    rot = dot( Euler, array([ x.ravel(), y.ravel(), z.ravel()]) )

    x_rot = rot[0,:].reshape( x.shape)
    y_rot = rot[1,:].reshape( y.shape)
    z_rot = rot[2,:].reshape( z.shape)

    # transformation
    x = x_rot + center[0]
    y = y_rot + center[1]
    z = z_rot + center[2]
    
    ax.plot_surface( x, y, z, alpha = 0.5, color = (0.1, 0.5, 1), linewidth=0.7, edgecolors='k' )

# Disk part in cylindre
def disk(ax, Ri, Rf, ang_i, ang_f, center, Euler, d):
    
    r = array( [ Ri, Rf] )
    t = linspace( ang_i, ang_f, 50)
    t, r = meshgrid(t, r)
    p, q = t.shape
    z = zeros([p,q], float) + d
    
    # cylindrical coordinates to Cartesian coordinate
    x, y, z = Pol2Cart(r,t,z)

    # Euler rotation
    rot = dot( Euler, array([ x.ravel(), y.ravel(), z.ravel()]) )

    x_rot = rot[0,:].reshape( x.shape)
    y_rot = rot[1,:].reshape( y.shape)
    z_rot = rot[2,:].reshape( z.shape)

    # transformation
    x = x_rot + center[0]
    y = y_rot + center[1]
    z = z_rot + center[2]
    
    ax.plot_surface( x, y, z, alpha = 0.5, color = (0.1, 0.5, 1), linewidth=0.7, edgecolors='k' )