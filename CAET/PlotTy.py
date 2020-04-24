#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append('/home/divyang/Documents/Python/CAE/CAE v1.2')

Path = r'/home/divyang/Documents/Python/CAE/CAE v1.2/Outputs/CAET output'

from TInputs import Const_B, Const_E, m, Phase_velocity
from numpy import array, float, pi

f = open( Path + '/x.txt', 'r')
X = array( f.read().split() ).astype(float)
f.close()

f = open( Path + '/y.txt', 'r')
Y = array( f.read().split() ).astype(float)
f.close()

f = open( Path + '/z.txt', 'r')
Z = array( f.read().split() ).astype(float)
f.close()

f  = open( Path + '/vx.txt', 'r')
Vx = array( f.read().split() ).astype(float)
f.close()

f  = open( Path + '/vy.txt', 'r')
Vy = array( f.read().split() ).astype(float)
f.close()

f  = open( Path + '/vz.txt', 'r')
Vz = array( f.read().split() ).astype(float)
f.close()

#%%

Px = m* array(Vx)
Py = m* array(Vy)
Pz = m* array(Vz)

#%% Plotting
from matplotlib import pyplot as plt
import mpl_toolkits.mplot3d as mp3d

#%% Phase space trajectory
if Phase_velocity == 1:
    
    fig = plt.figure('Phase space trajectory - X')
    plt.plot(X, Px)
    plt.title('X position versus Px momentum')
    plt.grid()
    
    fig = plt.figure('Phase space trajectory - Y')
    plt.plot(Y, Py)
    plt.title('Y position versus Py momentum')
    plt.grid()
    
    fig = plt.figure('Phase space trajectory - Z')
    plt.title('Z position versus Pz momentum')
    plt.plot(Z, Pz)
    plt.grid()

#%%
fig = plt.figure('Trajectory')
ax1 = fig.add_subplot( 111 , projection='3d')

ax1.plot(X, Y, Z, c = 'r', linewidth=2.0)
ax1.scatter(X[0], Y[0], Z[0])

ax1.set_xlabel('X', fontsize = 12, fontweight = 'bold')
ax1.set_ylabel('Y', fontsize = 12, fontweight = 'bold')
ax1.set_zlabel('Z', fontsize = 12, fontweight = 'bold')
plt.axis('square')
plt.show()

#%%
fig = plt.figure('Trajectory (scatter)')
ax2 = fig.add_subplot( 111 , projection='3d')

ax2.scatter(X, Y, Z, c = 'r')
ax2.scatter(X[0], Y[0], Z[0])

ax2.set_xlabel('X', fontsize = 12, fontweight = 'bold')
ax2.set_ylabel('Y', fontsize = 12, fontweight = 'bold')
ax2.set_zlabel('Z', fontsize = 12, fontweight = 'bold')

plt.show()

#%% Trajectory with Electro-magnets

if Const_B == 0:
    
    from Magnet_modules.B_Inputs import f_Nmbr, s_Nmbr, a_Nmbr, ra_Nmbr, \
                                        cy_Nmbr, need_curr_dir
    
    fig = plt.figure('Electromagnet')
    ax = fig.add_subplot( 111 , projection='3d')
    
    # Filament
    if f_Nmbr > 0:
        from Magnet_modules.B_Inputs import f_l, f_I, f_operator, f_centre
        from Magnet_modules.PlotCy import draw_line
        
        for i in range(0, f_Nmbr): 
            draw_line(ax, f_centre[i], f_operator[i], f_l[i], f_I[i], need_curr_dir)
    
    # Slab
    if s_Nmbr > 0:
        from Magnet_modules.B_Inputs import s_ver_org, s_J, s_operator, s_centre
        from Magnet_modules.PlotCy import draw_slab
        
        for i in range(0, s_Nmbr):
            draw_slab(ax, mp3d, s_centre[i], s_operator[i], s_ver_org[i], s_J[i], need_curr_dir)

    # Circular Arc
    if a_Nmbr > 0:
        from Magnet_modules.B_Inputs import a_a0, a_phy_1, a_phy_2, a_operator, a_cnt_cart
        from Magnet_modules.PlotCy import circle
        
        for i in range(0, a_Nmbr):
            circle(ax, a_a0[i], a_phy_1[i], a_phy_2[i], a_cnt_cart[i].reshape(3,1), a_operator[i])
    
    # Ractangular arcs
    if ra_Nmbr > 0:
        from Magnet_modules.B_Inputs import ra_r_inner, ra_r_outer, ra_phy_1, \
                                          ra_phy_2, ra_tckns, ra_operator, ra_cnt_cart
        from Magnet_modules.PlotCy import roll, disk
        
        for i in range(0, ra_Nmbr):
            # roll
            roll(ax, ra_r_inner[i], -ra_tckns[i]/2, ra_tckns[i]/2, ra_phy_1[i], ra_phy_2[i], ra_cnt_cart[i], ra_operator[i])
            roll(ax, ra_r_outer[i], -ra_tckns[i]/2, ra_tckns[i]/2, ra_phy_1[i], ra_phy_2[i], ra_cnt_cart[i], ra_operator[i])
            # disk
            disk(ax, ra_r_inner[i], ra_r_outer[i], ra_phy_1[i], ra_phy_2[i], ra_cnt_cart[i], ra_operator[i], ra_tckns[i]/2)
            disk(ax, ra_r_inner[i], ra_r_outer[i], ra_phy_1[i], ra_phy_2[i], ra_cnt_cart[i], ra_operator[i],-ra_tckns[i]/2)

    # Cylindrical wires
    if cy_Nmbr > 0:
        from Magnet_modules.B_Inputs import cy_ro, cy_l, cy_operator, cy_centre
        from Magnet_modules.PlotCy import roll, disk
        
        for i in range(0, cy_Nmbr):
            disk(ax, 0, cy_ro[i], 0, 2* pi, cy_centre[i], cy_operator[i], cy_l[i]/2 )
            disk(ax, 0, cy_ro[i], 0, 2* pi, cy_centre[i], cy_operator[i],-cy_l[i]/2 )
            roll(ax, cy_ro[i],-cy_l[i]/2, cy_l[i]/2, 0, 2* pi, cy_centre[i], cy_operator[i])

    #%%
    from Electric_modules.E_Inputs import torus_Nmbr

    # Electric Torus
    if torus_Nmbr > 0:
        from Electric_modules.E_Inputs import torus_Ro, torus_r, torus_centre, torus_operator
        from Electric_modules.E_PlotCy import Torus_plot
        
        for i in range(0, torus_Nmbr):
            Torus_plot(ax, torus_Ro[i], torus_r[i], torus_centre[i], torus_operator[i])
    
    ax.plot(X, Y, Z, c = 'r', linewidth = 2)
    ax.scatter(X[0], Y[0], Z[0])
    
    ax.set_xlabel('X', fontsize = 12, fontweight = 'bold')
    ax.set_ylabel('Y', fontsize = 12, fontweight = 'bold')
    ax.set_zlabel('Z', fontsize = 12, fontweight = 'bold')
    plt.axis('square')
    plt.show()
    
print('   -----------------------------------------------------------------------------------   ')