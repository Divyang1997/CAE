#!/usr/bin/env python3
# -*- coding: utf-8 -*-

''' CAE   - A Code for Arbitrary ElectroMagnets
    By    : Divyang R. Prajapati
    Guide : Gattu Ramesh Babu '''
    
print('   ===================================================================================   ')
print('                                         CAE v1.2                                        ')
print('   ===================================================================================   ')
print(' ')

#%%
import time
start_time = time.time()

import pandas as pd
from tqdm import tqdm
from numpy import pi, sin, cos, zeros, nan, empty, array, transpose, append, isnan, sqrt

from B_Inputs import f_Nmbr, s_Nmbr, a_Nmbr, ra_Nmbr, cy_Nmbr, Hl_Nmbr, Cart2Pol, \
                     need_plot, need_contour_plot, µ, point, need_curr_dir, \
                     fied_length_quiver, s, Cells

#%%
# dummy Magnetic field variables
f_B  = 0
s_B  = 0
a_B  = 0
ra_B = 0
cy_B = 0
Hl_B = 0

def Field( Short_Name, µ, Nmbr, Func, *args):
    print('   -----------------------------------------------------------------------------------   ')
    print('                                           ' + Short_Name)
    print(' ')
    
    B = empty((3,0), float)
    for j in tqdm(range(0, s)):
        b = Func(µ, Nmbr, *args, point[:,j])
        b = transpose( array([b]) )
        B = append(B, b, axis = 1)
    
    print('   Calculation completed.')
    print('   -----------------------------------------------------------------------------------   ')
    
    return B
    

if f_Nmbr > 0:
    from B_Inputs import f_l, f_I, f_operator, f_center
    print('   Successfully received filament data from Input.py')
    
    from Magnets import B_f
    print('   Calculating filaments...')
  
    f_B = Field('f', µ, f_Nmbr, B_f, f_l, f_I, f_operator, f_center)

if s_Nmbr > 0:
    from B_Inputs import s_ver_org, s_J, s_operator, s_center
    print('   Successfully received slab data from Input.py')
    
    from Magnets import B_s
    print('   Calculating Slabs...')
    
    s_B = Field('s', µ, s_Nmbr, B_s, s_ver_org, s_J, s_operator, s_center)

if a_Nmbr > 0:
    from B_Inputs import a_a0, a_phy_1, a_phy_2, a_I, a_operator, a_cnt_cart
    print('   Successfully received circular arcs data from Input.py')
    
    from Magnets import B_a
    print('   Calculating Circular arcs...')
    
    a_B = Field('a', µ, a_Nmbr, B_a, a_a0, a_phy_1, a_phy_2, a_I, a_operator, a_cnt_cart)

if ra_Nmbr > 0:
    from B_Inputs import ra_r_inner, ra_r_outer, ra_phy_1, ra_phy_2, ra_tckns, ra_J, ra_operator, ra_cnt_cart
    print('   Successfully received circalr arc of rectangular cross-section data from Input.py')
    
    from Magnets import B_ra
    print('   Calculating Circular arcs of rectangular cross-section...')
    
    ra_B = Field('ra', µ, ra_Nmbr, B_ra, ra_r_inner, ra_r_outer, ra_phy_1, \
                  ra_phy_2, ra_tckns, ra_J, ra_operator, ra_cnt_cart)

if cy_Nmbr > 0:
    from B_Inputs import cy_Nmbr, cy_ro, cy_l, cy_J, cy_operator, cy_center
    print('   Successfully received cylindrical wire data from Input.py')
    
    from Magnets import B_cy
    print('   Calculating Cylindrical wires...')
    
    cy_B = Field('cy', µ, cy_Nmbr, B_cy, cy_ro, cy_l, cy_J, cy_operator, cy_center)

if Hl_Nmbr > 0:
    from B_Inputs import Hl_a, Hl_L, Hl_p, Hl_I, Hl_operator, Hl_center
    print('   Successfully received helical wire data from Input.py')
    
    from Magnets import B_Helical
    print('   Calculating Helical wires...')
    
    Hl_B = Field('Hl', µ, Hl_Nmbr, B_Helical, Hl_a, Hl_L, Hl_p, Hl_I, Hl_operator, Hl_center)

    
#%% Total mangetic field

B = f_B + s_B + a_B + ra_B + cy_B + Hl_B

#%%

# cylindrical co-ordinates 
teta_ppol = zeros([1,s])
r_ppol    = zeros([1,s])
z_ppol    = zeros([1,s])
for i in range(0,s):
    r_ppol[0,i], teta_ppol[0,i], z_ppol[0,i] = Cart2Pol(point[0,i], point[1,i], point[2,i])

# cylindrical magnetic field
B_pol = zeros([3,s])
B_pol[0,:] = B[0,:]* cos(teta_ppol) + B[1,:]* sin(teta_ppol)
B_pol[1,:] =-B[0,:]* sin(teta_ppol) + B[1,:]* cos(teta_ppol)
B_pol[2,:] = B[2,:]
#
B_MAGNITUDE = (B[0,:]**2 + B[1,:]**2 + B[2,:]**2)**(0.5)

#%% Output file :-
print('   Generating Output file...')
output_data = pd.DataFrame({'x':( ['(m)'] + list(point[0,:]) ),
                            'y':( ['(m)'] + list(point[1,:]) ),
                            'z':( ['(m)'] + list(point[2,:]) ),
                            'empty1':nan,
                            'r':     ( ['(m)'] + list(r_ppol[0,:]) ),
                            'phi': ( ['(deg)'] + list(teta_ppol[0,:] *(180/pi)) ),
                            'zp':    ( ['(m)'] + list(point[2,:]) ),
                            'empty2':nan,
                            'Bx': ( ['(T)'] + list(B[0,:]) ),
                            'By': ( ['(T)'] + list(B[1,:]) ),
                            'Bz': ( ['(T)'] + list(B[2,:]) ),
                            'empty3':nan,
                            'Br':   ( ['(T)'] + list(B_pol[0,:]) ),
                            'Bphi': ( ['(T)'] + list(B_pol[1,:]) ),
                            'Bzp':  ( ['(T)'] + list(B_pol[2,:]) ),
                            'empty4':nan,
                            'B': ( ['(T)'] + list(B_MAGNITUDE) )}, 
                            columns = ['x','y','z','empty1','r','phi','z','empty2','Bx','By','Bz','empty3','Br','Bphi','Bz','empty4','B'] )

output_data = output_data.rename(columns={'empty1':'', 'empty2':'','empty3':'', 'empty4':''})
output_data.to_excel(r'/home/divyang/Documents/Python/CAE/CAE v1.2/Outputs/MAG_DATA.xlsx',index= False, sheet_name='Mag_field')

print('   Done.')
print('   Output file MAG_DATA.xlsx has been generated in the Outputs folder.')
print('   -----------------------------------------------------------------------------------   ')

#%% Electro-magnet

if need_plot == 1:
    
    print('   Drawing Electro-Magnets and field vectors...')
    
    import matplotlib.pyplot as plt
    import mpl_toolkits.mplot3d as mp3d
    
    fig = plt.figure('Electromagnet')
    ax = fig.add_subplot( 111 , projection='3d')
    
    # Filament
    if f_Nmbr > 0:
        from PlotCy import draw_line
        
        for i in range(0, f_Nmbr): 
            draw_line(ax, f_center[i], f_operator[i], f_l[i], f_I[i], need_curr_dir)
    
    # Slab
    if s_Nmbr > 0:
        from PlotCy import draw_slab
        
        for i in range(0, s_Nmbr):
            draw_slab(ax, mp3d, s_center[i], s_operator[i], s_ver_org[i], s_J[i], need_curr_dir)

    # Circular Arc
    if a_Nmbr > 0:
        from PlotCy import circle
        
        for i in range(0, a_Nmbr):
            circle(ax, a_a0[i], a_phy_1[i], a_phy_2[i], a_cnt_cart[i].reshape(3,1), a_operator[i])
    
    # Ractangular arcs
    if ra_Nmbr > 0:
        from PlotCy import roll, disk
        
        for i in range(0, ra_Nmbr):
            # roll
            roll(ax, ra_r_inner[i], -ra_tckns[i]/2, ra_tckns[i]/2, ra_phy_1[i], ra_phy_2[i], ra_cnt_cart[i], ra_operator[i])
            roll(ax, ra_r_outer[i], -ra_tckns[i]/2, ra_tckns[i]/2, ra_phy_1[i], ra_phy_2[i], ra_cnt_cart[i], ra_operator[i])
            # disk
            disk(ax, ra_r_inner[i], ra_r_outer[i], ra_phy_1[i], ra_phy_2[i], ra_cnt_cart[i], ra_operator[i], ra_tckns[i]/2)
            disk(ax, ra_r_inner[i], ra_r_outer[i], ra_phy_1[i], ra_phy_2[i], ra_cnt_cart[i], ra_operator[i],-ra_tckns[i]/2)

    # Cylindrical wires
    if cy_Nmbr > 0:
        from PlotCy import roll, disk
        
        for i in range(0, cy_Nmbr):
            disk(ax, 0, cy_ro[i], 0, 2* pi, cy_center[i], cy_operator[i], cy_l[i]/2 )
            disk(ax, 0, cy_ro[i], 0, 2* pi, cy_center[i], cy_operator[i],-cy_l[i]/2 )
            roll(ax, cy_ro[i],-cy_l[i]/2, cy_l[i]/2, 0, 2* pi, cy_center[i], cy_operator[i])


    # Sketch the field vectors
    length_B_vec = Cells(s)
    for j in range(0, s):
        if ( isnan( B[0,j]) or isnan( B[1,j]) or isnan( B[1,j]) ):
            length_B_vec[j] = 1
        else:
            length_B_vec[j] = sqrt( B[0,j]**2 + B[1,j]**2 + B[2,j]**2 )
    
    length_B_vec_maxx = max(length_B_vec)
    
    for j in range(0,s):
        ax.quiver(point[0,j], point[1,j], point[2,j], B[0,j], B[1,j], B[2,j], length = fied_length_quiver*length_B_vec[j] /length_B_vec_maxx, normalize = True, color = 'forestgreen')
    
    ax.set_xlabel('X (m)', fontsize = 12, fontweight = 'bold')
    ax.set_ylabel('Y (m)', fontsize = 12, fontweight = 'bold')
    ax.set_zlabel('Z (m)', fontsize = 12, fontweight = 'bold')
    plt.show()
    
    print('   Done.')
    print('   -----------------------------------------------------------------------------------   ')

#%% contour plot

if need_contour_plot == 1:

    print('   Generating contours...')
    
    from B_Inputs import cont_x, cont_y, cont_pp, cont_qq
    
    import matplotlib.pyplot as plt
    from matplotlib import cm
    
    plt.figure('Contour Plot |B|')
    cs = plt.contour( cont_x, cont_y, B_MAGNITUDE.reshape([cont_pp,cont_qq]), 100, cmap = cm.jet)
    plt.clabel(cs, inline = 1, fontsize = 10)
    plt.autoscale(enable = True, axis = 'x', tight = True)
    plt.xlabel('x (m)', fontsize = 20, fontweight = 'bold')
    plt.ylabel('y (m)', fontsize = 20, fontweight = 'bold')
    plt.colorbar()
    plt.axis('square')
    plt.grid(True)
    plt.show()
    
    plt.figure('Contour Plot (filled) |B|')
    cs = plt.contourf( cont_x, cont_y, B_MAGNITUDE.reshape([cont_pp,cont_qq]), 100, cmap = cm.jet)
    plt.autoscale(enable = True, axis = 'x', tight = True)
    plt.xlabel('x (m)', fontsize = 20, fontweight = 'bold')
    plt.ylabel('y (m)', fontsize = 20, fontweight = 'bold')
    plt.colorbar()
    plt.axis('square')
    plt.grid(False)
    plt.show()
    
    print('   Done.')
    print('   -----------------------------------------------------------------------------------   ')

#%%
elapsed_time = time.time() - start_time

print('   ')
print('   Process completed.')
print('   Elapsed times is ' + str(elapsed_time) + ' s.' )

print('   ===================================================================================   ')
print('                            Code Run completed successfully                              ')
print('   ===================================================================================   ')
print('   ')