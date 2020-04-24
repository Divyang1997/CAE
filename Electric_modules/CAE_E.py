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

from E_Inputs import fil_Nmbr, torus_Nmbr, Cart2Pol, \
                     need_plot, need_contour_plot, ε, point, need_curr_dir, \
                     fied_length_quiver, s, Cells

#%%
# dummy Magnetic field variables
f_E  = 0
tr_E = 0

def Field( Short_Name, ε, Nmbr, Func, *args):
    print('   -----------------------------------------------------------------------------------   ')
    print('                                           ' + Short_Name)
    print(' ')
    
    E = empty((3,0), float)
    for j in tqdm(range(0, s)):
        e = Func(ε, Nmbr, *args, point[:,j])
        e = transpose( array([e]) )
        E = append(E, e, axis = 1)
    
    print('   Calculation completed.')
    print('   -----------------------------------------------------------------------------------   ')
    
    return E

if fil_Nmbr > 0:
    from E_Inputs import f_l, f_lmbd, f_operator, f_center
    print('   Successfully received filament data from Input.py')
    
    from Electro import Electric_f
    print('   Calculating Filaments...')
    
    f_E =  Field('f', ε, fil_Nmbr, Electric_f, f_l, f_lmbd, f_operator, f_center)

if torus_Nmbr > 0:
    from E_Inputs import torus_Nmbr, torus_Ro, torus_r, torus_sgma, torus_operator, torus_center
    print('   Successfully received torus data from Input.py')
    
    from Electro import Electric_torus
    print('   Calculating Torus...')
    
    tr_E = Field('Tr', ε, torus_Nmbr, Electric_torus, torus_Ro, \
                 torus_r, torus_sgma, torus_operator, torus_center)

#%% Total mangetic field

E = f_E + tr_E

#%%

# cylindrical co-ordinates
r_ppol    = zeros([1,s])
teta_ppol = zeros([1,s])
z_ppol    = zeros([1,s])
for i in range(0,s):
    r_ppol[0,i], teta_ppol[0,i], z_ppol[0,i] = Cart2Pol(point[0,i], point[1,i], point[2,i])

# cylindrical magnetic field
E_pol = zeros([3,s])
E_pol[0,:] = E[0,:]* cos(teta_ppol) + E[1,:]* sin(teta_ppol)
E_pol[1,:] =-E[0,:]* sin(teta_ppol) + E[1,:]* cos(teta_ppol)
E_pol[2,:] = E[2,:]
#
E_MAGNITUDE = (E[0,:]**2 + E[1,:]**2 + E[2,:]**2)**(0.5)

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
                            'Ex': ( ['(V/m)'] + list(E[0,:]) ),
                            'Ey': ( ['(V/m)'] + list(E[1,:]) ),
                            'Ez': ( ['(V/m)'] + list(E[2,:]) ),
                            'empty3':nan,
                            'Er':   ( ['(V/m)'] + list(E_pol[0,:]) ),
                            'Ephi': ( ['(V/m)'] + list(E_pol[1,:]) ),
                            'Ezp':  ( ['(V/m)'] + list(E_pol[2,:]) ),
                            'empty4':nan,
                            'E': ( ['(V/m)'] + list(E_MAGNITUDE) )}, 
                            columns = ['x','y','z','empty1','r','phi','z','empty2','Ex','Ey','Ez','empty3','Er','Ephi','Ez','empty4','E'] )

output_data = output_data.rename(columns={'empty1':'', 'empty2':'','empty3':'', 'empty4':''})
output_data.to_excel(r'/home/divyang/Documents/Python/CAE/CAE v1.2/Outputs/ELC_DATA.xlsx',index= False, sheet_name='Elc_field')

print('   Done.')
print('   Output file ELC_DATA.xlsx is generated in the Outputs folder.')
print('   -----------------------------------------------------------------------------------   ')

#%% Electro-magnet

if need_plot == 1:
    
    print('   Drawing Electro-Magnets and field vectors...')
    
    import matplotlib.pyplot as plt
    import mpl_toolkits.mplot3d as mp3d
    
    fig = plt.figure('Electromagnet')
    ax = fig.add_subplot( 111 , projection='3d')
    
    # Filament
    if fil_Nmbr > 0:
        from E_PlotCy import draw_line
        
        for i in range(0, fil_Nmbr): 
            draw_line(ax, f_center[i], f_operator[i], f_l[i], f_lmbd[i], need_curr_dir)
    
    if torus_Nmbr > 0:
        from E_PlotCy import Torus_plot
        
        for i in range(0, torus_Nmbr):
            Torus_plot(ax, torus_Ro[i], torus_r[i], torus_center[i], torus_operator[i])

    # Sketch the field vectors
    length_E_vec = Cells(s)
    for j in range(0, s):
        if ( isnan( E[0,j]) or isnan( E[1,j]) or isnan( E[1,j]) ):
            length_E_vec[j] = 1
        else:
            length_E_vec[j] = sqrt( E[0,j]**2 + E[1,j]**2 + E[2,j]**2 )
    
    length_E_vec_maxx = max(length_E_vec)
    
    for j in range(0,s):
        ax.quiver(point[0,j], point[1,j], point[2,j], E[0,j], E[1,j], E[2,j], length = fied_length_quiver*length_E_vec[j] /length_E_vec_maxx, normalize = True, color = 'forestgreen')
    
    ax.set_xlabel('X (m)', fontsize = 12, fontweight = 'bold')
    ax.set_ylabel('Y (m)', fontsize = 12, fontweight = 'bold')
    ax.set_zlabel('Z (m)', fontsize = 12, fontweight = 'bold')
    plt.show()
    
    print('   Done.')
    print('   -----------------------------------------------------------------------------------   ')

#%% contour plot

if need_contour_plot == 1:

    print('   Generating contours...')
    
    from E_Inputs import cont_x, cont_y, cont_pp, cont_qq
    
    import matplotlib.pyplot as plt
    from matplotlib import cm
    
    plt.figure('Contour Plot |E|')
    cs = plt.contour( cont_x, cont_y, E_MAGNITUDE.reshape([cont_pp,cont_qq]), 100, cmap = cm.jet)
    plt.clabel(cs, inline = 1, fontsize = 10)
    plt.autoscale(enable = True, axis = 'x', tight = True)
    plt.xlabel('x (m)', fontsize = 20, fontweight = 'bold')
    plt.ylabel('y (m)', fontsize = 20, fontweight = 'bold')
    plt.colorbar()
    plt.axis('square')
    plt.grid(True)
    plt.show()
    
    plt.figure('Contour Plot (filled) |E|')
    cs = plt.contourf( cont_x, cont_y, E_MAGNITUDE.reshape([cont_pp,cont_qq]), 100, cmap = cm.jet)
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