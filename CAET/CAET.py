#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
start_time = time.time()

from numpy.linalg import inv
from numpy import pi, array, dot, savetxt, linspace
from tqdm import tqdm

π = pi

import sys
sys.path.append('/home/divyang/Documents/Python/CAE/CAE v1.2')

Path = r'/home/divyang/Documents/Python/CAE/CAE v1.2/Outputs/CAET output'

from TInputs import m, q, x, y, z, Vx, Vy, Vz, Const_B, Const_Bx, Const_By, \
                    Const_Bz, Const_E, Const_Ex, Const_Ey, Const_Ez, t_max, \
                    Points, B_x_extra, B_y_extra, B_z_extra

#%% Magnet Modules
 
if Const_B == 0:
    from Magnet_modules.B_Inputs import f_Nmbr, s_Nmbr, a_Nmbr, ra_Nmbr, cy_Nmbr, µ
    
    if f_Nmbr > 0:
        from Magnet_modules.B_Inputs import f_l, f_I, f_operator, f_center
        from Magnet_modules.Magnets import B_f
    
    if s_Nmbr > 0:
        from Magnet_modules.B_Inputs import s_ver_org, s_J, s_operator, s_center
        from Magnet_modules.Magnets import B_s
    
    if a_Nmbr > 0:
        from Magnet_modules.B_Inputs import a_a0, a_phy_1, a_phy_2, a_I, a_operator, a_cnt_cart
        from Magnet_modules.Magnets import B_a
    
    if ra_Nmbr > 0:
        from Magnet_modules.B_Inputs import ra_r_inner, ra_r_outer, ra_phy_1, ra_phy_2, \
                                            ra_tckns, ra_J, ra_operator, ra_cnt_cart
        from Magnet_modules.Magnets import B_ra
    
    if cy_Nmbr > 0:
        from Magnet_modules.B_Inputs import cy_ro, cy_l, cy_J, cy_operator, cy_center
        from Magnet_modules.Magnets import B_cy


def Magnet(x, y, z):
    
    point = array([x, y, z])
    
    # dummy Magnetic field variables
    f_B = 0
    s_B = 0
    a_B = 0
    ra_B = 0
    cy_B = 0
    
    if f_Nmbr > 0:
        f_B = B_f(µ, f_Nmbr, f_l, f_I, f_operator, f_center, point)    
    
    if s_Nmbr > 0:
        s_B = B_s(µ, s_Nmbr, s_ver_org, s_J, s_operator, s_center, point)
    
    if a_Nmbr > 0:
        a_B = B_a(µ, a_Nmbr, a_a0, a_phy_1, a_phy_2, a_I, a_operator, a_cnt_cart, point)
    
    if ra_Nmbr > 0:
        ra_B = B_ra(µ, ra_Nmbr, ra_r_inner, ra_r_outer, ra_phy_1, ra_phy_2, \
                    ra_tckns, ra_J, ra_operator, ra_cnt_cart, point)
    
    if cy_Nmbr > 0:
        cy_B = B_cy(µ, cy_Nmbr, cy_ro, cy_l, cy_J, cy_operator, cy_center, point)
    
    # Total mangetic field
    B = f_B + s_B + a_B + ra_B + cy_B
    
    return B

#%% Electric field

if Const_E == 0:
    from Electric_modules.E_Inputs import fil_Nmbr, torus_Nmbr, ε
        
    if fil_Nmbr > 0:
        from Electric_modules.E_Inputs import f_l as E_f_l
        from Electric_modules.E_Inputs import f_lmbd as E_f_lmbd
        from Electric_modules.E_Inputs import f_operator as E_f_operator
        from Electric_modules.E_Inputs import f_center as E_f_center
        from Electric_modules.Electro import Electric_f
    
    if torus_Nmbr > 0:
        from Electric_modules.E_Inputs import torus_Ro, torus_r, torus_sgma, torus_operator, torus_center
        from Electric_modules.Electro import Electric_torus

def Electric(x, y, z):
    
    point = array([x, y, z])
    
    # dummy Magnetic field variables
    f_E  = 0
    tr_E = 0
    
    if fil_Nmbr > 0:
        f_E  = Electric_f(ε, fil_Nmbr, E_f_l, E_f_lmbd, E_f_operator, E_f_center, point)    
    
    if torus_Nmbr > 0:
        tr_E = Electric_torus(ε, torus_Nmbr, torus_Ro, torus_r, torus_sgma, torus_operator, torus_center, point)
    
    # Total electric field
    E = f_E + tr_E
    
    return E


#%% Files

from FERSA import Rot_xy, FERSA

x_data  = open( Path + '/x.txt', 'w')
y_data  = open( Path + '/y.txt', 'w')
z_data  = open( Path + '/z.txt', 'w')

Vx_data = open( Path + '/vx.txt', 'w')
Vy_data = open( Path + '/vy.txt', 'w')
Vz_data = open( Path + '/vz.txt', 'w')

Ex_data = open( Path + '/Ex.txt', 'w')
Ey_data = open( Path + '/Ey.txt', 'w')
Ez_data = open( Path + '/Ez.txt', 'w')

Bx_data = open( Path + '/Bx.txt', 'w')
By_data = open( Path + '/By.txt', 'w')
Bz_data = open( Path + '/Bz.txt', 'w')

ro = array([  x,  y,  z ])
Vo = array([ Vx, Vy, Vz ])

#%% Initial values

# Position
savetxt( x_data, array([ ro[0] ]) )
savetxt( y_data, array([ ro[1] ]) )
savetxt( z_data, array([ ro[2] ]) )

# Velocity
savetxt( Vx_data, array([ Vo[0] ]) )
savetxt( Vy_data, array([ Vo[1] ]) )
savetxt( Vz_data, array([ Vo[2] ]) )

# Electric field
if Const_E == 0:
    E = Electric(ro[0], ro[1], ro[2])

elif Const_E == 1:
    E = array([ Const_Ex, Const_Ey, Const_Ez ])

savetxt( Ex_data, array([ E[0] ]) )
savetxt( Ey_data, array([ E[1] ]) )
savetxt( Ez_data, array([ E[2] ]) )

# Magnetic field
if Const_B == 0:
    B = Magnet(ro[0], ro[1], ro[2])
    
    # Add the extra field value
    B = B + array([ B_x_extra, B_y_extra, B_z_extra ])

elif Const_B == 1:
    B = array([ Const_Bx, Const_By, Const_Bz ])

savetxt( Bx_data, array([ B[0] ]) )
savetxt( By_data, array([ B[1] ]) )
savetxt( Bz_data, array([ B[2] ]) )

#%% Solver

Time  = linspace(0, t_max, Points)
dt    = Time[1]
count = 0

Ang = []

for t in tqdm(Time):
    
    #%% Magnetic field
    if Const_B == 0:
        B = Magnet(ro[0], ro[1], ro[2])
        B = B + array([ B_x_extra, B_y_extra, B_z_extra ])
    
    elif Const_B == 1:
        B = array([ Const_Bx, Const_By, Const_Bz ])
    
    #%% Electric field
    if Const_E == 0:
        E = Electric(ro[0], ro[1], ro[2])

    elif Const_E == 1:
        E = array([ Const_Ex, Const_Ey, Const_Ez ])
    
    #%% FERSA
    
    # Particle axes
    R, Bz = Rot_xy(B)
    E     = dot(R,  E)
    Vo    = dot(R, Vo)
    
    temp_ro = ro
    
    # New position and velocity in particle axes
    ro, Vo, theta = FERSA(m, q, temp_ro, Vo, E, Bz[2], dt)
    
    # Back to the reality
    ro = temp_ro + dot( inv(R), ro - temp_ro)
    Vo = dot( inv(R), Vo)
    B  = dot( inv(R), Bz)
    E  = dot( inv(R),  E)
    
    Ang.append(theta)
    
    # Save data
    if not( count % 1 ) and (t != Time[-1]):
        savetxt( x_data, array([ ro[0] ]) )
        savetxt( y_data, array([ ro[1] ]) )
        savetxt( z_data, array([ ro[2] ]) )
        
        savetxt( Vx_data, array([ Vo[0] ]) )
        savetxt( Vy_data, array([ Vo[1] ]) )
        savetxt( Vz_data, array([ Vo[2] ]) )
        
        savetxt( Ex_data, array([ E[0] ]) )
        savetxt( Ey_data, array([ E[1] ]) )
        savetxt( Ez_data, array([ E[2] ]) )
        
        savetxt( Bx_data, array([ B[0] ]) )
        savetxt( By_data, array([ B[1] ]) )
        savetxt( Bz_data, array([ B[2] ]) )

    count = count + 1

#%% Close all the files

x_data.close()
y_data.close()
z_data.close()

Vx_data.close()
Vy_data.close()
Vz_data.close()

Ex_data.close()
Ey_data.close()
Ez_data.close()

Bx_data.close()
By_data.close()
Bz_data.close()

#%%
elapsed_time = time.time() - start_time
print('   Calculation completed.')
print('   Elapsed times is ' + str(elapsed_time) + ' s.' )
print('   -----------------------------------------------------------------------------------   ')