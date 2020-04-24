#!/usr/bin/env python3
# -*- coding: utf-8 -*-

print('   =================================================================================== ')
print('                                          CAET.py                                      ')
print('   =================================================================================== ')
print(' ')
print('   ----------------------------------------------------------------------------------- ')
print('   Reading all the inputs...')

FILE = open(r'/home/divyang/Documents/Python/CAE/CAE v1.2/Input_files/CAET.txt', 'r')
Variables = FILE.read()
Variables = Variables.split()

#%%

# mass of electron
m = float( Variables[6] )

# charge of electron
q = float( Variables[9] )

# Initial position
x = float( Variables[14] )
y = float( Variables[16] )
z = float( Variables[18] )

# Initial velocity
Vx = float( Variables[23] )
Vy = float( Variables[25] )
Vz = float( Variables[27] )

# Constant B field
Const_B  = int( Variables[32] )

Const_Bx = float( Variables[34] )
Const_By = float( Variables[36] )
Const_Bz = float( Variables[38] )

# Constant E field
Const_E  = int( Variables[43] )

Const_Ex = float( Variables[45] )
Const_Ey = float( Variables[47] )
Const_Ez = float( Variables[49] )

# Simulation time
t_max  = float( Variables[53] )
Points = float( Variables[55] )

# Tokamak simulation
Cylinder_mode    = int( Variables[59] )
Cylinder_r_inner = float( Variables[62] )
Cylinder_r_outer = float( Variables[65] )
Cylinder_height  = float( Variables[67] )

# Add constant B field
# for example earth's magnetic field
B_x_extra = 0
B_y_extra = 0
B_z_extra = 0

B_extra_mode = int( Variables[73] )

if B_extra_mode == 1:
    B_x_extra = float( Variables[75] )
    B_y_extra = float( Variables[77] )
    B_z_extra = float( Variables[79] )

Phase_velocity = int( Variables[83] )

print('   ----------------------------------------------------------------------------------- ')