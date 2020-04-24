#!/usr/bin/env python3
# -*- coding: utf-8 -*-

print('   -----------------------------------------------------------------------------------   ')
print('                                        E_Inputs.py                                      ')
print(' ')
print('   Importing data from Input files...')

from numpy import seterr, sin, cos, sqrt, pi, array, linspace, meshgrid, arctan2, deg2rad
import sys
import xlrd

seterr( divide = 'ignore')
seterr( invalid= 'ignore')

#%%
# IsFloat will check that the data is floating point or not
def IsFloat(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
    
def MatAll(row_number, Nmbr_inCAE):
    
    ''' MatAll will replace data(x,:) operator
        this funtion is only for CAE_1E.xlsx not for CAE_2E.xlsx
        row_number is the row address of cell in CAE_1E.xlsx file
        start is the starting point which is 4 '''
    
    stop = Nmbr_inCAE + 3
    start = 4
    AA = []
    for i in range(start - 1, stop):
        AA.append( data.cell_value( row_number - 1, i))
    
    return AA

def Cells(q):
    ''' Cells will creat the vector of empty cells
        e.g. [[],[],[], ... ] '''
    
    A = []
    for i in range(q):
        A.append([])
    
    return A

# Cart2Pol will convert cartesian points into the polar points
def Cart2Pol(x, y, z):
    ρ = sqrt(x**2 + y**2)
    φ = arctan2(y, x)
    return ρ, φ, z

# Pol2Cart will convert polar points into the cartesian points
def Pol2Cart(ρ, φ, z):
    x = ρ* cos(φ)
    y = ρ* sin(φ)
    return x, y, z

def Euler(θ, φ, ψ):
    return array([[ cos(θ)* cos(φ), sin(ψ)* sin(θ)* cos(φ) - cos(ψ)* sin(φ), cos(ψ)* sin(θ)* cos(φ) + sin(ψ)* sin(φ)],
                  [ cos(θ)* sin(φ), sin(ψ)* sin(θ)* sin(φ) + cos(ψ)* cos(φ), cos(ψ)* sin(θ)* sin(φ) - sin(ψ)* cos(φ)],
                  [        -sin(θ),                          sin(ψ)* cos(θ),                          cos(ψ)* cos(θ)] ])

def Mr_Exit(*args):
    print('\033[1;31;40m ERROR!!! \033[0m : ' + ''.join(args) )
    sys.exit(1)

#%%
Path = r'/home/divyang/Documents/Python/CAE/CAE v1.2/Input_files'
    
FILE_CAE_1E = xlrd.open_workbook( Path + '/CAE_1E.xlsx')
FILE_CAE_2E = xlrd.open_workbook( Path + '/CAE_2E.xlsx')
data = FILE_CAE_1E.sheet_by_index(0)
data_2 = FILE_CAE_2E.sheet_by_index(0)

#
CAET_run = data_2.cell_value(7, 14)
if (type(CAET_run) == str):
    Mr_Exit('Please check the CAET option, 8 row, 15 column in CAE_2E.xlsx')
elif (CAET_run != 1) and (CAET_run != 0):
    Mr_Exit('Please check the CAET option, 8 row, 15 column in CAE_2E.xlsx')

εr = data_2.cell_value(19,6)
if (type(εr) == str):
    Mr_Exit('Please check the permitivity, 20 row, 7 column in CAE_2E.xlsx')

ε0 = 8.854187817e-12
ε  = εr * ε0

need_plot = data_2.cell_value(21,6)
if (type(need_plot) == str):
    Mr_Exit('Please check the plotting option, 22 row, 7 column in CAE_2E.xlsx')
elif (need_plot != 1) and (need_plot != 0):
    Mr_Exit('Please check the plotting option, 22 row, 7 column in CAE_2E.xlsx')

need_curr_dir = data_2.cell_value(25,6)
if (type(need_curr_dir) == str):
    Mr_Exit('Please check the current directon option, 26 row, 7 column in CAE_2E.xlsx')
elif (need_curr_dir != 1) and (need_curr_dir != 0):
    Mr_Exit('Please check the current directon option, 26 row, 7 column in CAE_2E.xlsx')

fied_length_quiver = data_2.cell_value(27,6)
if (type(fied_length_quiver) == str):
    Mr_Exit('Please check the field length option, 28 row, 7 column in CAE_2E.xlsx')


#%% filament

if (type(data.cell_value(2, 3)) == str):
    Mr_Exit('Please check the filament number option, 3 row, 4 column in CAE_1E.xlsx')
    
fil_Nmbr = int(data.cell_value(2, 3))

if fil_Nmbr > 0:
    print('   Importing Filaments...')
    
    f_cnt = array([MatAll(9, fil_Nmbr), MatAll(10, fil_Nmbr), MatAll(11, fil_Nmbr)])
    f_l = MatAll(13, fil_Nmbr)
    # euler transformation angles
    f_ppye = MatAll(17, fil_Nmbr)
    f_ttta = MatAll(18, fil_Nmbr)
    f_pphi = MatAll(19, fil_Nmbr)
    # current in filament
    f_lmbd = MatAll(21, fil_Nmbr)
    
    if not(isinstance(f_cnt[0, 0], float)):
        for i in range(0, fil_Nmbr):
            if not( IsFloat( f_cnt[0, i])):
                Mr_Exit('Please check the center of FILAMENT, 9 row, ', str(i + 4),' column in CAE_1E.xlsx')

            if not( IsFloat( f_cnt[1, i])):
                Mr_Exit('Please check the center of FILAMENT, 10 row, ', str(i + 4),' column in CAE_1E.xlsx')

            if not( IsFloat( f_cnt[2, i])):
                Mr_Exit('Please check the center of FILAMENT, 11 row, ', str(i + 4),' column in CAE_1E.xlsx')
    
    #
    f_operator = Cells(fil_Nmbr)
    f_center   = Cells(fil_Nmbr)
    
    for i in range(0, fil_Nmbr):
        if not( IsFloat( f_l[i])) or ( f_l[i] < 0):
            Mr_Exit('Please check the length of FILAMENT, 13 row, ', str(i + 4),' column in CAE_1E.xlsx')

        if not( IsFloat( f_pphi[i])):
            Mr_Exit('Please check the Euler angles of FILAMENT, 19 row, ', str(i + 4),' column in CAE_1E.xlsx')

        if not( IsFloat( f_ppye[i])):
            Mr_Exit('Please check the Euler angles of FILAMENT, 17 row, ', str(i + 4),' column in CAE_1E.xlsx')

        if not( IsFloat( f_ttta[i])):
            Mr_Exit('Please check the Euler angles of FILAMENT, 18 row, ', str(i + 4),' column in CAE_1E.xlsx')

        if not( IsFloat( f_lmbd[i])):
            Mr_Exit('Please check the charge density of FILAMENT, 21 row, ', str(i + 4),' column in CAE_1E.xlsx')

        # degree to radian converter
        f_ttta[i] = deg2rad( f_ttta[i] )
        f_pphi[i] = deg2rad( f_pphi[i] )
        f_ppye[i] = deg2rad( f_ppye[i] )
        
        # separation
        f_operator[i] = Euler(f_ttta[i], f_pphi[i], f_ppye[i])
        f_center[i]   = f_cnt[:,i]
        
    print('   Done.')

#%% slab
torus_Nmbr = int(data.cell_value(25, 3))

if torus_Nmbr > 0:
    print('   Importing Slabs...')
    
    torus_cnt = array([MatAll(32, torus_Nmbr), MatAll(33, torus_Nmbr), MatAll(34, torus_Nmbr)])
    torus_Ro = MatAll(36, torus_Nmbr)
    torus_r  = MatAll(37, torus_Nmbr)
    # euler transformation angels
    torus_ppye = MatAll(42, torus_Nmbr)
    torus_ttta = MatAll(43, torus_Nmbr)
    torus_pphi = MatAll(44, torus_Nmbr)
    # charge density in slab
    torus_sgma = MatAll(46, torus_Nmbr)

    if not(isinstance(torus_cnt[0, 0], float)):
        for i in range(0, torus_Nmbr):
            if not( IsFloat( torus_cnt[0, i])):
                Mr_Exit('Please check the center of RECTANGULAR SLAB, 32 row, ', str(i + 4), ' column in CAE_1E.xlsx')

            if not( IsFloat( torus_cnt[1, i])):
                Mr_Exit('Please check the center of RECTANGULAR SLAB, 33 row, ', str(i + 4), ' column in CAE_1E.xlsx')

            if not( IsFloat( torus_cnt[2, i])):
                Mr_Exit('Please check the center of RECTANGULAR SLAB, 34 row, ', str(i + 4), ' column in CAE_1E.xlsx')
    
    #
    torus_operator = Cells(torus_Nmbr)
    torus_center   = Cells(torus_Nmbr)
    
    for i in range(0, torus_Nmbr):
        if not( IsFloat( torus_Ro[i])) or ( torus_Ro[i] < 0):
            Mr_Exit('Please check the length of RECTANGULAR SLAB, 36 row, ', str(i + 4), ' column in CAE_1E.xlsx')

        if not( IsFloat( torus_r[i])) or ( torus_r[i] < 0):
            Mr_Exit('Please check the breadth of RECTANGULAR SLAB, 37 row, ', str(i + 4), ' column in CAE_1E.xlsx')

        if not( IsFloat( torus_pphi[i])):
            Mr_Exit('Please check the Euler angles of RECTANGULAR SLAB, 44 row, ', str(i + 4), ' column in CAE_1E.xlsx')

        if not( IsFloat( torus_ppye[i])):
            Mr_Exit('Please check the Euler angles of RECTANGULAR SLAB, 42 row, ', str(i + 4), ' column in CAE_1E.xlsx')

        if not( IsFloat( torus_ttta[i])):
            Mr_Exit('Please check the Euler angles of RECTANGULAR SLAB, 43 row, ', str(i + 4), ' column in CAE_1E.xlsx')

        if not( IsFloat( torus_sgma[i])):
            Mr_Exit('Please check the current density of RECTANGULAR SLAB, 46 row, ', str(i + 4), ' column in CAE_1E.xlsx')

        # degree to radian converter
        torus_ttta[i] = deg2rad( torus_ttta[i])
        torus_pphi[i] = deg2rad( torus_pphi[i])
        torus_ppye[i] = deg2rad( torus_ppye[i])
        
        # separation 
        torus_operator[i] = Euler(torus_ttta[i], torus_pphi[i], torus_ppye[i])
        torus_center[i] = torus_cnt[:,i]
        
    print('   Done.')


if ((fil_Nmbr == 0) and (torus_Nmbr == 0) ):
    Mr_Exit('Please enter the data for atleast an electromagnet, in CAE_1E.xlsx')
    
#%% Points at which electric field to be calculated :-

if CAET_run == 0:

    if int( data_2.cell_value(4,3)) == 1:
        cart_pnt = int(1)
        polar_pnt = int(0)
    elif int( data_2.cell_value(4,3)) == 0:
        cart_pnt = int(0)
        polar_pnt = int(1)
        
    print('   Imorting space point data...')
    
    if (cart_pnt == 1):
        if ( not( IsFloat( data_2.cell_value( 7,3))) or not( IsFloat( data_2.cell_value( 8,3))) or not( IsFloat( data_2.cell_value( 7,4))) or not( IsFloat( data_2.cell_value( 8,4))) or not( IsFloat( data_2.cell_value( 7,5))) or not( IsFloat( data_2.cell_value( 8,5))) or not( IsFloat( data_2.cell_value( 10,3))) or not( IsFloat( data_2.cell_value( 10,4))) or not( IsFloat( data_2.cell_value( 10,5))) ):
            Mr_Exit('Please check the space data in CAE_2E.xlsx')
    
        if ( ( data_2.cell_value(8,3) > data_2.cell_value(7,3)) or ( data_2.cell_value(8,4) > data_2.cell_value(7,4)) or ( data_2.cell_value(8,5) > data_2.cell_value(7,5)) ):
            Mr_Exit('Please enter the space data with min and max correctly in CAE_2E.xlsx')

        if ( ( data_2.cell_value(10,3) <= 0) or ( data_2.cell_value(10,4) <= 0) or ( data_2.cell_value(10,5) <= 0)):
            Mr_Exit('Please enter the space data correctly in CAE_2E.xlsx')
    
        xx = linspace( data_2.cell_value(8,3), data_2.cell_value(7,3), int(data_2.cell_value(10,3)) )
        yy = linspace( data_2.cell_value(8,4), data_2.cell_value(7,4), int(data_2.cell_value(10,4)) )
        zz = linspace( data_2.cell_value(8,5), data_2.cell_value(7,5), int(data_2.cell_value(10,5)) )
        #
        xx, yy, zz = meshgrid(xx, yy, zz)
    
    if (polar_pnt ==1):
        if ( not( IsFloat( data_2.cell_value( 7,3))) or not( IsFloat( data_2.cell_value( 8,3))) or not( IsFloat( data_2.cell_value( 7,4))) or not( IsFloat( data_2.cell_value( 8,4))) or not( IsFloat( data_2.cell_value( 7,5))) or not( IsFloat( data_2.cell_value( 8,5))) or not( IsFloat( data_2.cell_value( 10,3))) or not( IsFloat( data_2.cell_value( 10,4))) or not( IsFloat( data_2.cell_value( 10,5))) ):
            Mr_Exit('Please check the space data in CAE_2E.xlsx')
    
        if ( ( data_2.cell_value(8,3) > data_2.cell_value(7,3)) or ( data_2.cell_value(8,4) > data_2.cell_value(7,4)) or ( data_2.cell_value(8,5) > data_2.cell_value(7,5)) ):
            Mr_Exit('Please enter the space data with min and max correctly in CAE_2E.xlsx')
    
        if ( ( data_2.cell_value(10,3) <= 0) or ( data_2.cell_value(10,4) <= 0) or ( data_2.cell_value(10,5) <= 0)):
            Mr_Exit('Please enter the space data correctly in CAE_2E.xlsx')
    
        # degree to radian converter
        degdata1 = deg2rad(data_2.cell_value(8,4))
        degdata2 = deg2rad(data_2.cell_value(7,4))
        degpoint = int(data_2.cell_value(10,4))
        #
        rr     = linspace( data_2.cell_value(8,3), data_2.cell_value(7,3), int(data_2.cell_value(10,3)) )
        phyphy = linspace( degdata1, degdata2, degpoint)
        zz     = linspace( data_2.cell_value(8,5), data_2.cell_value(7,5), int(data_2.cell_value(10,5)) )
        #
        rr, phyphy, zz = meshgrid(rr,phyphy,zz)
        xx,   yy,   zz = Pol2Cart(rr, phyphy, zz)
    
    print('   Done.')
    
    #%% Contour Plot option
    need_contour_plot = data_2.cell_value(30,6)
    
    if (type(need_contour_plot) == str):
        Mr_Exit('Please check the contour plot option, 31 row, 7 column in CAE_2E.xlsx')
    elif (need_contour_plot != 1.0) and (need_contour_plot != 0.0):
        Mr_Exit('Please check the contour plot option, 31 row, 7 column in CAE_2E.xlsx')
    
    if need_contour_plot == 1:
        print('   Checking conditions for contour plots...')
        
        if cart_pnt == 1:
            logic_x = ( data_2.cell_value(8,3) == data_2.cell_value(7,3) )
            logic_y = ( data_2.cell_value(8,4) == data_2.cell_value(7,4) )
            logic_z = ( data_2.cell_value(8,5) == data_2.cell_value(7,5) )
            if ( (logic_x and logic_y) or (logic_y and logic_z) or (logic_x and logic_z) or not(logic_x or logic_y or logic_z) ):
                Mr_Exit('To draw a contour plot one parameter must be constant, 8 and 9 row in CAE_2E.xlsx')
            
        elif polar_pnt == 1:
            logic_r = ( data_2.cell_value(8,3) == data_2.cell_value(7,3) )
            logic_t = ( data_2.cell_value(8,4) == data_2.cell_value(7,4) )
            logic_z = ( data_2.cell_value(8,5) == data_2.cell_value(7,5) )
            if not( not(logic_r) and not(logic_r) and logic_z):
                Mr_Exit('To draw a contour plot in cylindrical system z must be constant, 8 and 9 row in CAE_2E.xlsx')
    
    if need_contour_plot == 1:
        if cart_pnt == 1:
            
            if logic_x == True:        # Here the x is constant
                cont_xx = linspace( data_2.cell_value(8,4), data_2.cell_value(7,4), int(data_2.cell_value(10,4)) ) # y
                cont_yy = linspace( data_2.cell_value(8,5), data_2.cell_value(7,5), int(data_2.cell_value(10,5)) ) # z
                cont_x, cont_y = meshgrid(cont_xx, cont_yy)
                #
                cont_pp, cont_qq = cont_x.shape
            
            elif logic_y == True:      # Here the y is constant
                cont_xx = linspace( data_2.cell_value(8,5), data_2.cell_value(7,5), int(data_2.cell_value(10,5)) ) # z
                cont_yy = linspace( data_2.cell_value(8,3), data_2.cell_value(7,3), int(data_2.cell_value(10,3)) ) # x
                cont_x, cont_y = meshgrid(cont_xx, cont_yy)
                #
                cont_pp, cont_qq = cont_x.shape
            
            elif logic_z == True:      # Here the z is constant
                cont_xx = linspace( data_2.cell_value(8,3), data_2.cell_value(7,3), int(data_2.cell_value(10,3)) ) # x
                cont_yy = linspace( data_2.cell_value(8,4), data_2.cell_value(7,4), int(data_2.cell_value(10,4)) ) # y
                cont_x, cont_y = meshgrid(cont_xx, cont_yy)
                #
                cont_pp, cont_qq = cont_x.shape
                
        elif polar_pnt == 1:
            if logic_z == True:        # Here the z is constant
                cont_RR = linspace( data_2.cell_value(8,3), data_2.cell_value(7,3), int(data_2.cell_value(10,3)) )      # Same as 'rr'
                cont_TT = linspace( degdata1, degdata2, degpoint)                                                       # Same as 'phyphy'
                cont_RR, cont_TT = meshgrid(cont_RR, cont_TT)
                #
                cont_x = cont_RR * cos(cont_TT)
                cont_y = cont_RR * sin(cont_TT)
                #
                cont_pp, cont_qq = cont_x.shape
        
        print('   Done.')
    
    q, w, e = xx.shape
    s = q* w* e
    point_1 = Cells(s)
    point_2 = Cells(s)
    point_3 = Cells(s)
    
    counter = 0
    for i in range(0,q):
        for j in range(0,w):
            for k in range(0,e):
                point_1[counter] = xx[i, j, k]
                point_2[counter] = yy[i, j, k]
                point_3[counter] = zz[i, j, k]
                counter = counter + 1
    
    point = array([point_1,point_2,point_3])
else:
    # magnet codes will directly receive the data from bridge.txt
    pass

print('   All the input parameters are validated...')
print('   -----------------------------------------------------------------------------------   ')
