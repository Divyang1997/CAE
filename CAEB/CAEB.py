#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pandas as pd
from numpy import pi, sin, cos, arange, array, nan

sys.path.append('/home/divyang/Documents/Python/CAE/CAE v1.2/Magnet_modules')

def Cells(q):

    A = []
    for i in range(q):
        A.append([])
    
    return A

def Mr_Exit(*args):
    print('\033[1;31;40m ERROR!!! \033[0m : ' + ''.join(args) )
    sys.exit(1)

Num = int( input('Please enter the number of coils: '))

#%%
def Cent_generator(cent, Num):
    # cent is the center's co-ordinate of magnet
    # Num is number of coil
    
    if cent[1] != 0:
        Mr_Exit('CAE BUILDER needs y zero, Please check CAE_1.xlsx')
    
    # as y == 0, x - co-ordinate is the radius in cylindrical system
    rn = cent[0]
    
    # define angles at which we need the coil
    ang = arange(0, 2* pi, 2* pi/ Num)
    
    # New co-ordinates
    xc = ( rn* cos(ang) ).tolist()
    yc = ( rn* sin(ang) ).tolist()
    zc = [cent[2]]* len(xc)
    
    return xc, yc, zc


from B_Inputs import f_Nmbr, s_Nmbr, a_Nmbr, ra_Nmbr, cy_Nmbr

# filament
if f_Nmbr > 0:
    from B_Inputs import f_cnt, f_l, f_ppye, f_ttta, f_pphi, f_I
    
    f_cnt_x = Cells(len(f_cnt[0,:]))
    f_cnt_y = Cells(len(f_cnt[0,:]))
    f_cnt_z = Cells(len(f_cnt[0,:]))
    
    for i in range(0, len(f_cnt[0,:])):
        f_cnt_x[i], f_cnt_y[i], f_cnt_z[i] = Cent_generator( [f_cnt[0,i], f_cnt[1,i], f_cnt[2,i]], Num)

# slabs  
if s_Nmbr > 0:
    from B_Inputs import s_cnt, s_l, s_b, s_h, s_ppye, s_ttta, s_pphi, s_J
    
    s_cnt_x = Cells(len(s_cnt[0,:]))
    s_cnt_y = Cells(len(s_cnt[0,:]))
    s_cnt_z = Cells(len(s_cnt[0,:]))
    
    for i in range(0, len(s_cnt[0,:])):
        s_cnt_x[i], s_cnt_y[i], s_cnt_z[i] = Cent_generator( [s_cnt[0,i], s_cnt[1,i], s_cnt[2,i]], Num)
        
# arc 
if a_Nmbr > 0:
    from B_Inputs import a_cnt, a_a0, a_phy_1, a_phy_2, a_ppye, a_ttta, a_pphi, a_I
    
    a_cnt_x = Cells(len(a_cnt[0,:]))
    a_cnt_y = Cells(len(a_cnt[0,:]))
    a_cnt_z = Cells(len(a_cnt[0,:]))
    
    for i in range(0, len(a_cnt[0,:])):
        a_cnt_x[i], a_cnt_y[i], a_cnt_z[i] = Cent_generator( [a_cnt[0,i], a_cnt[1,i], a_cnt[2,i]], Num)
        
# ractangular arc
if ra_Nmbr > 0:
    from B_Inputs import ra_cnt, ra_r_inner, ra_r_outer, ra_phy_1, ra_phy_2, ra_tckns, ra_ppye, ra_ttta, ra_pphi, ra_J
    
    ra_cnt_x = Cells(len(ra_cnt[0,:]))
    ra_cnt_y = Cells(len(ra_cnt[0,:]))
    ra_cnt_z = Cells(len(ra_cnt[0,:]))
    
    for i in range(0, len(ra_cnt[0,:])):
        ra_cnt_x[i], ra_cnt_y[i], ra_cnt_z[i] = Cent_generator( [ra_cnt[0,i], ra_cnt[1,i], ra_cnt[2,i]], Num)
        
# ractangular arc
if cy_Nmbr > 0:
    from B_Inputs import cy_cnt, cy_ro, cy_l, cy_ppye, cy_ttta, cy_pphi, cy_J
    
    cy_cnt_x = Cells(len(cy_cnt[0,:]))
    cy_cnt_y = Cells(len(cy_cnt[0,:]))
    cy_cnt_z = Cells(len(cy_cnt[0,:]))
    
    for i in range(0, len(cy_cnt[0,:])):
        cy_cnt_x[i], cy_cnt_y[i], cy_cnt_z[i] = Cent_generator( [cy_cnt[0,i], cy_cnt[1,i], cy_cnt[2,i]], Num)

#%%
def duplicate(chunk, times):
    return list(chunk)*times

def z_rotations(pphi, Nm):
    if pphi != 0:
        Mr_Exit('CAE BUILDER needs zero z ratation, Please check CAE_1.xlsx')

    ang = arange(0, 2* pi, 2* pi/ Nm)
    return ang

def Rearrange( List, Nm):
    
    ''' Nm is number of coils '''
    
    New_list = Cells(Nm)
    
    for i in range(0, Nm):
        New_list[i] = List[ i : len(List) : Nm ]
        
    New_list = array(New_list).reshape([1,len(List)])
        
    return New_list

#%%
# filament
if f_Nmbr > 0:
    f_dub_l    = Cells(len(f_l))
    f_dub_ppye = Cells(len(f_l))
    f_dub_ttta = Cells(len(f_l))
    f_dub_pphi = Cells(len(f_l))
    f_dub_I    = Cells(len(f_l))
    
    for i in range(0, len(f_l)):
        f_dub_l[i]    = duplicate([f_l[i]], Num)
        f_dub_ppye[i] = duplicate([f_ppye[i]], Num) 
        f_dub_ttta[i] = duplicate([f_ttta[i]], Num)
        f_dub_pphi[i] = z_rotations(f_pphi[i], Num)
        f_dub_I[i]    = duplicate([f_I[i]], Num)

# slabs  
if s_Nmbr > 0:
    s_dub_l    = Cells(len(s_l))
    s_dub_b    = Cells(len(s_l))
    s_dub_h    = Cells(len(s_l))
    s_dub_ppye = Cells(len(s_l))
    s_dub_ttta = Cells(len(s_l))
    s_dub_pphi = Cells(len(s_l))
    s_dub_J    = Cells(len(s_l))
    
    for i in range(0, len(s_l)):
        s_dub_l[i]    = duplicate([s_l[i]], Num)
        s_dub_b[i]    = duplicate([s_b[i]], Num)
        s_dub_h[i]    = duplicate([s_h[i]], Num)
        s_dub_ppye[i] = duplicate([s_ppye[i]], Num)
        s_dub_ttta[i] = duplicate([s_ttta[i]], Num)
        s_dub_pphi[i] = z_rotations(s_pphi[i], Num)
        s_dub_J[i]    = duplicate([s_J[i]], Num)

# arc 
if a_Nmbr > 0:
    a_dub_a0    = Cells(len(a_a0))
    a_dub_phy_1 = Cells(len(a_a0))
    a_dub_phy_2 = Cells(len(a_a0))
    a_dub_ppye  = Cells(len(a_a0))
    a_dub_ttta  = Cells(len(a_a0))
    a_dub_pphi  = Cells(len(a_a0))
    a_dub_I     = Cells(len(a_a0))
    
    for i in range(0, len(a_a0)):
        a_dub_a0[i]    = duplicate([a_a0[i]], Num)
        a_dub_phy_1[i] = duplicate([a_phy_1[i]], Num)
        a_dub_phy_2[i] = duplicate([a_phy_2[i]], Num)
        a_dub_ppye[i]  = duplicate([a_ppye[i]], Num)
        a_dub_ttta[i]  = duplicate([a_ttta[i]], Num)
        a_dub_pphi[i]  = z_rotations(a_pphi[i], Num)
        a_dub_I[i]     = duplicate([a_I[i]], Num)
        
# ractangular arc
if ra_Nmbr > 0:
     ra_dub_r_inner = Cells(len(ra_J))
     ra_dub_r_outer = Cells(len(ra_J))
     ra_dub_phy_1   = Cells(len(ra_J))
     ra_dub_phy_2   = Cells(len(ra_J))
     ra_dub_tckns   = Cells(len(ra_J))
     ra_dub_ppye    = Cells(len(ra_J))
     ra_dub_ttta    = Cells(len(ra_J))
     ra_dub_pphi    = Cells(len(ra_J))
     ra_dub_J       = Cells(len(ra_J))
     
     for i in range(0, len(ra_J)):
         ra_dub_r_inner[i] = duplicate([ra_r_inner[i]], Num)
         ra_dub_r_outer[i] = duplicate([ra_r_outer[i]], Num)
         ra_dub_phy_1[i]   = duplicate([ra_phy_1[i]], Num)
         ra_dub_phy_2[i]   = duplicate([ra_phy_2[i]], Num)
         ra_dub_tckns[i]   = duplicate([ra_tckns[i]], Num)
         ra_dub_ppye[i]    = duplicate([ra_ppye[i]], Num)
         ra_dub_ttta[i]    = duplicate([ra_ttta[i]], Num)
         ra_dub_pphi[i]    = z_rotations(ra_pphi[i], Num)
         ra_dub_J[i]       = duplicate([ra_J[i]], Num)

# Cylindrical wire
if cy_Nmbr > 0:
    cy_dub_ro   = Cells(len(cy_ro))
    cy_dub_l    = Cells(len(cy_ro))
    cy_dub_ppye = Cells(len(cy_ro))
    cy_dub_ttta = Cells(len(cy_ro))
    cy_dub_pphi = Cells(len(cy_ro))
    cy_dub_J    = Cells(len(cy_ro))
    
    for i in range(0, len(cy_ro)):
        cy_dub_ro[i]   = duplicate([cy_ro[i]], Num)
        cy_dub_l[i]    = duplicate([cy_l[i]], Num)
        cy_dub_ppye[i] = duplicate([cy_ppye[i]], Num)
        cy_dub_ttta[i] = duplicate([cy_ttta[i]], Num)
        ra_dub_pphi[i] = z_rotations(cy_pphi[i], Num)
        cy_dub_J[i]    = duplicate([cy_J[i]], Num)     

#%%
def one_line(data):
    data = array(data)
    if len(data) == 1:
        s = 1
        var = data
    else:
        s = data.shape[0]* data.shape[1]
        var = data.reshape([1,s])
        
    return var

# filament
if f_Nmbr > 0:
    f_cnt_x    = one_line(f_cnt_x)
    f_cnt_y    = one_line(f_cnt_y)
    f_cnt_z    = one_line(f_cnt_z)
    f_dub_l    = one_line(f_dub_l)
    f_dub_ppye = one_line(f_dub_ppye)* (180/ pi)
    f_dub_ttta = one_line(f_dub_ttta)* (180/ pi)
    f_dub_pphi = one_line(f_dub_pphi)* (180/ pi)
    f_dub_I    = one_line(f_dub_I)
    
    f_cnt_x    = Rearrange(f_cnt_x[0], Num)
    f_cnt_y    = Rearrange(f_cnt_y[0], Num)
    f_cnt_z    = Rearrange(f_cnt_z[0], Num)
    f_dub_l    = Rearrange(f_dub_l[0], Num)
    f_dub_ppye = Rearrange(f_dub_ppye[0], Num)
    f_dub_ttta = Rearrange(f_dub_ttta[0], Num)
    f_dub_pphi = Rearrange(f_dub_pphi[0], Num)
    f_dub_I    = Rearrange(f_dub_I[0], Num)
    
    output_f_data = pd.DataFrame({'x': ( ['X center'] + list(f_cnt_x[0]) ),
                                  'y': ( ['Y center'] + list(f_cnt_y[0]) ),
                                  'z': ( ['Z center'] + list(f_cnt_z[0]) ),
                                  'empty1': nan,
                                  'l': ( ['Length'] + list(f_dub_l[0]) ),
                                  'empty2': nan,
                                  'ppye': ( ['X rotation'] + list(f_dub_ppye[0]) ),
                                  'ttta': ( ['Y rotation'] + list(f_dub_ttta[0]) ),
                                  'pphi': ( ['Z rotation'] + list(f_dub_pphi[0]) ),
                                  'empty3': nan,
                                  'I':   ( ['Current'] + list(f_dub_I[0]) )}, columns=['x','y','z','empty1','l','empty2','ppye','ttta','pphi','empty3','I'] )

    output_f_data = output_f_data.rename(columns={'empty1':'', 'empty2':'','empty3':''})
    output_f_data = output_f_data.transpose()

# slabs  
if s_Nmbr > 0:
    s_cnt_x    = one_line(s_cnt_x)
    s_cnt_y    = one_line(s_cnt_y)
    s_cnt_z    = one_line(s_cnt_z)
    s_dub_l    = one_line(s_dub_l)
    s_dub_b    = one_line(s_dub_b)
    s_dub_h    = one_line(s_dub_h)
    s_dub_ppye = one_line(s_dub_ppye)* (180/ pi)
    s_dub_ttta = one_line(s_dub_ttta)* (180/ pi)
    s_dub_pphi = one_line(s_dub_pphi)* (180/ pi)
    s_dub_J    = one_line(s_dub_J)
    
    s_cnt_x    = Rearrange(s_cnt_x[0], Num)
    s_cnt_y    = Rearrange(s_cnt_y[0], Num)
    s_cnt_z    = Rearrange(s_cnt_z[0], Num)
    s_dub_l    = Rearrange(s_dub_l[0], Num)
    s_dub_b    = Rearrange(s_dub_b[0], Num)
    s_dub_h    = Rearrange(s_dub_h[0], Num)
    s_dub_ppye = Rearrange(s_dub_ppye[0], Num)
    s_dub_ttta = Rearrange(s_dub_ttta[0], Num)
    s_dub_pphi = Rearrange(s_dub_pphi[0], Num)
    s_dub_J    = Rearrange(s_dub_J[0], Num)

    output_s_data = pd.DataFrame({'x': ( ['X center'] + list(s_cnt_x[0]) ),
                                  'y': ( ['Y center'] + list(s_cnt_y[0]) ),
                                  'z': ( ['Z center'] + list(s_cnt_z[0]) ),
                                  'empty1': nan,
                                  'l': ( ['Length'] + list(s_dub_l[0]) ),
                                  'b': ( ['Breadth'] + list(s_dub_b[0]) ),
                                  'h': ( ['Height'] + list(s_dub_h[0]) ),
                                  'empty2': nan,
                                  'ppye': ( ['X rotation'] + list(s_dub_ppye[0]) ),
                                  'ttta': ( ['Y rotation'] + list(s_dub_ttta[0]) ),
                                  'pphi': ( ['Z rotation'] + list(s_dub_pphi[0]) ),
                                  'empty3': nan,
                                  'J':   ( ['Current density'] + list(s_dub_J[0]) )}, columns=['x','y','z','empty1','l','b','h','empty2','ppye','ttta','pphi','empty3','J'] )

    output_s_data = output_s_data.rename(columns={'empty1':'', 'empty2':'','empty3':''})
    output_s_data = output_s_data.transpose()

# arcs
if a_Nmbr > 0:
    a_cnt_x     = one_line(a_cnt_x)
    a_cnt_y     = one_line(a_cnt_y)
    a_cnt_z     = one_line(a_cnt_z)
    a_dub_a0    = one_line(a_dub_a0)
    a_dub_phy_1 = one_line(a_dub_phy_1)* (180/ pi)
    a_dub_phy_2 = one_line(a_dub_phy_2)* (180/ pi)
    a_dub_ppye  = one_line(a_dub_ppye)* (180/ pi)
    a_dub_ttta  = one_line(a_dub_ttta)* (180/ pi)
    a_dub_pphi  = one_line(a_dub_pphi)* (180/ pi)
    a_dub_I     = one_line(a_dub_I)
    
    a_cnt_x     = Rearrange(a_cnt_x[0], Num)
    a_cnt_y     = Rearrange(a_cnt_y[0], Num)
    a_cnt_z     = Rearrange(a_cnt_z[0], Num)
    a_dub_a0    = Rearrange(a_dub_a0[0], Num)
    a_dub_phy_1 = Rearrange(a_dub_phy_1[0], Num)
    a_dub_phy_2 = Rearrange(a_dub_phy_2[0], Num)
    a_dub_ppye  = Rearrange(a_dub_ppye[0], Num)
    a_dub_ttta  = Rearrange(a_dub_ttta[0], Num)
    a_dub_pphi  = Rearrange(a_dub_pphi[0], Num)
    a_dub_I     = Rearrange(a_dub_I[0], Num)
    
    output_a_data = pd.DataFrame({'x': ( ['X center'] + list(a_cnt_x[0]) ),
                                  'y': ( ['Y center'] + list(a_cnt_y[0]) ),
                                  'z': ( ['Z center'] + list(a_cnt_z[0]) ),
                                  'empty1': nan,
                                  'r': ( ['Radius'] + list(a_dub_a0[0]) ),
                                  'phy1': ( ['Angle 1'] + list(a_dub_phy_1[0]) ),
                                  'phy2': ( ['Angle 2'] + list(a_dub_phy_2[0]) ),
                                  'empty2': nan,
                                  'ppye': ( ['X rotation'] + list(a_dub_ppye[0]) ),
                                  'ttta': ( ['Y rotation'] + list(a_dub_ttta[0]) ),
                                  'pphi': ( ['Z rotation'] + list(a_dub_pphi[0]) ),
                                  'empty3': nan,
                                  'I':   ( ['Current'] + list(a_dub_I[0]) )}, columns=['x','y','z','empty1','r','phy1','phy2','empty2','ppye','ttta','pphi','empty3','I'] )

    output_a_data = output_a_data.rename(columns={'empty1':'', 'empty2':'','empty3':''})
    output_a_data = output_a_data.transpose()
    
# ractangular arc
if ra_Nmbr > 0:
    ra_cnt_x       = one_line(ra_cnt_x)
    ra_cnt_y       = one_line(ra_cnt_y)
    ra_cnt_z       = one_line(ra_cnt_z)
    ra_dub_r_inner = one_line(ra_dub_r_inner)
    ra_dub_r_outer = one_line(ra_dub_r_outer)
    ra_dub_phy_1   = one_line(ra_dub_phy_1)* (180/ pi)
    ra_dub_phy_2   = one_line(ra_dub_phy_2)* (180/ pi)
    ra_dub_tckns   = one_line(ra_dub_tckns)
    ra_dub_ppye    = one_line(ra_dub_ppye)* (180/ pi)
    ra_dub_ttta    = one_line(ra_dub_ttta)* (180/ pi)
    ra_dub_pphi    = one_line(ra_dub_pphi)* (180/ pi)
    ra_dub_J       = one_line(ra_dub_J)
    
    ra_cnt_x       = Rearrange(ra_cnt_x[0], Num)
    ra_cnt_y       = Rearrange(ra_cnt_y[0], Num)
    ra_cnt_z       = Rearrange(ra_cnt_z[0], Num)
    ra_dub_r_inner = Rearrange(ra_dub_r_inner[0], Num)
    ra_dub_r_outer = Rearrange(ra_dub_r_outer[0], Num)
    ra_dub_phy_1   = Rearrange(ra_dub_phy_1[0], Num)
    ra_dub_phy_2   = Rearrange(ra_dub_phy_2[0], Num)
    ra_dub_tckns   = Rearrange(ra_dub_tckns[0], Num)
    ra_dub_ppye    = Rearrange(ra_dub_ppye[0], Num)
    ra_dub_ttta    = Rearrange(ra_dub_ttta[0], Num)
    ra_dub_pphi    = Rearrange(ra_dub_pphi[0], Num)
    ra_dub_J       = Rearrange(ra_dub_J[0], Num)
    
    output_ra_data = pd.DataFrame({'x': ( ['X center'] + list(ra_cnt_x[0]) ),
                                   'y': ( ['Y center'] + list(ra_cnt_y[0]) ),
                                   'z': ( ['Z center'] + list(ra_cnt_z[0]) ),
                                   'empty1': nan,
                                   'r_inner': ( ['Inner radius'] + list(ra_dub_r_inner[0]) ),
                                   'r_outer': ( ['Outer radius'] + list(ra_dub_r_outer[0]) ),
                                   'phy1': ( ['Angle 1'] + list(ra_dub_phy_1[0]) ),
                                   'phy2': ( ['Angle 2'] + list(ra_dub_phy_2[0]) ),
                                   'thickness': ( ['Thickness'] + list(ra_dub_tckns[0]) ),
                                   'empty2': nan,
                                   'ppye': ( ['X rotation'] + list(ra_dub_ppye[0]) ),
                                   'ttta': ( ['Y rotation'] + list(ra_dub_ttta[0]) ),
                                   'pphi': ( ['Z rotation'] + list(ra_dub_pphi[0]) ),
                                   'empty3': nan,
                                   'J':   ( ['Current density'] + list(ra_dub_J[0]) )}, columns=['x','y','z','empty1','r_inner','r_outer','phy1','phy2','thickness','empty2','ppye','ttta','pphi','empty3','J'] )

    output_ra_data = output_ra_data.rename(columns={'empty1':'', 'empty2':'','empty3':''})
    output_ra_data = output_ra_data.transpose()
     
# Cylindrical wire
if cy_Nmbr > 0:
    cy_cnt_x    = one_line(cy_cnt_x)
    cy_cnt_y    = one_line(cy_cnt_y)
    cy_cnt_z    = one_line(cy_cnt_z)
    cy_dub_ro   = one_line(cy_dub_ro)
    cy_dub_l    = one_line(cy_dub_l)
    cy_dub_ppye = one_line(cy_dub_ppye)* (180/ pi)
    cy_dub_ttta = one_line(cy_dub_ttta)* (180/ pi)
    cy_dub_pphi = one_line(cy_dub_pphi)* (180/ pi)
    cy_dub_J    = one_line(cy_dub_J)
    
    cy_cnt_x    = Rearrange(cy_cnt_x[0], Num)
    cy_cnt_y    = Rearrange(cy_cnt_y[0], Num)
    cy_cnt_z    = Rearrange(cy_cnt_z[0], Num)
    cy_dub_ro   = Rearrange(cy_dub_ro[0], Num)
    cy_dub_l    = Rearrange(cy_dub_l[0], Num)
    cy_dub_ppye = Rearrange(cy_dub_ppye[0], Num)
    cy_dub_ttta = Rearrange(cy_dub_ttta[0], Num)
    cy_dub_pphi = Rearrange(cy_dub_pphi[0], Num)
    cy_dub_J    = Rearrange(cy_dub_J[0], Num)

    output_cy_data = pd.DataFrame({'x': ( ['X center'] + list(cy_cnt_x[0]) ),
                                   'y': ( ['Y center'] + list(cy_cnt_y[0]) ),
                                   'z': ( ['Z center'] + list(cy_cnt_z[0]) ),
                                   'empty1': nan,
                                   'r': ( ['Radius'] + list(cy_dub_ro[0]) ),
                                   'l': ( ['Length'] + list(cy_dub_l[0]) ),
                                   'empty2': nan,
                                   'ppye': ( ['X rotation'] + list(cy_dub_ppye[0]) ),
                                   'ttta': ( ['Y rotation'] + list(cy_dub_ttta[0]) ),
                                   'pphi': ( ['Z rotation'] + list(cy_dub_pphi[0]) ),
                                   'empty3': nan,
                                   'J':   ( ['Current density'] + list(cy_dub_J[0]) )}, columns=['x','y','z','empty1','r','l','empty2','ppye','ttta','pphi','empty3','J'] )

    output_cy_data = output_cy_data.rename(columns={'empty1':'', 'empty2':'','empty3':''})
    output_cy_data = output_cy_data.transpose()

#%%
# Create a Pandas Excel writer using XlsxWriter as the engine.
writer = pd.ExcelWriter('Builder.xlsx', engine='xlsxwriter')

if f_Nmbr > 0:
    output_f_data.to_excel(writer, index = False, sheet_name = 'Filaments')
if s_Nmbr > 0:
    output_s_data.to_excel(writer, index = False, sheet_name = 'Slabs')
if a_Nmbr > 0:
    output_a_data.to_excel(writer, index = False, sheet_name = 'Arcs')
if ra_Nmbr > 0:
    output_ra_data.to_excel(writer, index = False, sheet_name = 'Rectangular arcs')
if cy_Nmbr > 0:
    output_cy_data.to_excel(writer, index = False, sheet_name = 'Cylindrical wire')

# Close the Pandas Excel writer and output the Excel file.
writer.save()