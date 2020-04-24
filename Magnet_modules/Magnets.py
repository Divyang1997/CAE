#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from scipy.integrate import simps
from numpy.linalg import inv
from numpy import pi, sqrt, sin, cos, log, arctan, arcsinh, arctan2, trapz, \
                  array, arange, meshgrid, dot

def Cart2Pol(x, y, z):
    ρ = sqrt(x**2 + y**2)
    φ = arctan2(y, x)
    return ρ, φ, z

#%% Filament
def B_f(µ, Nmbr, l, I, operator, center, point):

    B = 0
    for i in range(0, Nmbr):
        pnt = dot( inv(operator[i]), point - center[i] )
        
        # (xo,yo,z1) and (xo,yo,z2) here xo = 0 and yo = 0
        z1 =-l[i]/2
        z2 = l[i]/2
        x  = 0
        y  = 0
        
        # Magnetic field calculation
        u   = pnt[0] - x
        v   = pnt[1] - y
        w1  = pnt[2] - z1
        w2  = pnt[2] - z2
        auv = sqrt(u**2 + v**2)

        Hx = ( v/ (auv**2))* ( -w1/ sqrt( auv**2 + w1**2) + w2/ sqrt( auv**2 + w2**2))
        Hy =-( u/ (auv**2))* ( -w1/ sqrt( auv**2 + w1**2) + w2/ sqrt( auv**2 + w2**2))
        Hz = 0

        B_i = (µ* I[i]/ (4*pi))* array([ Hx, Hy, Hz ])

        # Transformation 2 and Addition of magnetic field
        B = B + dot( operator[i], B_i)
    
    return B

#%% Slab
def tau(ui, vj, wk):
    return -(  wk* arcsinh( ui/ sqrt(wk**2 + vj**2) ) \
             + ui* arcsinh( wk/ sqrt(ui**2 + vj**2) ) \
             - vj* arctan( (ui/ vj)* (wk/ sqrt(ui**2 + vj**2 + wk**2) ) ) )

def B_s(µ, Nmbr, ver_org, J, operator, center, point):

    B = 0
    for i in range(0, Nmbr):
        pnt = dot( inv(operator[i]), point - center[i] )
    
        # Geometrically opposite points
        opp_pnt = array([ver_org[i][0], ver_org[i][6] ]).transpose()
        x = opp_pnt[0,:]
        y = opp_pnt[1,:]
        z = opp_pnt[2,:]
    
        # Magnetic field calculation
        u1 = pnt[0] - x[0]
        u2 = pnt[0] - x[1]
        v1 = pnt[1] - y[0]
        v2 = pnt[1] - y[1]
        w1 = pnt[2] - z[0]
        w2 = pnt[2] - z[1]
        #
        Hx = (- tau(u1, v1, w1) + tau(u1, v1, w2) + tau(u1, v2, w1) + tau(u2, v1, w1) 
              - tau(u1, v2, w2) - tau(u2, v1, w2) - tau(u2, v2, w1) + tau(u2, v2, w2) )
        
        Hy =-(- tau(v1, u1, w1) + tau(v1, u1, w2) + tau(v1, u2, w1) + tau(v2, u1, w1) 
              - tau(v1, u2, w2) - tau(v2, u1, w2) - tau(v2, u2, w1) + tau(v2, u2, w2) )
        Hz = 0
        
        B_i = (µ* J[i]/ (4* pi))* array([ Hx, Hy, Hz ])
    
        # Transformation 2 and Addition of magnetic field
        B = B + dot(operator[i] , B_i)
        
    return B

#%% Arc
def B_a(µ, Nmbr, a0, phy_1, phy_2, I, operator, cnt_cart, point):
    
    B = 0
    for i in range(0, Nmbr):
        pnt = dot( inv(operator[i]) , point - cnt_cart[i] )
        #
        r, phy, z = Cart2Pol( pnt[0], pnt[1], pnt[2] )
        
        # Constants
        b = sqrt( a0[i]**2 + r**2 + z**2 )
        phi_1 = phy_1[i] - phy
        phi_2 = phy_2[i] - phy
    
        # Integrals
        phi_n = arange( phi_1, phi_2 + 0.001, 0.001)
        #
        In_1 = simps(cos(phi_n)/ (b**2 - 2* r* a0[i]* cos(phi_n) )**(3/2) , phi_n)
        In_3 = simps(1/ (b**2 - 2* r* a0[i]* cos(phi_n) )**(3/2) , phi_n)
        if r == 0:
            In_2 = (1/ (a0[i]**2 + z**2)**(3/2) ) * ( cos(phi_1) - cos(phi_2) )
        else:
            In_2 =-(1/ (r* a0[i]) ) * ( 1/ sqrt(b**2 - 2* r* a0[i]* cos(phi_2) ) - 1/sqrt( b**2 - 2* r* a0[i]* cos(phi_1) ) )
    
        # Magnetic field calculation
        Hx = a0[i]* z* ( cos(phy)* In_1 - sin(phy)* In_2)
        Hy = a0[i]* z* ( sin(phy)* In_1 + cos(phy)* In_2)
        Hz = a0[i]* (a0[i]* In_3 - r* In_1)
        
        B_i = (µ* I[i]/ (4* pi))* array([ Hx, Hy, Hz ])
    
        # Transformation 2 and Addition of magnetic field
        B = B + dot(operator[i] , B_i)

    return B


#%% Ractangular Arc
def ra_I_Rr1( r,  ri, phi1, phi2, wk):
    phi   = arange( phi1, phi2 + 0.001, 0.001) 
    return simps( cos(phi)* sqrt( r**2 + ri**2 + wk**2 - 2* r* ri* cos(phi) ) +
                ( r* ( cos(phi)**2)* log( ri - r* cos(phi) + sqrt( ri**2 + r**2 + wk**2 - 2* ri* r* cos(phi) ) ) ), phi)

def ra_I_Rr2( r, ri, phi_i, wk):
    if r == 0:
        I_Rr2 = -cos(phi_i)* sqrt( ri**2 + wk**2)
    else:
        I_Rr2 = (1/ (2*r))* ( (ri - r* cos(phi_i) )* sqrt( r**2 + ri**2 + wk**2 - 2* r* ri* cos(phi_i) ) + 
                            ( wk**2 + (r* sin(phi_i) )**2 ) *
                            log( ri - r* cos(phi_i) + sqrt( r**2 + ri**2 + wk**2 - 2* r* ri* cos(phi_i) ) ) )
    return I_Rr2

def ra_I_Zr1( r, ri, phi, wk):
    return wk* simps( log( ri - r* cos(phi) + sqrt( r**2 + ri**2 + wk**2 - 2*r*ri* cos(phi) ) ), phi)

def ra_I_Zr2( r, r_ND, r_nd, phi_ND, phi_nd, wk):
    return wk* simps( simps( ( r**2 - r* r_ND* cos(phi_ND) ) / 
                           ( ( r**2 + r_ND**2 - 2* r* r_ND* cos(phi_ND) ) * 
                           sqrt( r**2 + r_ND**2 + wk**2 - 2* r* r_ND* cos(phi_ND) ) ), r_nd), phi_nd)

#
def B_ra(µ, Nmbr, r_inner, r_outer, phy_1, phy_2, tckns, J, operator, cnt_cart, point):
    
    B = 0
    for i in range(0, Nmbr):
        pnt = dot( inv(operator[i]) , point - cnt_cart[i] )
        #
        r, phy, z = Cart2Pol( pnt[0], pnt[1], pnt[2] )
        
        # Constants
        z_1 =-tckns[i]/2
        z_2 = tckns[i]/2
        r1  = r_inner[i]
        r2  = r_outer[i]
        
        phi_1 = phy_1[i] - phy
        phi_2 = phy_2[i] - phy
        w1    = z - z_1
        w2    = z - z_2
    
        # Integrals
        r_nd   = arange( r1, r2 + 0.01, 0.01)
        phi_n  = arange( phi_1, phi_2 + 0.01, 0.01)
        phi_nd = arange( phi_1, phi_2 + 0.01, 0.01)
        r_ND, phi_ND = meshgrid( r_nd, phi_nd)
        #
        I_Rr1 =   ra_I_Rr1(r, r1, phi_1, phi_2, w1) - ra_I_Rr1(r, r1, phi_1, phi_2, w2) \
                - ra_I_Rr1(r, r2, phi_1, phi_2, w1) + ra_I_Rr1(r, r2, phi_1, phi_2, w2)
        
        I_Rr2 = - ra_I_Rr2( r, r1, phi_1, w1) + ra_I_Rr2( r, r1, phi_1, w2) \
                + ra_I_Rr2( r, r1, phi_2, w1) + ra_I_Rr2( r, r2, phi_1, w1) \
                - ra_I_Rr2( r, r1, phi_2, w2) - ra_I_Rr2( r, r2, phi_1, w2) \
                - ra_I_Rr2( r, r2, phi_2, w1) + ra_I_Rr2( r, r2, phi_2, w2)
        #
        I_Zr1 =   ra_I_Zr1( r, r1, phi_n, w1) - ra_I_Zr1( r, r1, phi_n, w2) \
                - ra_I_Zr1( r, r2, phi_n, w1) + ra_I_Zr1( r, r2, phi_n, w2)
        
        I_Zr2 =   ra_I_Zr2( r, r_ND, r_nd, phi_ND, phi_nd, w2) \
                - ra_I_Zr2( r, r_ND, r_nd, phi_ND, phi_nd, w1)

        # magnetic field calculation
        Hx = cos(phy)* I_Rr1 - sin(phy)* I_Rr2
        Hy = cos(phy)* I_Rr2 + sin(phy)* I_Rr2
        Hz = -I_Zr1 + I_Zr2
        
        B_i = (µ* J[i]/ (4* pi))* array([ Hx, Hy, Hz ])
    
        # Transformation 2 and Addition of magnetic field
        B = B + dot( operator[i], B_i)
    
    return B

#%% Cylindrical wire
    
def cy_I_r2(r, r_n, ro, fai_j, wk):
    if r == 0:
        I_r2 = - cos(fai_j)* wk* arcsinh(ro/ wk)
    else:
        I_r2 = - (1/ r)* trapz( r_n* log( wk + sqrt(r**2 + r_n**2 + wk**2 - 2* r* r_n* cos(fai_j) )), r_n)
    return I_r2

def cy_I_r13(r, r_ND, r_nd, fai_ND, fai_nd, wk):
    return wk* trapz( trapz( ( r* r_ND - (r_ND**2)* cos(fai_ND) ) / 
                           ( ( r**2 + r_ND**2 - 2* r* r_ND* cos(fai_ND) ) * 
                           sqrt( r**2 + r_ND**2 + wk**2 - 2* r* r_ND* cos(fai_ND) ) ), r_nd), fai_nd)

#
def B_cy(µ, Nmbr, ro, l, J, operator, center, point):
    
    B = 0
    for i in range(0, Nmbr):
        pnt = dot( inv(operator[i]) , point - center[i] )
        #
        r, phy, z = Cart2Pol( pnt[0], pnt[1], pnt[2] )
        
        # constants
        fai_1 =-pi - phy
        fai_2 = pi - phy
        w1    =-(l[i]/2) - z
        w2    = (l[i]/2) - z
    
        # Integrals        
        r_n  = arange( 0, ro[i] + 0.001, 0.001 )
        r_nd = arange( 0, ro[i] + 0.01, 0.01 )
        fai_nd = arange( fai_1, fai_2 + 0.01, 0.01)
        
        r_ND, fai_ND = meshgrid(r_nd, fai_nd)
        #
        I_r13 = cy_I_r13(r, r_ND, r_nd, fai_ND, fai_nd, w2) - cy_I_r13(r, r_ND, r_nd, fai_ND, fai_nd, w1)
        
        I_r2 =  cy_I_r2(r, r_n, ro[i], fai_1, w1) - cy_I_r2(r, r_n, ro[i], fai_1, w2) \
              - cy_I_r2(r, r_n, ro[i], fai_2, w1) + cy_I_r2(r, r_n, ro[i], fai_2, w2)
    
        # Magnetic field calculation
        Hx = cos(phy)* I_r2 - sin(phy)* I_r13
        Hy = sin(phy)* I_r2 + cos(phy)* I_r13
        Hz = 0
        
        B_i = (µ* J[i]/ (4*pi) )* array([ Hx, Hy, Hz ])
    
        # Transformation 2 and Addition of magnetic field
        B = B + dot( operator[i], B_i)

    return B
  
#%% Helical wire
    
def B_Helical(µ, Nmbr, a, L, p, I, operator, center, point):

    B = 0
    for i in range(0, Nmbr):
        pnt = dot( inv(operator[i]), point - center[i] )
        
        # Test point
        x = pnt[0]
        y = pnt[1]
        z = pnt[2]
        
        # Magnetic field calculation 
        phi_n = arange(0, 2* pi* L[i]/ p[i] + pi/100, pi/100)
        D = (x - a[i]* cos(phi_n))**2 + (y - a[i]* sin(phi_n))**2 + (z - p[i]* phi_n/ (2* pi))**2

        Hx = simps( ( a[i]*cos(phi_n)* (z - p[i]* phi_n/ (2*pi)) - (y - a[i]* sin(phi_n)) * p[i]/ (2*pi) )/ D**1.5, phi_n)
        Hy = simps( ( a[i]*sin(phi_n)* (z - p[i]* phi_n/ (2*pi)) - (x - a[i]* cos(phi_n)) * p[i]/ (2*pi) )/ D**1.5, phi_n)
        Hz = simps( ( a[i]**2 - a[i]* ( x* cos(phi_n) + y* sin(phi_n) ) )/ D**1.5, phi_n)

        B_i = (µ* I[i]/ (4*pi))* array([ Hx, Hy, Hz ])

        # Transformation 2 and Addition of magnetic field
        B = B + dot( operator[i], B_i)
    
    return B
