#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from scipy.integrate import simps
from numpy.linalg import inv
from numpy import pi, sqrt, sin, cos, array, dot, arctan2, arange

def Cart2Pol(x, y, z):
    rho = sqrt(x**2 + y**2)
    phi = arctan2(y, x)
    return( rho, phi, z)

#%% Arc1

from scipy.special import ellipkinc as K_elp
from scipy.special import ellipeinc as E_elp

''' Jacobian elliptic functions
    ---> Amplitude: t | Parameter: m
    ---> elliptic modulus: k = sqrt(m) | Modular angle : m = sin(angle)**2
    ---> Complements of k and m:
         m' = 1 - m | k' = sqrt(m') '''

sn_u  = lambda t : sin(t)
cn_u  = lambda t : cos(t)
nd_u  = lambda t, m : (1 - m* sin(t)**2 )**(-0.5)
cd_u  = lambda t, m : cn_u(t)* nd_u(t, m)
W_elp = lambda t, m : E_elp(t, m) - m* sn_u(t)* cd_u(t, m)

sgn = lambda x : 1 if x >= 0 else -1
mod = lambda x : sgn(x)* x

def B_a1(myeu, Nmbr, r1, phy_1, phy_2, I, operator, cnt_cart, point):
    
    # Expressions for H
    Hr_exp = lambda t, alp : gma* sgn(alp)* (kp**2 * K_elp(t, k**2) - (1 - k**2/2)* W_elp(t, k**2) )
    Hp_exp = lambda t : -gma* kp**2 * nd_u(t, k**2) 
    Hz_exp = lambda t, alp : sgn(alp)* (r* kp**2 * K_elp(t, k**2) - (r - b* k**2/2)* W_elp(t, k**2) )
    
    Hr_exp2 = lambda t, alp : gma* sgn(alp)* (kp**2 * ( 2* K_elp(pi/2, k**2) - K_elp(t, k**2) ) - (1 - k**2/2)* ( 2* W_elp(pi/2, k**2) - W_elp(t, k**2) ) )
    Hz_exp2 = lambda t, alp : sgn(alp)* (r* kp**2 * ( 2* K_elp(pi/2, k**2) - K_elp(t, k**2) ) - (r - b* k**2/2)* ( 2* W_elp(pi/2, k**2) - W_elp(t, k**2) )  )
    
    B = 0
    for i in range(0, Nmbr):
        pnt = dot( inv(operator[i]) , point - cnt_cart[i] )
        #
        r, phy, z = Cart2Pol( pnt[0], pnt[1], pnt[2] )
        phy = phy if (phy >= 0) else (phy + 2* pi)
        
        if r == 0:
            # Constants
            b = sqrt( r1[i]**2 + r**2 + z**2 )
            phi_1 = phy_1[i] - phy
            phi_2 = phy_2[i] - phy
            
            # Integrals
            phi_n = arange( phi_1, phi_2 + 0.001, 0.001)
            #
            In_1 = simps(cos(phi_n)* (b**2 - 2* r* r1[i]* cos(phi_n) )**(-1.5) , phi_n)
            In_2 = ( r1[i]**2 + z**2 )**(-1.5) * ( cos(phi_1) - cos(phi_2) )
            In_3 = simps( (b**2 - 2* r* r1[i]* cos(phi_n) )**(-1.5) , phi_n)

            # Magnetic field calculation
            Hx = r1[i]* z* ( cos(phy)* In_1 - sin(phy)* In_2)
            Hy = r1[i]* z* ( sin(phy)* In_1 + cos(phy)* In_2)
            Hz = r1[i]* (r1[i]* In_3 - r* In_1)

        else:
            # Constants
            gma = 0 - z
            a = sqrt( gma**2 + (r1[i] + r)**2 )
            b = r1[i] + r
            k = sqrt( 4* r* r1[i]/ a**2 )
            kp = sqrt( 1 - k**2 )
            alp1 = (pi + phy - phy_1[i])/2
            alp2 = (pi + phy - phy_2[i])/2
            
            if (phy_1 <= phy ) and (phy <= phy_2 ):
                
                tht1 = mod(alp1)
                tht2 = mod(alp2)
                
                Hr = (1/(r*a*kp**2))* ( Hr_exp(tht2, alp2) - Hr_exp(tht1, alp1) )
                Hp = (1/(r*a*kp**2))* ( Hp_exp(tht2) - Hp_exp(tht1) )
                Hz = (1/(r*a*kp**2))* ( Hz_exp(tht2, alp2) - Hz_exp(tht1, alp1) )
                
            elif (phy_1 > phy ) or (phy > phy_2 ):
                
                tht1 = pi - mod(alp1)
                tht2 = pi - mod(alp2)
#                alpp = (pi/2 + phy)/2
                
#                Hr = (1/(r*a*kp**2))* (  sgn(alp2)* (2* Hr_exp(pi/2, alpp) - Hr_exp(tht2, alp2)) \
#                                       - sgn(alp1)* (2* Hr_exp(pi/2, alpp) - Hr_exp(tht1, alp1)) )
                Hr = (1/(r*a*kp**2))* ( Hr_exp2(tht2, alp2) - Hr_exp2(tht1, alp1) )
                Hp = (1/(r*a*kp**2))* ( Hp_exp(tht2) - Hp_exp(tht1) )
                Hz = (1/(r*a*kp**2))* ( Hz_exp2(tht2, alp2) - Hz_exp2(tht1, alp1) )
#                Hz = (1/(r*a*kp**2))* (  sgn(alp2)* (2* Hz_exp(pi/2, alpp) - Hz_exp(tht2, alp2)) \
#                                       - sgn(alp1)* (2* Hz_exp(pi/2, alpp) - Hz_exp(tht1, alp1)) )

            # Cartesian field
            Hx =-( Hr* cos(phy) + Hp* sin(phy) )
            Hy =- Hr* sin(phy) + Hp* cos(phy)
            Hz =- Hz

        # Magnetic field calculation
        B_i = (myeu* I[i]/ (4*pi))* array([ Hx, Hy, Hz ])
    
        # Transformation 2 and Addition of magnetic field
        B = B + dot(operator[i] , B_i)

    return B
