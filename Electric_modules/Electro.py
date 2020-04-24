#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from scipy.integrate import simps
from numpy.linalg import inv
from numpy import dot, sqrt, pi, array, arange, meshgrid, sin, cos

#%% CAE Electric filament

def Electric_f(ε, Nmbr, l, lmbd, operator, center, point):
    
    E = 0
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
        
        Ex = - ( u/ (auv**2))* (-w1/ sqrt( auv**2 + w1**2) + w2/ sqrt( auv**2 + w2**2))
        Ey = - ( v/ (auv**2))* (-w1/ sqrt( auv**2 + w1**2) + w2/ sqrt( auv**2 + w2**2))
        Ez =   (-1/ sqrt( auv**2 + w1**2) + 1/ sqrt( auv**2 + w2**2))
        
        E_i = ( lmbd[i]/ (4*pi*ε))* array([ Ex, Ey, Ez ])
        
        # Transformation 2 and Addition of magnetic field
        E = E + dot( operator[i], E_i)
        
    return E

#%% CAE Electric Arc
    
def Electric_arc(ε, Nmbr, R, lmbd, operator, center, point):
    
    E = 0
    for i in range(0, Nmbr):
        pnt = dot( inv(operator[i]), point - center[i] )
        
        # phi is varing from 0 to 2 pi
        phi = arange(0, 2* pi + pi/100, pi/100)
        
        # Magnetic field calculation
        u   = pnt[0] - x
        v   = pnt[1] - y
        w1  = pnt[2] - z1
        w2  = pnt[2] - z2
        auv = sqrt(u**2 + v**2)
        
        Ex = - ( u/ (auv**2))* (-w1/ sqrt( auv**2 + w1**2) + w2/ sqrt( auv**2 + w2**2))
        Ey = - ( v/ (auv**2))* (-w1/ sqrt( auv**2 + w1**2) + w2/ sqrt( auv**2 + w2**2))
        Ez =   (-1/ sqrt( auv**2 + w1**2) + 1/ sqrt( auv**2 + w2**2))
        
        E_i = ( lmbd[i]/ (4*pi*ε))* array([ Ex, Ey, Ez ])
        
        # Transformation 2 and Addition of magnetic field
        E = E + dot( operator[i], E_i)
        
    return E
  
#%% CAE Electric Torus
    
def Electric_torus(ε, Nmbr, Ro, r, sgma, operator, center, point):
    
    E = 0
    for i in range(0, Nmbr):
        pnt = dot( inv(operator[i]), point - center[i] )
        
        tht1 = arange(0, 2* pi + pi/100, pi/100)
        phi1 = arange(0, 2* pi + pi/100, pi/100)
        
        tht, phi = meshgrid( tht1, phi1 )
        
        A = (pnt[0]**2 + pnt[1]**2 + pnt[2]**2 + Ro[i]**2 + r[i]**2)**0.5
        D = ( A**2 + 2* r[i]* (Ro[i]* cos(tht) - pnt[2]* sin(tht)) \
             + 2* (Ro[i] + r[i]* cos(tht))* (pnt[1]* sin(phi) - pnt[0]* cos(phi)) )**0.5
             
        X = r[i]* (Ro[i] + r[i]* cos(tht))* (pnt[0] - (Ro[i] + r[i]* cos(tht))* cos(phi))
        Y = r[i]* (Ro[i] + r[i]* cos(tht))* (pnt[1] + (Ro[i] + r[i]* cos(tht))* sin(phi))
        Z = r[i]* (Ro[i] + r[i]* cos(tht))* (pnt[2] - r[i]* sin(tht))
        
        Ex = simps( simps(X/ D**3, tht1), phi1)
        Ey = simps( simps(Y/ D**3, tht1), phi1)
        Ez = simps( simps(Z/ D**3, tht1), phi1)
        
        E_i = ( sgma[i]/ (4*pi*ε))* array([ Ex, Ey, Ez ])
        
        # Transformation 2 and Addition of magnetic field
        E = E + dot( operator[i], E_i)
        
    return E