#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from numpy import sin, cos, sqrt, arctan, pi, array, arcsin, dot              
π = pi

def Rot_xy(A):
    
    D = sqrt(A[0]**2 + A[1]**2 + A[2]**2)
    
    Ψ = 0
    θ =-arcsin(A[0]/ D)
    
    if A[1] != 0:
        if   A[2] == 0:
             Ψ = π/2
        elif A[2] >  0:
             Ψ = arctan(A[1]/ A[2])
        elif A[2] <  0:
             Ψ = arctan(A[1]/ A[2]) + π
    
    R = array([ [ cos(θ), sin(Ψ)* sin(θ), cos(Ψ)* sin(θ)],
                [     0 ,         cos(Ψ),        -sin(Ψ)],
                [-sin(θ), sin(Ψ)* cos(θ), cos(Ψ)* cos(θ)] ])
    
    B = dot(R, A)
    
    return R, B

#%%
 
def FERSA(m, q, ro, Vo, E, Bz, t):
    ''' Field Euler Rotatonal Symmetry Algorithm'''
    
    xo = ro[0]
    yo = ro[1]
    zo = ro[2]
    
    Vxo = Vo[0]
    Vyo = Vo[1]
    Vzo = Vo[2]
    
    Ex = E[0]
    Ey = E[1]
    Ez = E[2]

    # Can be plus or minus
    ω = q* Bz/ m
    
    # Integration constants
    A = sqrt( (Vxo - Ey/Bz)**2 + (Vyo + Ex/Bz)**2 )

    if   (Vyo + Ex/Bz) == 0:
         Ф = π/2
    elif (Vyo + Ex/Bz) >  0:
         Ф = arctan( (Vxo - Ey/Bz)/ (Vyo + Ex/Bz) )
    elif (Vyo + Ex/Bz) <  0:
         Ф = arctan( (Vxo - Ey/Bz)/ (Vyo + Ex/Bz) ) + π
    
    # Velocity
    Vx = A* sin(ω* t + Ф) + Ey/ Bz
    Vy = A* cos(ω* t + Ф) - Ex/ Bz
    Vz = Vzo + (ω* Ez/Bz)* t
    
    # Guiding centre
    xc = xo + (A/ ω)* cos(Ф)
    yc = yo - (A/ ω)* sin(Ф)
    
    # Position
    x = xc + (Ey/ Bz)* t - (A/ ω)* cos(ω* t + Ф)
    y = yc - (Ex/ Bz)* t + (A/ ω)* sin(ω* t + Ф)
    z = zo +      Vzo* t + 0.5* (ω* Ez/Bz)* t**2
    
    # Position and Velocity
    r = array([x,   y,  z])
    V = array([Vx, Vy, Vz])
    
    return r, V, Ф
