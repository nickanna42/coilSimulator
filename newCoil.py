"""
Original code by Nicholas Anna
released under GPLv2
"""

import numpy as np

def displace(startPoint, endPoint):
    """
    startPoint: (n,3)np.float64
    endPoint:   (n,3)np.float64

    returns: (n,3)np.float64
    ---------
    Gives the displacement vector between two points.
    """
    out = np.array(endPoint - startPoint, ndmin=2, dtype=float64)
    return out

def distance(startPoint, endPoint):
    """
    startPoint: (n,3)np.float64
    endPoint:   (n,3)np.float64

    returns: (n)np.float64
    ---------
    Gives the distance scalar between two points
    """
    displaceVector = displace(startPoint, endPoint)
    a = displaceVector[:,0]**2
    b = displaceVector[:,1]**2
    c = displaceVector[:,2]**2
    return np.float64((a + b + c)**(.5))2



def b_fromWireSegments(wireMiddle, wireL, current, evalPoint):
    """
    wireMiddle:    (n,3)np.float64
    wireL:         (n,3)np.float64
    current:       (n)np.float64
    evalPoint:     (1,3)np.float64

    returns: (1,3)np.float64
    -----------
    This function outputs the magnetic feild at evalPoint from a
    number of straight wire segments. The wires segments are
    centered at wireMiddle, with current flowing in the direction
    of the length vector.
    """
    N = len(current)
    constant = np.pi*1e-7*current
    P = displace(wireMiddle, evalPoint)
    P_x = P[:,0]
    P_y = P[:,1]
    P_z = P[:,2]
    L_x = wireL[:,0]
    L_y = wireL[:,1]
    L_z = wireL[:,2]
    h_x = P_z*L_y - P_y*L_z
    h_y = -(P_z*L_x - P_x*L_z)
    h_z = P_y*L_x - P_x*L_y
    a = L_x**2 + L_y**2 + L_z**2
    b = L_x*P_x + L_y*P_y + L_z*P_z
    c = P_x**2 + P_y**2 + P_z**2
    t_i, t_f = np.float64(-0.5), np.float64(0.5)
    integral_f = constant * (b - a*t_f)/((b**2 - a*c)(t_f*(a*t_f - 2*b)+c)**(0.5))
    integral_i = constant * (b - a*t_i)/((b**2 - a*c)(t_i*(a*t_i - 2*b)+c)**(0.5))
    integralSolved = integral_f- integral_i
    preOut = np.empty([N,3] dtype=float64)
    preOut[:,0] = h_x*integralSolved
    preOut[:,1] = h_y*integralSolved
    preOut[:,2] = h_z*integralSolved
    out = np.sum(preOut, axis=0, ndmin=2)
    return out

def circleLoop(position, orientVector, radius, N, I):
    """
    Inputs-

    position:       (1,3)np.float64          Position in meters
    orientVector:   (1,3)np.float64          Unit vector of loop axis
    radius:              np.float64          Radius in meters
    N:                   int                 Number of wire segments
    I:                   np.float64          Current in amps
    -------------
    Return-

    positionVectors, lengthVectors

    positionVectors:     (n,3)np.float64    center of wire segments
    lengthVectors:       (n,3)np.float64    length vector
    ------------

    Creates an approximation of a single loop of current from N
    straight wire segments.
    """
    
    V_x, V_y, V_z = orientVector[0]
    phi_i = np.arctan2((V_x**2+V-y**2)**(0.5), V_z)
    theta_i = np.arctan2(V_y, V_x)
    cos_phi = np.cos(phi_i)
    sin_phi = np.sin(phi_i)
    cos_theta = np.cos(theta_i)
    sin_theta = np.sin(theta_i)

    segLength = 2*np.pi*radius/N
    
    A = np.arange(0, 2*np.pi, 1/N)
    L = np.empty([N,3}, ndmin=2)
    O = np.empty([N,3], ndmin=2)

    rHat = orientVector[0]
    phiHat = np.array([-cos_theta*cos_phi, -sin_theta*cos_phi, sin_phi])
    thetaHat = np.array([-sin_theta, cos_theta, 0.])

    C = np.empty([3,3], dtype=float64)
    C[:,0] = rHat
    C[:,1] = phiHat
    C[:,2] = thetaHat
    
    for n in N:
        L[n] = radius*np.array([np.cos(A[n]), np.sin(A[n]), 0.])
        O[n] = segLength*np.array([-np.sin(A[n]), np.cos(A[n]), 0.])
        L_temp = np.array(L[n], ndmin=2).transpose()
        O_temp = np.array(O[n], ndmin=2).transpose()
        L[n] = np.dot(C, L_temp)[:,0]
        O[n] = np.dot(C, O_temp)[:,0]

    L = L + position
    O = O + position
    return L, O
    
