"""
Original code by Nicholas Anna
released under GPLv2

list of functions in this file

b_fromWireSegments
circleLoop
displace
distance
makeCircleHHC
makeEvalPoints3D
makeSquareHHC
solve
squareLoop
"""

import numpy as np
#import matplotlib.pyplot as plt

def displace(startPoint, endPoint):
    """
    startPoint: (n,3)np.float64
    endPoint:   (n,3)np.float64

    returns: (n,3)np.float64
    ---------
    Gives the displacement vector between two sets of points.
    """
    out = np.array(endPoint - startPoint, ndmin=2, dtype="float64")
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
    return np.float64((a + b + c)**(.5))

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
    integral_f = constant * (b - a*t_f)/((b**2 - a*c)*(t_f*(a*t_f - 2*b)+c)**(0.5))
    integral_i = constant * (b - a*t_i)/((b**2 - a*c)*(t_i*(a*t_i - 2*b)+c)**(0.5))
    integralSolved = integral_f- integral_i
    preOut = np.empty([N,3], dtype="float64")
    preOut[:,0] = h_x*integralSolved
    preOut[:,1] = h_y*integralSolved
    preOut[:,2] = h_z*integralSolved
    preOut = np.sum(preOut, axis=0)
    out = np.array(preOut, ndmin=2)
    return out

def circleLoop(position, radius, orientVector, N, I):
    """
    Inputs-

    position:       (1,3)np.float64          Position in meters
    radius:              np.float64          Radius in meters
    orientVector:   (1,3)np.float64          Unit vector of loop axis
    N:                   int                 Number of wire segments
    I:                   np.float64          Current in amps
    -------------
    Return-

    positionVectors, lengthVectors, I_out

    positionVectors:     (n,3)np.float64    center of wire segments
    lengthVectors:       (n,3)np.float64    length vector
    I_out                (n)np.float64      current in amps
    ------------

    Creates an approximation of a single loop of current from N
    straight wire segments.
    """
    
    V_x, V_y, V_z = orientVector[0]
    phi_i = np.arctan2((V_x**2+V_y**2)**(0.5), V_z)
    theta_i = np.arctan2(V_y, V_x)
    cos_phi = np.cos(phi_i)
    sin_phi = np.sin(phi_i)
    cos_theta = np.cos(theta_i)
    sin_theta = np.sin(theta_i)

    segLength = 2*np.pi*radius/N
    
    A = np.arange(0, 2*np.pi, 2*np.pi/N)
    L = np.empty([N,3], dtype="float64")
    O = np.empty([N,3], dtype="float64")
    I_out = np.empty([N], dtype="float64")

    rHat = orientVector[0]
    phiHat = np.array([-cos_theta*cos_phi, -sin_theta*cos_phi, sin_phi])
    thetaHat = np.array([-sin_theta, cos_theta, 0.])

    C = np.empty([3,3], dtype="float64")
    C[:,0] = rHat
    C[:,1] = phiHat
    C[:,2] = thetaHat

    N_list = range(0, N)
    for n in N_list:
        L[n] = radius*np.array([np.cos(A[n]), np.sin(A[n]), 0.])
        O[n] = segLength*np.array([-np.sin(A[n]), np.cos(A[n]), 0.])
        I_out[n] = I
        L_temp = np.array(L[n], ndmin=2).transpose()
        O_temp = np.array(O[n], ndmin=2).transpose()
        L[n] = np.dot(C, L_temp)[:,0]
        O[n] = np.dot(C, O_temp)[:,0]

    L = L + position
    return L, O, I_out

def squareLoop(position, sideLength, orientVector, I):
    """
    position:            (1,3)np.float64     Position in meters
    sideLength:               np.float64     in meters
    orientVector:        (1,3)np.float64     unit vector of loop axis
    I:                        np.float64     current in amps
    -------------
    Return-

    positionVectors, lengthVectors, I_out

    positionVectors:     (4,3)np.float64    center of wire segments
    lengthVectors:       (4,3)np.float64    length vector
    I_out                (4)np.float64      current in amps
    ------------

    Creates a square coil of 4 wire segments
    """
    V_x, V_y, V_z = orientVector[0]
    phi_i = np.arctan2((V_x**2+V_y**2)**(0.5), V_z)
    theta_i = np.arctan2(V_y, V_x)
    cos_phi = np.cos(phi_i)
    sin_phi = np.sin(phi_i)
    cos_theta = np.cos(theta_i)
    sin_theta = np.sin(theta_i)
    
    rHat = orientVector[0]
    phiHat = np.array([-cos_theta*cos_phi, -sin_theta*cos_phi, sin_phi])
    thetaHat = np.array([-sin_theta, cos_theta, 0.])
    
    C = np.empty([3,3], dtype="float64")
    C[:,0] = rHat
    C[:,1] = phiHat
    C[:,2] = thetaHat
    
    A = np.arange(0, 2*np.pi, 2*np.pi/4)
    L = np.empty([4,3], dtype="float64")
    O = np.empty([4,3], dtype="float64")
    I_out = np.empty([4], dtype="float64")
    
    N_list = range(0, 4)
    for n in N_list:
        L[n] = (sideLength/2)*np.array([np.cos(A[n]), np.sin(A[n]), 0.])
        O[n] = sideLength*np.array([-np.sin(A[n]), np.cos(A[n]), 0.])
        I_out[n] = I
        L_temp = np.array(L[n], ndmin=2).transpose()
        O_temp = np.array(O[n], ndmin=2).transpose()
        L[n] = np.dot(C, L_temp)[:,0]
        O[n] = np.dot(C, O_temp)[:,0]
    
    L = L + position
    
    return L, O, I_out

def makeCircleHHC(coilCenter, coilOrientation, coilSpecs, R, I, N=10000):
    """
    coilCenter:         (1,3np.float64     position (m)
    coilOrientation:    (1,3)np.float64    unit vector. orientation of field

    coilSpecs            list              coilSpec is a list containing the
                                           following aggregate coil properties,
                                           in the order. wire diameter (m),
                                           # of lengthwise wraps, # of wraps
                                           deep.
                            
    R                   np.float64         radius in meters
    I                   np.float64         current in amps
    N                   int                number of line segments per wrap
                                           (default 10,000)
    ----------------
    Returns-
    
    segmentPosition, segmentLength, I_out
    
    segmentPosition    (n,3)np.float64
    segmentLength      (n,3)np.float64
    I_out              (n)np.float64
    ----------------
    Makes a helmholtz coil. The coil wraps are centered, coaxially and radially,
    around the coil wrap who's center is R/2 from coilCenter and has radius R.
    """
    segmentPosition = np.array([[0,0,0]], dtype="float64")
    segmentLength = np.array([[0,0,0]], dtype="float64")
    I_out = np.array([0], dtype="float64")
    tempL = coilCenter + (R/2)*coilOrientation + (coilSpecs[1]-1)*coilSpecs[0]*coilOrientation/2
    for i in range(0,coilSpecs[1]):
        tempR = R + (coilSpecs[2]-1)*coilSpecs[0]/2
        for j in range(0,coilSpecs[2]):
            temp_1, temp_2, temp_3 = circleLoop(tempL, tempR, coilOrientation, N, I)
            segmentPosition = np.append(segmentPosition, temp_1, axis=0)
            segmentLength = np.append(segmentLength, temp_2, axis=0)
            I_out = np.append(I_out, temp_3)
            
            temp_1, temp_2, temp_3 = circleLoop(-tempL, tempR, coilOrientation, N, I)
            segmentPosition = np.append(segmentPosition, temp_1, axis=0)
            segmentLength = np.append(segmentLength, temp_2, axis=0)
            I_out = np.append(I_out, temp_3)
            
            tempR = tempR - coilSpecs[0]
        tempL = tempL - coilOrientation*coilSpecs[0]
    segmentPosition = np.delete(segmentPosition, 0, axis=0)
    segmentLength = np.delete(segmentLength, 0, axis=0)
    I_out = np.delete(I_out, 0, axis=0)
    return segmentPosition, segmentLength, I_out

def makeSquareHHC(coilCenter, coilOrientation, coilSpecs, sideLength, I):
    """
    coilCenter:         (1,3np.float64     position (m)
    coilOrientation:    (1,3)np.float64    unit vector. orientation of field

    coilSpecs            list              coilSpec is a list containing the
                                           following aggregate coil properties,
                                           in the order. wire diameter (m),
                                           # of lengthwise wraps, # of wraps
                                           deep.
                            
    sideLength          np.float64         side length of coil in meters
    I                   np.float64         current in amps
    ----------------
    Returns-
    
    segmentPosition, segmentLength, I_out
    
    segmentPosition    (n,3)np.float64
    segmentLength      (n,3)np.float64
    I_out              (n)np.float64
    -------------
    makes a sqaure one
    """
    segmentPosition = np.array([[0,0,0]], dtype="float64")
    segmentLength = np.array([[0,0,0]], dtype="float64")
    I_out = np.array([], dtype="float64")
    tempL = coilCenter + sideLength*1.00*coilOrientation + (coilSpec[1]-1)*coilSpec[0]*coilOrientation/2
    for i in range(0,coilSpecs[1]):
        tempR = sideLength/2 + (coilSpec[2]-1)*coilSpec[0]/2
        for j in range(0,coilSpecs[2]):
            temp_1, temp_2, temp_3 = squareLoop(tempL, tempR, coilOrientation, I)
            segmentPosition = np.append(segmentPosition, temp_1, axis=0)
            segmentLength = np.append(segmentLength, temp_2, axis=0)
            I_out = np.append(I_out, temp_3)
            
            temp_1, temp_2, temp_3 = squareLoop(-tempL, tempR, coilOrientation, I)
            segmentPosition = np.append(segmentPosition, temp_1, axis=0)
            segmentLength = np.append(segmentLength, temp_2, axis=0)
            I_out = np.append(I_out, temp_3)
            
            tempR = tempR - coilSpecs[0]
        tempL = tempL - coilOrientation*coilSpecs[0]
    segmentPosition = np.delete(segmentPosition, 0, axis=0)
    segmentLength = np.delete(segmentLength, 0, axis=0)
    I_out = np.delete(I_out, 0, axis=0)
    return segmentPosition, segmentLength, I_out

def makeEvalPoints3D(postion, size_xyz, steps_xyz):
    """
    position:      (1,3)np.float64
    size_xyz:      (1,3)np.float64     must be > 0
    steps_xyz:     (1,3)np.float64     must be > 1
    
    ------------
    Returns-
    
    pointsArray:     (n,3)np.float64
    ------------
    Returns an array of points, which are arranged in an orthogonal grid.
    This grid is centered at position, is size_xyz, and goes 
    """
    X = size_xyz[0,0]
    Y = size_xyz[0,1]
    Z = size_xyz[0,2]
    I = steps_xyz[0,0]+1
    J = steps_xyz[0,1]+1
    K = steps_xyz[0,2]+1
    N = I*J*K
    n = 0
    pointsArray = np.empty([N,3], dtype="float64")
    for i in range(0,I):
        for j in range(0,J):
            for k in range(0,K):
                pointsArray[n] = np.array([X/2 - i*X/(I-1), Y/2 - j*Y/(J-1), Z/2 - k*Z/(K-1)])
                n = n + 1
    
    return pointsArray

def solve(segmentPosition, segmentLength, segmentI, evalPoints):
    """
    segmentPosition:    (n_1, 3)np.float64
    segmentLength:      (n_1, 3)np.float64
    segmentI:              (n_1)np.float64
    evalPoints          (n_2, 3)np.float64
    ------------
    Returns-
    
    magneticFeild:      (n_2, 3)np.float64
    """
    points_N = len(evalPoints)
    magneticFeild = np.zeros( [points_N, 3], dtype="float64")
    for n in range(0, points_N):
        magneticFeild[n] = b_fromWireSegments(segmentPosition, segmentLength, segmentI, np.array([evalPoints[n]]))
    return magneticFeild