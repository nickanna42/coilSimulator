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
    return endPoint - startPoint
    

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



def b_fromWireSegments(wireMiddle, wireTip, currentVector, evalPoint):
    """
    wireMiddle:      (n,3)np.float64
    wireTip:         (n,3)np.float64
    currentVector:   (n,3)np.float64
    evalPoint:       (1,3)np.float64

    returns: (n,3)np.float64
    -----------
    This function outputs the magnetic feild at evalPoint of a number
    of straight wire segments centered at wireMiddle, with the end
    the current is flowing towards centered at wireTip. The magnitude
    and 3D orientation in the wire is given by currentVector
    """
    P_x, P_y, P_z = distance(wireMiddle, evalPoint)
    I
    halfL_x, halfL_y, halfL_z = displace(wireMiddle, wireTip)
    theta
