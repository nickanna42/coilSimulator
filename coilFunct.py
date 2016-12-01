import numpy as np
def evalList(xWidth, xCount, yWidth, yCount, zWidth, zCount):
    """Returns a nested array of [x,y,z] points within the specified space"""
    x, xStep = -xWidth/2, xWidth / (xCount - 1)
    y, yStep = -yWidth/2, yWidth / (yCount - 1)
    z, zStep = -zWidth/2, zWidth / (zCount - 1)
    evalPoints = np.empty([(xCount * yCount * zCount), 3])
    n = 0
    while x <= xWidth/2:
        y = -yWidth/2
        while y <= yWidth/2:
            z = -zWidth/2
            while z <= zWidth/2:
                evalPoints[n] = [x,y,z]
                z = z + zStep
                n = n + 1
            y = y + yStep
        x = x + xStep
    return evalPoints
#
def sqaureLoopZ(evalPoint, a, h, n, I):
    xPrime = evalPoint[0]
    yPrime = evalPoint[1]
    zPrime = evalPoint[2]
    dl = a * 4 / n
    x = a/2
    y = -a/2
    z = h/2
    mu = (4 * np.pi)*10**(-7)
    pointField = np.array([0, 0, 0])
    while y <= a/2:
        r = np.array([xPrime - x, yPrime - y, zPrime - z])
        pointField = np.add(pointField, np.cross( np.array([0, dl, 0]), r) * mu * I / (4 * np.pi * np.linalg.norm(r)**3))
        y = y + dl
    while x >= -a/2:
        r = np.array([xPrime - x, yPrime - y, zPrime - z])
        pointField = np.add(pointField, np.cross( np.array([-dl, 0, 0]), r) * mu * I / (4 * np.pi * np.linalg.norm(r)**3))
        x = x - dl
    while y >= -a/2:
        r = np.array([xPrime - x, yPrime - y, zPrime - z])
        pointField = np.add(pointField, np.cross( np.array([0, -dl, 0]), r) * mu * I / (4 * np.pi * np.linalg.norm(r)**3))
        y = y - dl
    while x <= a/2:
        r = np.array([xPrime - x, yPrime - y, zPrime - z])
        pointField = np.add(pointField, np.cross( np.array([dl, 0, 0]), r) * mu * I / (4 * np.pi * np.linalg.norm(r)**3))
        x = x + dl
    return pointField
#
def sqaureLoopX(evalPoint, a, h, n, I):
    xPrime = evalPoint[0]
    yPrime = evalPoint[1]
    zPrime = evalPoint[2]
    dl = a*4/n
    mu = (4 * np.pi)*10**(-7)
    x = h/2
    y = -a/2
    z = a/2
    pointField = np.array([0, 0, 0])
    while z >= -a/2:
        r = np.array([xPrime - x, yPrime - y, zPrime - z])
        pointField = np.add(pointField, np.cross( np.array([0, 0, -dl]), r) * mu * I / (4 * np.pi * np.linalg.norm(r)**3))
        z = z - dl
    while y <= a/2:
        r = np.array([xPrime - x, yPrime - y, zPrime - z])
        pointField = np.add(pointField, np.cross( np.array([0, dl, 0]), r) * mu * I / (4 * np.pi * np.linalg.norm(r)**3))
        y = y + dl
    while z <= a/2:
        r = np.array([xPrime - x, yPrime - y, zPrime - z])
        pointField = np.add(pointField, np.cross( np.array([0, 0, dl]), r) * mu * I / (4 * np.pi * np.linalg.norm(r)**3))
        z = z + dl
    while y >= -a/2:
        r = np.array([xPrime - x, yPrime - y, zPrime - z])
        pointField = np.add(pointField, np.cross( np.array([0, -dl, 0]), r) * mu * I / (4 * np.pi * np.linalg.norm(r)**3))
        y = y - dl
    return pointField
#
def sqaureLoopY(evalPoint, a, h, n, I):
    xPrime = evalPoint[0]
    yPrime = evalPoint[1]
    zPrime = evalPoint[2]
    dl = a*4/n
    mu = (4 * np.pi)*10**(-7)
    x = a/2
    y = h/2
    z = a/2
    pointField = np.array([0, 0, 0])
    while z >= -a/2:
        r = np.array([xPrime - x, yPrime - y, zPrime - z])
        pointField = np.add(pointField, np.cross( np.array([0, 0, -dl]), r) * mu * I / (4 * np.pi * np.linalg.norm(r)**3))
        z = z - dl
    while x >= -a/2:
        r = np.array([xPrime - x, yPrime - y, zPrime - z])
        pointField = np.add(pointField, np.cross( np.array([-dl, 0, 0]), r) * mu * I / (4 * np.pi * np.linalg.norm(r)**3))
        x = x - dl
    while z <= a/2:
        r = np.array([xPrime - x, yPrime - y, zPrime - z])
        pointField = np.add(pointField, np.cross( np.array([0, 0, dl]), r) * mu * I / (4 * np.pi * np.linalg.norm(r)**3))
        z = z + dl
    while x <= a/2:
        r = np.array([xPrime - x, yPrime - y, zPrime - z])
        pointField = np.add(pointField, np.cross( np.array([dl, 0, 0]), r) * mu * I / (4 * np.pi * np.linalg.norm(r)**3))
        x = x + dl
    return pointField
#
def circleLoopZ(evalPoint, a, h, n, I):    
    xPrime = evalPoint[0]
    yPrime = evalPoint[1]
    zPrime = evalPoint[2]
    dPhi = (2 * np.pi) / n
    mu = (4 * np.pi)*10**(-7)
    phi = 0
    pointField = np.array([0, 0, 0])    
    while phi < 2 * np.pi:
        x = -np.sin(phi)*a
        y = np.cos(phi)*a
        z = h/2
        r = np.array([xPrime - x, yPrime - y, zPrime - z])
        dl = np.array([-np.cos(phi) * dPhi * a, -np.sin(phi) * dPhi * a, 0])
        pointField = np.add(pointField, np.cross(dl, r ) * mu * I / (4 * np.pi * np.linalg.norm(r)**3))
        phi = phi + dPhi
    return pointField
#
def circleLoopX(evalPoint, a, h, n, I):
    xPrime = evalPoint[0]
    yPrime = evalPoint[1]
    zPrime = evalPoint[2]
    dPhi = (2 * np.pi) / n
    mu = (4 * np.pi)*10**(-7)
    phi = 0
    pointField = np.array([0, 0, 0])
    while phi < 2 * np.pi:
        x = h/2
        y = -np.sin(phi) * a
        z = np.cos(phi) * a
        r = np.array([xPrime - x, yPrime - y, zPrime - z])
        dl = np.array([0, -np.cos(phi) * dPhi * a, -np.sin(phi) * dPhi * a])
        pointField = np.add(pointField, np.cross(dl, r) * mu * I / (4 * np.pi *np.linalg.norm(r)**3))
        phi = phi + dPhi
    return pointField
#
def circleLoopY(evalPoint, a, h, n, I):
    xPrime = evalPoint[0]
    yPrime = evalPoint[1]
    zPrime = evalPoint[2]
    dPhi = (2 * np.pi) / n
    mu = 1.25663706 * 10**(-6)
    phi = 0
    pointField = np.array([0, 0, 0])
    while phi < 2 * np.pi:
        x = np.sin(phi) * a
        y = h / 2
        z = np.cos(phi) * a
        r = np.array([xPrime - x, yPrime - y, zPrime - z])
        dl = np.array([np.cos(phi) * dPhi * a, 0, -np.sin(phi) * dPhi *a])
        pointField = np.add(pointField, np.cross(dl, r) * mu * I / (4 * np.pi *np.linalg.norm(r)**3))
        phi = phi + dPhi
    return pointField