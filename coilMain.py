import numpy as np
from coilFunct import *
from mayavi import mlab
# defines the size of the area to be evaluated and specifies
# the number of points per side to evaluate
xWidth = .30
xCount = 3
yWidth = .30
yCount = 3
zWidth = .30
zCount = 3
nPoints = xCount * yCount * zCount #number of points where B is evaluated
#
# This section defines the characterstics of the wire guide
hX = .884 # gap between the two halfes of the coil set, x-axis
aX = 1.626 # length of one side of the wire guide, x-axis
hY = 1.071 # gap between the two halfes of the coil set, y-axis
aY = 1.975 # length of one side of the wire guide, y-axis
hZ = .36 # gap between the two halfes of the coil, z-axis
aZ = .36 # radius of the wire guide, z-axis
#
# defines the coil bundles by wraps per layer & # of layers
layerX = 0
thickX = 0
layerY = 0
thickY = 0
layerZ = 1
thickZ = 1
#
# the diameter of the wire used on each set of coils
xDia = 0.0008788
yDia = 0.001095
zDia = 0.0008788
#
# current in each set of coils
xI = 0
yI = 0
zI = 36
#
#creates an array of all points at which B will be evaluated.
evalPoints = evalList(xWidth, xCount, yWidth, yCount, zWidth, zCount)
#
#sums the feild from all wraps on the x-axis coil pair
#at all evaluation points
fieldPoints = np.zeros([nPoints, 3])
n = 0
while n < nPoints:
    i = 0    
    while i < layerX:
        j = 0        
        while j < thickX:
            fieldPoints[n] = fieldPoints[n] + sqaureLoopX(evalPoints[n], aX + xDia/2 + j*xDia, hX + xDia/2 + i*xDia, 5000, xI)
            fieldPoints[n] = fieldPoints[n] + sqaureLoopX(evalPoints[n], aX + xDia/2 + j*xDia, -(hX + xDia/2 + i*xDia), 5000, xI)      
            j = j + 1
        i = i + 1
    n = n + 1
#
#
n = 0
while n < nPoints:
    i = 0    
    while i < layerY:
        j = 0        
        while j < thickY:
            fieldPoints[n] = fieldPoints[n] + sqaureLoopY(evalPoints[n], aY + yDia/2 + j*yDia, hY + yDia/2 + i*xDia, 5000, yI)
            fieldPoints[n] = fieldPoints[n] + sqaureLoopY(evalPoints[n], aY + yDia/2 + j*yDia, -(hY + yDia/2 + i*xDia), 5000, yI)
            j = j + 1
        i = i + 1
    n = n + 1
#
#
n = 0
while n < nPoints:
    i = 0
    while i < layerZ:
        j = 0
        while j < thickZ:
            fieldPoints[n] = fieldPoints[n] + circleLoopZ(evalPoints[n], aZ + zDia/2 + j*zDia, hZ + zDia/2 + i*zDia, 10000, zI)
            fieldPoints[n] = fieldPoints[n] + circleLoopZ(evalPoints[n], aZ + zDia/2 + j*zDia, -(hZ + zDia/2 + i*zDia), 10000, zI)
            j = j + 1
        i = i + 1
    n = n + 1
#
#This section converts our results to a format readable by quiver3d
n = 0
x = np.empty([nPoints, 1])
while n < nPoints:
    x[n] = evalPoints[n,0]
    n = n + 1
n = 0
y = np.empty([nPoints, 1])
while n < nPoints:
    y[n] = evalPoints[n, 1]
    n = n + 1
n = 0
z = np.empty([nPoints, 1])
while n < nPoints:
    z[n] = evalPoints[n, 2]
    n = n + 1
n = 0
u = np.empty([nPoints, 1])
while n < nPoints:
    u[n] = fieldPoints[n, 0]
    n = n + 1
n = 0
v = np.empty([nPoints, 1])
while n < nPoints:
    v[n] = fieldPoints[n, 1]
    n = n + 1
n = 0
w = np.empty([nPoints, 1])
while n < nPoints:
    w[n] = fieldPoints[n, 2]
    n = n + 1
#
# displays results
mlab.quiver3d(x, y, z, u, v, w)