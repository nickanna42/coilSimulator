# -*- coding: utf-8 -*-
"""
Created on Sat Nov  8 09:23:29 2014

@author: nick
"""
import numpy as np
from mayavi import mlab
pointList = np.empty([40, 3])
vectorList = np.empty([40, 3])
dPhi = 2 * np.pi / 40
n = 0
phi = 0
while phi < 2 * np.pi:
    pointList[n] = [np.cos(phi)*20, np.sin(phi)*20, 0]
    phi = phi + dPhi
    n = n + 1
n = 0
phi = 0
while phi < 2 * np.pi:
    vectorList[n] = [-np.sin(phi)*4, np.cos(phi)*4, 4]
    phi = phi + dPhi
    n = n + 1
x = np.empty([40, 1])
n = 0
while n < 40:
    x[n] = pointList[n,0]
    n  = n + 1
y = np.empty([40, 1])
n = 0
while n < 40:
    y[n] = pointList[n, 1]
    n = n + 1
z = np.empty([40, 1])
n = 0
while n < 40:
    z[n] = pointList[n, 2]
    n = n + 1
u = np.empty([40,1])
n = 0
while n < 40:
    u[n] = vectorList[n, 0]
    n = n + 1
v = np.empty([40,1])
n = 0
while n < 40:
    v[n] = vectorList[n, 1]
    n = n + 1
w = np.empty([40, 1])
n = 0
while n < 40:
    w[n] = vectorList[n, 2]
    n = n + 1
mlab.quiver3d(x, y, z, u, v, w)