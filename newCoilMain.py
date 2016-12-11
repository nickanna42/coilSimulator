"""
Original code by Nicholas Anna
Released under GPL v2
----------
Release Notes.

mayavi is only available on python 2.7 as of most recent update.
make sure to run this code on a 2.7 kernel

----------
Magnetic Coil Simulator

Creates a set of arrays which represent a series of straight wires with current.
Creates an array which represents series of points in 3D space.
Solves for the magnetic feild at a point in space, given a series of wires with current
Displays the magnetic field over space.

The code is distrubuted set up to display a create and display the field from a simple
Helmholtz coil.
"""

import numpy as np
import matplotlib.pyplot as plt
import newCoilFunct as cf
from mayavi import mlab

coilCenter = np.array([[0.,0.,0.]], dtype="float64")
coilRadius = np.float64(1)
coilOrient = np.array([[1.,0.,0.]], dtype='float64')
coilSpec = [np.float64(.002), 1, 1]

gridCenter = np.array([[0.,0.,0.]], dtype ="float64")
gridSize = np.array([[2.,2.,2.]], dtype ="float64")
gridResolution = np.array([[5,5,5]])

evalGrid = cf.makeEvalPoints3D(gridCenter, gridSize, gridResolution)

segmentLocation, segmentLength, segmentI = cf.makeCircleHHC(coilCenter, coilOrient, coilSpec, coilRadius, np.float64(60), 10000)

magneticGrid = cf.solve(segmentLocation, segmentLength, segmentI, evalGrid)

mlab.quiver3d(evalGrid[:,0], evalGrid[:,1], evalGrid[:,2], magneticGrid[:,0], magneticGrid[:,1], magneticGrid[:,2])