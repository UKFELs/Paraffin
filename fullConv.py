
# Copyright 2012-2018, University of Strathclyde
# Authors: Lawrence T. Campbell and Piotr Traczykowski
# License: BSD-3-Clause

"""
Simple example script for using different packages in Paraffin.
"""

import sys
import numpy as np
from paraffin.matching import SU2Matched  # for beam matching
from paraffin.puffin import puffData      # For setting up the scaled Puffin frame
from paraffin.puffin import undulator     # For specifying the undulator for matching
from paraffin.upsample import CDFupsampler  # For upsampling the beam to a finer resolution
from paraffin.SU2Puffin import SU2Puffin  # For scaling the SU particle data to Puffin units

##### EXAMPLE DRIVER:
def fullConv(fnamein):

    puffVars = puffData()

# Setup the Puffin scaled frame:
#     This example is for CLARA

    puffVars.aw = 0.8667*np.sqrt(2)   # The PEAK!!!
    puffVars.gamma0 = 469.0
    puffVars.lw = 0.025
    puffVars.rho = 0.005
    puffVars.undtype = ''
    puffVars.ux = 1.
    puffVars.uy = 0.

# Generate the rest of the Puffin scaling from the above

    puffVars.genParams()  # generate rest of scaled params

# generate undulator object

    undmod = undulator(puffVars, undtype = '', Nw = 26, ux = 1., uy = 0.)

# Focusing factor of quads, given by F = (rho B) / (g L), where (rho B)
# is the magnetic rigidity, g is the gradient of the quad field (i.e. dBy/dx)
# and L is the length of the quad. (thin quad approximation)

    qf = 3.14 * puffVars.lg

# Drift length between undulators in meters

    DL = 24. * puffVars.lw # Drift lengths

####    For matching a segment of the beam 'by eye', specify these
####    as the Twiss parameters in the region you want matched
####    In format (geometric emittance, beta, alpha)

    #twxA = [2e-10, 0.002, -0.8]  # (Twiss parameters we have)
    #twyA = [2e-10, 0.002, -0.6]

    #twxB = [7.5e-10, 3.76, -0.189]  # (Twiss parameters to match to)
    #twyB = [7.5e-10, 1.44, 0.0]

####   The routine performs a transfom on the beam contained in fnamein to
####   match the transverse x and y beam phase space to a periodic
####   lattice. The periodic lattice is specified by passing in the
####   undulators, quads, and drifts between undulator modules 
####   (undmod, qf, and DL respectively). Similar to the beam matching
####   in Genesis, it assumes the first element is HALF a quad, focusing
####   in x (defocusing in y). The beam should be in SU units, and contained
####   in SU format in the file with name 'fnamein' below.
####   The Twiss parameters input here (twx1, twy1, twx2, twy2) are optional!
####   The beam is matched from parameters twx1, twy1 to twx2, twy2.
####   The initial Twiss parameters are calculated in the routine from the beam
####   if they are not present. (NOTE that they will then be the Twiss parameters
####   for the FULL beam, incuding the bits you are maybe not bothered about lasing.)
####   The parameters twx2 and twy2 are calculated from the lattice, if not included.
####   The optional inputs are then overrides to what would 'normally' be done.
####   Result is then written to fname2 in SU format.

    fname2 = SU2Matched(fnamein, puffVars, undmod, qf, DL) # , twx1 = twxA, twy1 = twyA, twx2 = twxB, twy2 = twyB)

#   Generate fine beam from sparse beam

    fname3 = CDFupsampler(fname2, puffVars.lr, iloadtype = 0)


#    Convert fine beam to Puffin vars:

    fname4 = SU2Puffin(fname3, puffVars)

    return fname4


if __name__ == '__main__':

    if len(sys.argv)==2:
        fname = sys.argv[1]
        print 'Processing file:', fname
        oname = fullConv(fname)
        print 'Data written to ', oname
    else:
        print 'Usage: SU2Puffin <FileName> \n'
        sys.exit(1)



