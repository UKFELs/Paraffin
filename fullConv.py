# Copyright 2012-2018, University of Strathclyde
# Authors: Lawrence T. Campbell and Piotr Traczykowski
# License: BSD-3-Clause

"""
Simple example script for using different packages in Paraffin.
"""

import sys
import numpy as np
from paraffin import SU2Matched
from paraffin import puffData
from paraffin import undulator
from paraffin import CDFupsampler
from paraffin import SU2Puffin

##### FOR NEW OUTER DRIVER:
def fullConv(fnamein):

# Create an instance of the Puffin class

    puffVars = puffData()

# Initialize CLARA base parameters

    puffVars.aw = 0.8667*np.sqrt(2.)   # The PEAK!!!
    puffVars.gamma0 = 469.0
    puffVars.lw = 0.025
    puffVars.rho = 0.005
    puffVars.undtype = ''
    puffVars.ux = 1.
    puffVars.uy = 0.

#    emitx = 1.022e-9
#    emity = 1.022e-9

# Generate the rest of the Puffin scaling from the above

    puffVars.genParams()  # generate rest of scaled params


# generate undulator

    undmod = undulator(puffVars, undtype = '', Nw = 26, ux = 1., uy = 0.)


# Focusing factor of quads, given by F = (rho B) / (g L), where (rho B)
# is the magnetic rigidity, g is the gradient of the quad field (i.e. dBy/dx)
# and L is the length of the quad. (thin quad approximation)

    qf = 3.14 * puffVars.lg

# Drift length between undulators in meters


    DL = 24. * puffVars.lw # Drift lengths

####    For matching a segment of the beam 'by eye', specify these
####    as the Twiss parameters in the region you want matched
####
####    Twiss parameters are in format: [emittance, beta, alpha] (geometric emmitance)

    twxA = [2e-10, 0.002, -0.8]  # (Twiss parameters we have)
    twyA = [2e-10, 0.002, -0.6]
    
    twxB = [2e-10, 2.3, 0.0]  # (Twiss parameters to match to)
    twyB = [2e-10, 2.3, 0.0]

## The function SU2Matched transforms the distribution in fnamein (in SU format)
## FROM twx1, twy1 to twx2, twy2. 
##
## twx1 are the Twiss parameters in x, and are calculated in the function from 
## the beam. However, this calculation takes into account the WHOLE beam, but
## you may actually only wish to match a section of the beam. twx1 and twy1 
## can be overridden as inputs to the function
##
## Similarly, twx2 are the matched Twiss parameters, matched to the FODO lattice
## specified in the input, but can be overridden in the input to the function.
## If the qf quad factor is = 0, then the Twiss parameters to match to are
## calculated from the undulator focusing ONLY (no FODO).

## Here, there is no twx1 or twy1 specified, so it will be calculated from the
## beam file, global across the whole beam.

    fname2 = SU2Matched(fnamein, puffVars, undmod, qf, DL, twx2 = twxB, twy2 = twyB)

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



