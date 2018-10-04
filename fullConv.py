import sys
import numpy as np
from paraffin.matching import SU2Matched
from paraffin.puffin import puffData
from paraffin.puffin import undulator
from paraffin.upsample import CDFupsampler
from paraffin.SU2Puffin import SU2Puffin

##### FOR NEW OUTER DRIVER:
def fullConv(fnamein):

    puffVars = puffData()

# Initialize CLARA base parameters

    puffVars.aw = 1.41   # The PEAK!!!
    puffVars.gamma0 = 680.0
    puffVars.lw = 0.015
    puffVars.rho = 0.018
    puffVars.undtype = 'curved'
    puffVars.ux = 1.
    puffVars.uy = 0.

#    emitx = 1.022e-9
#    emity = 1.022e-9

# Generate the rest of the Puffin scaling from the above

    puffVars.genParams()  # generate rest of scaled params


# generate undulator

    undmod = undulator(puffVars, undtype = '', Nw = 120, ux = 1., uy = 0.)


# Focusing factor of quads, given by F = (rho B) / (g L), where (rho B)
# is the magnetic rigidity, g is the gradient of the quad field (i.e. dBy/dx)
# and L is the length of the quad. (thin quad approximation)

    qf = 0.0 # 3.14 * puffVars.lg

# Drift length between undulators in meters


    DL = 0. # 24. * puffVars.lw # Drift lengths

####    For matching a segment of the beam 'by eye', specify these
####    as the Twiss parameters in the region you want matched

    twxA = [2e-10, 0.002, -0.8]  # (Twiss parameters we have)
    twyA = [2e-10, 0.002, -0.6]

    twxB = [2e-10, 2.3, 0.0]  # (Twiss parameters to match to)
    twyB = [2e-10, 2.3, 0.0]

    fname2 = SU2Matched(fnamein, puffVars, undmod, qf, DL, twx1 = twxA, twy1 = twyA, twx2 = twxB, twy2 = twyB)



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



