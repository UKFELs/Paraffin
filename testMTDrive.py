import sys
import numpy as np
from paraffin.matching import SU2Matched as TMT
from paraffin.puffin import puffData
from paraffin.puffin import undulator

# import SUMatch 


##### FOR NEW OUTER DRIVER:
def MTdriver(fnamein):

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

    undmod = undulator(puffVars, undtype = '', Nw = 26, ux = 1., uy = 0.)


# Focusing factor of quads, given by F = (rho B) / (g L), where (rho B)
# is the magnetic rigidity, g is the gradient of the quad field (i.e. dBy/dx)
# and L is the length of the quad. (thin quad approximation)

    qf = 3.14 * puffVars.lg

# Drift length between undulators in meters


    DL = 24. * puffVars.lw # Drift lengths

    #fnamein = 'short_CLARA_001_A2SU.h5'

#    TMT.SU2Matched(fnamein, puffVars, undmod, qf, DL)

####    For matching a segment of the beam 'by eye', specify these
####    as the Twis parameters in the region you want matched

    twx = [7.5e-10, 3.76, -0.189]  # (Twiss parameters to match to)
    twy = [7.5e-10, 1.44, 0.]

    fout = TMT(fnamein, puffVars, undmod, qf, DL, twx2 = twx, twy2 = twy)

if __name__ == '__main__':

    if len(sys.argv)==2:
        fname = sys.argv[1]
        print 'Processing file:', fname
        MTdriver(fname)
    else:
        print 'Usage: SU2Puffin <FileName> \n'
        sys.exit(1)
