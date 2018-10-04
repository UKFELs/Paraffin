# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 17:03:14 2015

@author: piotrt
"""
# Import necessary libraries
import tables
import numpy as np
import gc
import sys
from .. import SU  # SU Format
#import getTwiss
#import matchTwiss
from ..puffin import puffData
#from undulator import undulator

def SU2Puffin(fnamein, puffVars = False):

#    Store basename of file

    file_name_base  = (fnamein.split('.')[0]).strip()

# Initialize CLARA base parameters

    if (puffVars == False):
        puffVars = puffData()
        puffVars.aw = 0.8667*np.sqrt(2.)   # The PEAK!!!
        puffVars.gamma0 = 469.0
        puffVars.lw = 0.025
        puffVars.rho = 0.005
        puffVars.undtype = ''
        puffVars.ux = 1.
        puffVars.uy = 0.

        # Generate the rest of the Puffin scaling from the above

        puffVars.genParams()  # generate rest of scaled params



#    Get beam distribution and scale

    MPs = SU.readSUF(fnamein)

    x  = MPs[:,0]
    px = MPs[:,1]
    y  = MPs[:,2]
    py = MPs[:,3]
    z  = MPs[:,4]
    pz = MPs[:,5]
    wghts = MPs[:,6]

    p_tot=np.sqrt((px[:]**2)+(py[:]**2)+(pz[:]**2))


    gamma = (np.sqrt(1+(p_tot)**2))
    gamma = gamma / puffVars.gamma0

    xb = x / (np.sqrt(puffVars.lg * puffVars.lc))
    yb = y / (np.sqrt(puffVars.lg * puffVars.lc))

    pxb = px[:]/(puffVars.aw)
    pyb = -1.0 * py[:]/(puffVars.aw)

    z2 = -z / puffVars.lc   # (ct - z) / lc
    z2 = z2 - min(z2) + 0.001

    wghts = wghts / puffVars.npkbar


    oname = file_name_base+'_Puffin.h5'

    MPs=np.vstack((xb, yb, z2, pxb, pyb, gamma, wghts)).T

    output_file=tables.open_file(oname,'w')


    ParticleGroup=output_file.create_array('/','electrons', MPs)
    ParticleGroup._v_attrs.vsType='variableWithMesh'
    ParticleGroup._v_attrs.vsTimeGroup='time'
    ParticleGroup._v_attrs.vsNumSpatialDims = 3
    ParticleGroup._v_attrs.vsLimits='globalGridGlobalLimits'
    ParticleGroup._v_attrs.vsLabels='xbar,ybar,z2,pxbar,pybar,Gamma,NE'

    boundsGroup=output_file.create_group('/','globalGridGlobalLimits','')
    boundsGroup._v_attrs.vsType='limits'
    boundsGroup._v_attrs.vsKind='Cartesian'

    timeGroup=output_file.create_group('/','time','time')
    timeGroup._v_attrs.vsType='time'

    output_file.close()
    return oname



# Close the file
    #output_file.close()


if __name__ == '__main__':

    if len(sys.argv)==2:
        fname = sys.argv[1]
        print 'Processing file:', fname
        SU2Puffin(fname)
    else:
        print 'Usage: SU2Puffin <FileName> \n'
        sys.exit(1)
        
        

