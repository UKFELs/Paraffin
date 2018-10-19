# Copyright 2012-2018, University of Strathclyde
# Authors: Piotr Traczykowski & Lawrence T. Campbell
# License: BSD-3-Clause

####
####
##  'SU' format:
##
## x, y, z - x, y, and z macroparticle coordinates, in metres
## px, py, pz - (px, py, pz)/mc, normalized macroparticle momenta (normalize to 
##              mc, where m is rest electron mass, and c is the speed of light)
## wghts - macroparticle weights, in units if the number of real electrons the
##         the macroparticle represents.


import tables  # for reading and writing hdf5 files

def readSUF(fname):
    f=tables.open_file(fname,'r')
    MPs=f.root.Particles.read()
    f.close()
    return MPs

def writeSUF(fname, MPs):

    output_file=tables.open_file(fname ,'w')

# Particle data:

    ParticleGroup=output_file.create_array('/','Particles', MPs)

# VizSchema metadata    
    
    boundsGroup=output_file.create_group('/','globalGridGlobalLimits','')
    boundsGroup._v_attrs.vsType='limits'
    boundsGroup._v_attrs.vsKind='Cartesian'
    timeGroup=output_file.create_group('/','time','time')
    timeGroup._v_attrs.vsType='time'
    ParticleGroup._v_attrs.vsType='variableWithMesh'
    ParticleGroup._v_attrs.vsTimeGroup='time'
    ParticleGroup._v_attrs.vsNumSpatialDims = 3
    ParticleGroup._v_attrs.vsLimits='globalGridGlobalLimits'
    ParticleGroup._v_attrs.vsLabels='x,px,y,py,z,pz,NE'
    output_file.close()
    