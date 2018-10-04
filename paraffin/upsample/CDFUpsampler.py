# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 15:07:49 2016

@author: piotrt
"""
# This scripts uses Cumulative Distribution Function to increase number of particles
# used in FEL simulations.
# The first what has to be done is to sort particles along Z-axis, the purpose of doing so is 
# to easy access the ordered particles in loop 'walking' along the Z-axis. This method saves us
# from reaching low particle number in slice what could happen if we just do the slicing
# along Z-axis using length values (splitting the beam to smaller parts of certain length). Each slice
# has same number of particles.
# The next step is to generate density histogram - this is done by using 'histogramdd' for
# each axis separatyly and also using particle number per macroparticle as weights - I assume
# that sometimes source file will have macroparticle of where electron numbers per particle is not the same
# Then the shape is approximated by linear function and next, proper CDF function is created
# Finally, the CDF function is projected on randomly placed x,y,z values - the Z-axis values are supposed to
# be uniformely spaced that is why 'linspace' is used and not 'random.uniform'.
# The above is done in loop for each slice. The results are appended to larger common array.
# The last step for generating the particle beam is to project momentum onto new set of particles.
# This is done using 'interpolate.griddata' what projects exact shape of momentum. The number of
# electrons per macroparticles is equalized - so the total sum of electrons is same as at the beginnig
# just divided over larger number of macroparticles.
# The output file is saved as HDF5 file with VizSchema metadata added

#############################################################


#################################################################################
#################################################################################
#DISABLE WARNINGS - COMMENT IF YOU WISH TO SEE THE WARNINGS FROM BELOW CODE !!!!!
#################################################################################
#################################################################################
import warnings
warnings.filterwarnings('ignore')
#################################################################################
#################################################################################



import numpy as np
import tables
import sys
import time
from multiprocessing.pool import ThreadPool
from scipy import interpolate
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
import datetime
from ..SU import writeSUF
from ..SU import readSUF


from utilities import getDens
from utilities import upsample2DDens
from utilities import currentSmoothing
from utilities import halton

def elapsed():
        return time.time() - start

def SliceCalculate(slicez, dz, fxz, fyz, fz, xsample, ysample, xg, yg, zsample):

# Calculate value (in meters) of current slice
    #slicez=(slice_number*Step_Size)+(np.min(m_Z)-S_factor*size_z)
    

# Interpolate density curev for current slice and remove values below 0
      
    Dens_XZ = fxz(xg, slicez)  # xsample was new_x
    Dens_YZ = fyz(yg, slicez)
    Dens_XZ = Dens_XZ.clip(min=0)
    Dens_YZ = Dens_YZ.clip(min=0)

    sumx = np.sum(Dens_XZ)
    sumy = np.sum(Dens_YZ)

# Check if the density is above 0, we don't want to create slices
# with electrons charge = 0
# Note that due to used algorithm macroparticle might have charge less than
# charge of single electron - Puffin accepts this type of data
    if (sumx > 0 and sumy > 0): 

# Scale density profile so that the max value is 1
        Dens_XZ = Dens_XZ / sumx
        Dens_YZ = Dens_YZ / sumy
   
# Calculate CDF
        cumulative_XZ = np.cumsum(Dens_XZ)
        cumulative_YZ = np.cumsum(Dens_YZ)
    
# Sort the CDF and remove duplicates
#        cumulative_nq_XZ=sorted(set(cumulative_XZ))
#        cumulative_nq_YZ=sorted(set(cumulative_YZ))  

        cumulative_nq_XZ = cumulative_XZ
        cumulative_nq_YZ = cumulative_YZ

# Create linear array for interpolation of CDF onto new X,Y,Z
        xx_0_XZ = np.linspace(np.amin(xg), np.amax(xg), len(cumulative_nq_XZ))
        xx_0_YZ = np.linspace(np.amin(yg), np.amax(yg), len(cumulative_nq_YZ))

# Create CDF interpolation function using  interp1d 
        ff_XZ = interpolate.interp1d(cumulative_nq_XZ, xx_0_XZ)
        ff_YZ = interpolate.interp1d(cumulative_nq_YZ, xx_0_YZ)
        
# Calculate the charge for current slice taking into account that there were some
# slices with charge 0 added by using S_factor        
        Slice_Ne = np.ones(len(xsample)) * fz(slicez) * dz
        #print 'slice ch', Slice_Ne

# If the charge of slice is > 0 then apply calculated above CDF and position the electrons
# according to CDF      

        if (np.sum(Slice_Ne))>0:        
            Full_Xl=ff_XZ(xsample)
            Full_Yl=ff_YZ(ysample)
            #if (not zsample):
            #    Full_Zl=slicez*np.ones(len(xsample))
            #else:
            Full_Zl=slicez*np.ones(len(xsample)) + dz*(zsample-0.5)
            Full_Nel=Slice_Ne/len(Slice_Ne)
    #print zsample[1]
    return Full_Xl, Full_Yl, Full_Zl, Full_Nel



def main(nslices, startz, dz, fxz, fyz, fz, xsample, ysample, newx, newy, zsample):
    import multiprocessing
    pool = multiprocessing.Pool()
    #if __name__ == '__main__':
    result = []
    for slice_number in range(0,nslices):
        ProgressValue=100.0*float(slice_number)/float(nslices)
        #print 'slice', slice_number
        print 'Completed = ',ProgressValue,' [%]\r',
        slicez = (slice_number*dz) + startz
        #result = SliceCalculate(slicez, dz, fxz, fyz, fz, xsample, ysample, newx, newy, zsample=zsample)
        pool.apply_async(SliceCalculate, (slicez, dz, fxz, fyz, fz, xsample, ysample, newx, newy, zsample), callback=result.append)
    pool.close()
    pool.join()
    return result


def getP(x, y, z, p, xn, yn, zn):

    newp = interpolate.griddata((x.ravel(), y.ravel(), z.ravel()), p.ravel(), (xn, yn, zn), method='linear',rescale=True)

    return newp


def CDFupsampler(file_name_in, charlen, partsPerSlice = 1000, slicesPerLen = 8, bnx = 40, bny = 40, bnz = 40, S_factor = 0.0, iloadtype = 1):

    file_name_base  = (file_name_in.split('.')[0]).strip()

    #*************************************************************
    # The below section calculates size of the bin according
    # to value of Lc (binnumbers=total_length/Lc)
    # USER DATA - MODIFY ACCORDING TO REQUIREMENTS
    Pi=np.pi                    # Pi number taken from 'numpy' as more precise than just 3.1415
    c=3.0e+8                    # Speed of light
    m=9.11e-31                  # mass of electron
    e_0=8.854E-12               # vacuum permitivity
    e_ch=1.602e-19              # charge of one electron

    start = time.time()
    


# Pre-allocate array and load data into array



    MPs = readSUF(file_name_in)

    x  = MPs[:,0]
    px = MPs[:,1]
    y  = MPs[:,2]
    py = MPs[:,3]
    z  = MPs[:,4]
    pz = MPs[:,5]
    wghts = MPs[:,6]




    # Stack the particles and charge in common array to allow sorting of data along Z axis

    xyzW = np.vstack((x.flat,y.flat,z.flat,wghts.flat)).T
    xyzW=xyzW[xyzW[:,2].argsort()]


    #Calculate total charge
    TotalNumberOfElectrons=np.sum(wghts)

    # Assign the sorted data to separate arrays

    xs = xyzW[:,0].flat
    ys = xyzW[:,1].flat
    zs = xyzW[:,2].flat
    wghts_s = xyzW[:,3].flat

    NumberOfSourceParticles=len(xs)
    InitialParticleCharge=TotalNumberOfElectrons*e_ch
    # Print some user useful informations

    print 'Initial charge of particles = ',InitialParticleCharge
    print 'Total number of electrons = ',TotalNumberOfElectrons
    print 'Number of source macroparticles = ',len(xs)
    print 'Electrons/macroparticles in source data = ',round(TotalNumberOfElectrons/len(xs))
    print 'Macroparticles in job = ',NumberOfSourceParticles




    # Set the factor to extend histogram with ZERO values to smoothen the edges. Set to 0 if not needed.
    # The value of 0.15 means that the histogram will grow 30% in each direction (from -1.30*size to +1.13*size)

    lenx = np.amax(xs) - np.amin(xs)
    leny = np.amax(ys) - np.amin(ys)
    lenz = np.amax(zs) - np.amin(zs)


# ******* 1D upsampling of current

    zst = min(z)-S_factor*lenz
    zed = max(z)+S_factor*lenz
    f_Z = currentSmoothing(z, wghts, zst, zed, bnz, S_factor, plotname = 'current.png')


# *******  END 1D upsampling of current


    fdx = upsample2DDens(x, z, wghts, bnx, bnz, S_factor, plotbasename = 'XZ')
    fdy = upsample2DDens(y, z, wghts, bny, bnz, S_factor, plotbasename = 'YZ')

    #Non_Zero_Z=float(np.count_nonzero(Hz))
    #print 'Non-zero histogram values in Z = ',Non_Zero_Z

    # Calculate the number of slices using the slicesPerLen and characteristic length charlen
    
    NumberOfSlices = int(slicesPerLen * ((max(z)+S_factor*lenz)-(min(z)-S_factor*lenz))/(charlen))
    print 'Number of slices = ',NumberOfSlices


    print 'Number of particles in each slice = ',partsPerSlice



    # Calculate the step size - this is just size of samples divided over number of calculated Slices
    Step_Size=((max(z)+S_factor*lenz)-(min(z)-S_factor*lenz))/NumberOfSlices



    # Initiate data for fitting density profile in each loop
    Fit_Particle_Number=1000

    New_X=np.linspace(min(x)-S_factor*lenx,max(x)+S_factor*lenx,Fit_Particle_Number)
    New_Y=np.linspace(min(y)-S_factor*leny,max(y)+S_factor*leny,Fit_Particle_Number)


    # Initiate empty array for Z positions or particles
    density_Z=np.zeros(partsPerSlice)

    #Initiate random placed particles with values 0-1 which will be next projected using CDF onto X,Y positions
    # The 0-1 is CDF value not the X,Y


    if (iloadtype == 1):
        HALTON_XY=halton(2, partsPerSlice)
        #print 'Halton shape is = ',np.shape(HALTON_XY)
        xseq=HALTON_XY[:,0]
        yseq=HALTON_XY[:,1]
        zseq = np.ones(partsPerSlice) * 0.5
    else:
        HALTON_XY=halton(3, partsPerSlice)
        #print 'Halton shape is = ',np.shape(HALTON_XY)
        xseq=HALTON_XY[:,0]
        yseq=HALTON_XY[:,1]
        zseq=HALTON_XY[:,2]
    


    Slice_Ne=np.zeros(partsPerSlice)

    # Choose below between cubic or linear approximation of shape
    #Smoothen the data
    #from statsmodels.nonparametric.smoothers_lowess import lowess
    #smoothing_factor=0.05
    #mmax_X_smth = lowess(mmax_X, mm_Z, frac=smoothing_factor)
    #mmin_X_smth = lowess(mmin_X, mm_Z, frac=smoothing_factor)
    #mmax_Y_smth = lowess(mmax_Y, mm_Z, frac=smoothing_factor)
    #mmin_Y_smth = lowess(mmin_Y, mm_Z, frac=smoothing_factor)
    #==============================================================================
    #f_mmax_X=interpolate.Rbf(mm_Z,mmax_X,function='linear')
    #f_mmin_X=interpolate.Rbf(mm_Z,mmin_X,function='linear')
    #f_mmax_Y=interpolate.Rbf(mm_Z,mmax_Y,function='linear')
    #f_mmin_Y=interpolate.Rbf(mm_Z,mmin_Y,function='linear')
    #==============================================================================

    # Create values for maximum values of X/Y of particle set
    OldRange_X=max(New_X)-min(New_X)
    OldRange_Y=max(New_Y)-min(New_Y)
    OldMin_X=min(New_X)
    OldMax_X=max(New_X)
    OldMin_Y=min(New_Y)
    OldMax_Y=max(New_Y)


#***Start multiprocessing of 'SliceCalculate' procedure
#==============================================================================
    stz_eff = (np.min(zs)-S_factor*lenz)
    result = main(NumberOfSlices, stz_eff, Step_Size, fdx, fdy, f_Z, xseq, yseq, New_X, New_Y, zsample = zseq)
#==============================================================================
#***End of multiprocessing

    #Rearrange data from multiprocessing into more convenient arrays

    Total_Number_Of_Particles=partsPerSlice*len(result)
    print 'Total particles in set = ',Total_Number_Of_Particles

    #Inititate empty arrays for positions, momentum and charge of macroparticles 
    Full_X=np.zeros(Total_Number_Of_Particles)
    Full_PX=np.zeros(0)
    Full_Y=np.zeros(Total_Number_Of_Particles)
    Full_PY=np.zeros(0)
    Full_Z=np.zeros(Total_Number_Of_Particles)
    Full_PZ=np.zeros(0)
    Full_Ne=np.zeros(Total_Number_Of_Particles)

#    print 'y', Full_Y[0:10]
#    print 'z', Full_Z[0:10]
#    print 'Ne', Full_Ne[0:10]

    # Loopt that rearranges array 'results' into separate X,Y,Z,Ne arrays
    print 'Starting rearranging array...'
    counter=0
    for j in range(0,partsPerSlice):
        for i in range(0,len(result)):
            Full_X[counter]=result[i][0][j]
            Full_Y[counter]=result[i][1][j]
            Full_Z[counter]=result[i][2][j]
            Full_Ne[counter]=result[i][3][j]
            counter=counter+1


    # print 'Distance between slices = ',Step_Size
    #print 'Average noise shift = ',np.mean((Rand_Z/np.sqrt(Full_Ne)))
    #print 'Min and Max noise shift = ',np.min((Rand_Z/np.sqrt(Full_Ne))),' ',np.max((Rand_Z/np.sqrt(Full_Ne)))
    #print 'Noise standard deviation = ',np.std((Rand_Z/np.sqrt(Full_Ne)))   
    # Interpolate momentum onto new macroparticles using the momentum map from initial data   
    print 'Starting to interpolate momentum data... - takes time' 

    # Start three simultaneous processes for griddata interpolation
    pool = ThreadPool(processes=3)
    async_result_PX = pool.apply_async(getP, (x, y, z, px, Full_X, Full_Y, Full_Z)) 
    async_result_PY = pool.apply_async(getP, (x, y, z, py, Full_X, Full_Y, Full_Z))
    async_result_PZ = pool.apply_async(getP, (x, y, z, pz, Full_X, Full_Y, Full_Z))
    Full_PX = async_result_PX.get()
    Full_PY = async_result_PY.get()
    Full_PZ = async_result_PZ.get()
    # End of interpolation






    # Merge all data into one array
    MPs2 = np.vstack([Full_X.flat,Full_PX.flat,Full_Y.flat,Full_PY.flat,Full_Z.flat,Full_PZ.flat,Full_Ne.flat]).T






    # Add noise
    print 'Adding noise...'

    # print 'Add Poisson noise to particle weights...'  (AJTC)
    MPs2[:,6] = np.random.poisson(MPs2[:,6])

    # Remove all particles with weights <= zero
    MPs2 = MPs2[MPs2[:,6] > 0]

    # Add noise to particle positions
    MPs2[:,4] = MPs2[:,4] + (Step_Size*(np.random.random(len(MPs2[:,4])) - 0.50))/np.sqrt(MPs2[:,6])



    # Check for NaN calues in data - happens when using 'linear' option in momentum interpolations (griddata parameter)
    NaN_Mask=~np.any(np.isnan(MPs2), axis=1)
    MPs2=MPs2[NaN_Mask]

    print 'len of parts', MPs2.shape

    # Rescale the charge of new particle set (needed due to S_factor usage)
    ChargeFactor = InitialParticleCharge / np.sum(MPs2[:,6]*e_ch)
    print 'Charge scaling factor = ',ChargeFactor
    MPs2[:,6]=MPs2[:,6]*ChargeFactor

    print 'Final charge of particles = ',np.sum(MPs2[:,6]*e_ch)
    print 'Saving the output to files...'


    # Write to output file
    oname = file_name_base+'_CDF.h5'
    writeSUF(oname, MPs2)
    return oname


if __name__ == '__main__':

    if len(sys.argv) == 3:
        
        fname = sys.argv[1]
        charlen = np.float(sys.argv[2])
        print 'Processing file:', fname
        print 'char length is:', charlen
    
        CDFupsampler(fname, charlen)

    elif (len(sys.argv) == 4):
        fname = sys.argv[1]
        charlen = np.float(sys.argv[2])
        ilt = np.int(sys.argv[3])
        print 'Processing file:', fname
        print 'char length is:', charlen
        print 'Loading type is:', ilt
    
        CDFupsampler(fname, charlen, iloadtype = ilt)


    elif (len(sys.argv) == 8):
        
        fname = sys.argv[1]
        charlen = sys.argv[2]
        partsPerSlice = sys.argv[3]
        slicesPerLen = sys.argv[4]
        bnx = sys.argv[5]
        bny = sys.argv[6]
        bnz = sys.argv[7]
        S_factor = sys.argv[8]
        ilt = sys.argv[9]

        print 'Processing file:', fname

        CDFupsampler(fname, charlen, partsPerSlice = partsPerSlice, slicesPerLen = slicesPerLen, bnx = bnx, bny = bny, bnz = bnz, S_factor = S_factor, iloadtype = ilt)

    else:
        print 'Usage: SU2CDF <FileName> <characteristic length>\n'
        sys.exit(1)





