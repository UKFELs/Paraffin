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
# The output file is saved as HDF5 file with VizSchema metadata added.

#################################################################################
#################################################################################
# THIS SCRIPT CREATES PARTICLES WITH EQUAL WEIGHT - USEFUL WHEN GOING TO SOFTWARE
# WHICH DEMANDS SUCH INPUT E.G. GENESIS OR ELEGANT. DO NOT USE FOR PUFFIN !!!!!!!
#################################################################################
#################################################################################


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
import datetime
import matplotlib.pyplot as plt
from scipy import interpolate

now = datetime.datetime.now()
print 'Current time: ',now.strftime("%Y-%m-%d %H:%M:%S")



if len(sys.argv)==2:
   file_name_in=sys.argv[1]
   print 'Processing file:', file_name_in
else:
   print 'Usage: SU2CDF <Input File Name>\n'
   sys.exit(1)  
file_name_base  = (file_name_in.split('.')[0]).strip()

import PARAMS_CDF

#  Set default parameters and read from parameters file

try:
    k_u=PARAMS_CDF.k_u
except:
    k_u = 228.3636
try:
    a_u=PARAMS_CDF.a_u
except:
    a_u = 1.0121809
try:
    DensityFactor=PARAMS_CDF.DensityFactor
except:
    DensityFactor = 10
try:
    SlicesMultiplyFactor=PARAMS_CDF.SlicesMultiplyFactor
except:
    SlicesMultiplyFactor = 20
try:
    binnumber=PARAMS_CDF.DensitySampling
except:
    binnumber = 15
try:    
    binnumber_X=PARAMS_CDF.X_DensitySampling 
except:
    binnumber_X = 50
try:
    binnumber_Y=PARAMS_CDF.Y_DensitySampling
except:
    binnumber_Y = 50
try:
    binnumber_Z=PARAMS_CDF.Z_DensitySampling
except:
    binnumber_Z = 25
try:
    S_factor=PARAMS_CDF.BeamStretchFactor
except:
    S_factor = 0.0
try:
    NumShapeSlices=PARAMS_CDF.ShapeSampling
except:
    NumShapeSlices = 100


# pre-allocate array and load data into array



f=tables.open_file(file_name_in,'r')
Particles=f.root.Particles.read()

Particles=Particles[Particles[:,4].argsort()]
#Filter the particles with z values
#z_low=-100
#z_high=100

#mA_X=Particles[:,0][(Particles[:,4]>=z_low) & (Particles[:,4]<z_high)]
#mA_Y=Particles[:,2][(Particles[:,4]>=z_low) & (Particles[:,4]<z_high)]
#mA_Z=Particles[:,4][(Particles[:,4]>=z_low) & (Particles[:,4]<z_high)]
#mA_PX=Particles[:,1][(Particles[:,4]>=z_low) & (Particles[:,4]<z_high)]
#mA_PY=Particles[:,3][(Particles[:,4]>=z_low) & (Particles[:,4]<z_high)]
#mA_PZ=Particles[:,5][(Particles[:,4]>=z_low) & (Particles[:,4]<z_high)]
#mA_WGHT=Particles[:,6][(Particles[:,4]>=z_low) & (Particles[:,4]<z_high)]


#*************************************************************
# The below section calculates size of the bin according
# to value of Lc (binnumbers=total_length/Lc)
# USER DATA - MODIFY ACCORDING TO REQUIREMENTS
Pi=np.pi                    # Pi number taken from 'numpy' as more precise than just 3.1415
c=3.0e+8                    # Speed of light
m=9.11e-31                  # mass of electron
e_0=8.854E-12               # vacuum permitivity
e_ch=1.602e-19              # charge of one electron

#*************************************************************


mA_X = Particles[:,0]
mA_Y = Particles[:,2]    
mA_Z = Particles[:,4]

# Descale the momentum units from p/mc to SI - to keep compatibility with rest of calculations
mA_PX = Particles[:,1]*(m*c)
mA_PY = Particles[:,3]*(m*c)  
mA_PZ = Particles[:,5]*(m*c)


mA_WGHT = Particles[:,6]





# The below section calculate some initial data - 4*Pi*Rho is the one most desired
xyz = np.vstack([mA_X,mA_Y,mA_Z]).T
size_x=max(mA_X)-min(mA_X)
size_y=max(mA_Y)-min(mA_Y)
size_z=max(mA_Z)-min(mA_Z)

# This is number of bins just to calculate initial data - don't change if not sure
cube_volume=(size_x*size_y*size_z)/float(binnumber**3)
print 'Volume of sample cube: ',cube_volume
print 'Size of sample X,Y,Z: ', size_x,size_y,size_z
H, edges = np.histogramdd(xyz, bins = (binnumber,binnumber,binnumber),normed=False,weights=mA_WGHT.flat)
n_p=float(np.amax(H))/cube_volume
print 'Particle density Np= ',n_p

# Calculate some needed values: Omega,Rho,Lc,Lambda_U,Lambda_R,
p_tot=np.sqrt((mA_PX[:]**2)+(mA_PY[:]**2)+(mA_PZ[:]**2))
gamma=(np.sqrt(1+(p_tot/(m*c))**2))
gamma_0=np.mean(gamma)
omega_p=np.sqrt((e_ch*e_ch*n_p)/(e_0*m))


rho=(1/gamma_0)*(((a_u*omega_p)/(4*c*k_u))**(2.0/3.0))
lambda_u=(2.0*Pi)/k_u
# Calculated for planar undulator -> (1+(a^2)/2)
lambda_r=(lambda_u/(2.0*gamma_0**2.0))*(1+(a_u**2.0)/2.0)
Lc=lambda_r/(4.0*Pi*rho)
print 'Lc = ',Lc
print 'Rho = ',rho
print 'lambda_u = ',lambda_u
print 'lambda_r = ',lambda_r

# End of inital data calculations


#*************************************************************

# Stack the particles and charge in common array to allow sorting of data along Z axis
xyzW = np.vstack((mA_X.flat,mA_Y.flat,mA_Z.flat,mA_WGHT.flat)).T
# Sort the data along Z axis
xyzW=xyzW[xyzW[:,2].argsort()]


#Calculate total charge
TotalNumberOfElectrons=np.sum(mA_WGHT)

# Assign the sorted data to separate arrays
m_X=xyzW[:,0].flat
m_Y=xyzW[:,1].flat
m_Z=xyzW[:,2].flat
m_WGHT=xyzW[:,3].flat
NumberOfSourceParticles=len(m_X)
InitialParticleCharge=TotalNumberOfElectrons*e_ch

# Print some user useful informations
print 'Initial charge of particles = ',InitialParticleCharge
print 'Total number of electrons = ',TotalNumberOfElectrons
print 'Number of source macroparticles = ',len(m_X)
print 'Electrons/macroparticles in source data = ',round(TotalNumberOfElectrons/len(m_X))
print 'Desired Electrons/macroparticles in output data = ',int(round(TotalNumberOfElectrons/(len(m_X)*DensityFactor)))
print 'Macroparticles in job = ',NumberOfSourceParticles


# Merge XZ and YZ into common array to allow creation of 2D histogram data for XZ and YZ planes          
m_Xm_Z=np.vstack((mA_X.flat,mA_Z.flat)).T
m_Ym_Z=np.vstack((mA_Y.flat,mA_Z.flat)).T


# Create histogram for Z direction and stretch it using S_factor
Hz, edges_Z = np.histogramdd(m_Z, bins = binnumber_Z,range=((min(mA_Z)-S_factor*size_z,max(mA_Z)+S_factor*size_z),(min(mA_Z)-S_factor*size_z,max(mA_Z)+S_factor*size_z)),normed=False,weights=m_WGHT)
edges_YZ = np.histogramdd(m_Ym_Z, bins = (binnumber_Y,binnumber_Z),normed=False,weights=m_WGHT)

# Crate histogram for XZ and YZ planes and stretch it using S_factor variable
HxHz,edges_XZ = np.histogramdd(m_Xm_Z, bins = (binnumber_X,binnumber_Z),range=((min(mA_X)-S_factor*size_x,max(mA_X)+S_factor*size_x),(min(mA_Z)-S_factor*size_z,max(mA_Z)+S_factor*size_z)),normed=False,weights=m_WGHT)
HyHz,edges_YZ = np.histogramdd(m_Ym_Z, bins = (binnumber_Y,binnumber_Z),range=((min(mA_Y)-S_factor*size_y,max(mA_Y)+S_factor*size_y),(min(mA_Z)-S_factor*size_z,max(mA_Z)+S_factor*size_z)),normed=False,weights=m_WGHT)


# Gaussian smoothing of denisty profiles in X/Z and Y/Z
#Hz=ndimage.gaussian_filter(Hz,1.50)
#HxHz=ndimage.gaussian_filter(HxHz,1.50)
#HyHz=ndimage.gaussian_filter(HyHz,1.50)



# Initiate empty array for data points in histogram (move the histogram to array like: Value,X,Y,Z)
XZarr=np.zeros(((len(edges_XZ[0])-1)*(len(edges_XZ[1])-1),3))
YZarr=np.zeros(((len(edges_YZ[0])-1)*(len(edges_YZ[1])-1),3))


#*******Interpolate density*************

# Z axis interpolation 
# Create equispaced points 
x0_Z = np.linspace(0.5*(edges_Z[0][0]+edges_Z[0][1]),0.5*(edges_Z[0][binnumber_Z]+edges_Z[0][binnumber_Z-1]),binnumber_Z)
y0_Z = Hz

#Smoothen the data
from statsmodels.nonparametric.smoothers_lowess import lowess
smoothing_factor=0.01
y0_Z_smth = lowess(y0_Z, x0_Z, frac=smoothing_factor)

# If user want to use LSQ interpolation then unhash next three lines and comment the line where RBF is used 
#z_hst_lngth=np.max(x0_Z)-np.min(x0_Z)
#t_knots_z=np.linspace(np.min(x0_Z)+0.1*z_hst_lngth,np.max(x0_Z)-0.1*z_hst_lngth,5)
#f_Z = interpolate.LSQUnivariateSpline(x0_Z, y0_Z,t_knots_z)

# Use RBF interpolation for Z-axis, hash next lines and unhash 3 lines for LSQ interpolation above  
f_Z = interpolate.Rbf(y0_Z_smth[:,0], y0_Z_smth[:,1])
#f_Z = interpolate.UnivariateSpline(x0_Z, y0_Z_smth[:,1])
#f_Z = interpolate.interp1d(x0_Z, y0_Z,fill_value='extrapolate')
#****Below is just for plotting
m_Z_plt=np.linspace(min(mA_Z)-S_factor*size_z,max(mA_Z)+S_factor*size_z,100)
plt.plot(m_Z_plt,f_Z(m_Z_plt))
plt.show()
#****End of plotting Z axis density profile

# Print the value of the bins with non-zero values in histogram
# This tells us how many histogram values have charge and also informs us how many zero-values we added
# with S_factor variable
Non_Zero_Z=float(np.count_nonzero(Hz))
print 'Non-zero histogram values in Z = ',Non_Zero_Z

# Convert XZ/YZ density histograms to XZ_Density/YZ_Density arrays (move the histogram to array like: Value,X,Y,Z)
for zz in range(1,len(edges_XZ[1])):
  for xx in range(1,len(edges_XZ[0])):
    XZarr[(xx-1)+(zz-1)*(len(edges_XZ[0])-1),0]=(edges_XZ[0][xx]+edges_XZ[0][xx-1])*0.5
    XZarr[(xx-1)+(zz-1)*(len(edges_XZ[0])-1),1]=(edges_XZ[1][zz]+edges_XZ[1][zz-1])*0.5
    XZarr[(xx-1)+(zz-1)*(len(edges_XZ[0])-1),2]=HxHz[xx-1,zz-1]

for zz in range(1,len(edges_YZ[1])):
  for yy in range(1,len(edges_YZ[0])):
    YZarr[(yy-1)+(zz-1)*(len(edges_YZ[0])-1),0]=(edges_YZ[0][yy]+edges_YZ[0][yy-1])*0.5
    YZarr[(yy-1)+(zz-1)*(len(edges_YZ[0])-1),1]=(edges_YZ[1][zz]+edges_YZ[1][zz-1])*0.5
    YZarr[(yy-1)+(zz-1)*(len(edges_YZ[0])-1),2]=HyHz[yy-1,zz-1]

# Calculate the number of slices using the SlicesMultiplyFactor and 4*Pi*rho calculated before
# Do NOT change this line unless you want to set some number of slices not binded with 4*Pi*rho

NumberOfSlices=int(SlicesMultiplyFactor*((max(mA_Z)+S_factor*size_z)-(min(mA_Z)-S_factor*size_z))/(4.0*Pi*rho*Lc))
print 'Number of slices = ',NumberOfSlices

# Calculate the step size - this is just size of samples divided over number of calculated Slices
Step_Size=((max(mA_Z)+S_factor*size_z)-(min(mA_Z)-S_factor*size_z))/NumberOfSlices

#Step_Size=(np.max(m_Z)-np.min(m_Z))/NumberOfSlices

#*** INTERPOLATE XZ AND YZ PLANES USING 2D FUNCTION
# Calculate the length of X,Y,Z histogram for fitting
x_hst_lngth=np.max(XZarr[:,0])-np.min(XZarr[:,0])
y_hst_lngth=np.max(YZarr[:,0])-np.min(YZarr[:,0])
z_hst_lngth=np.max(XZarr[:,1])-np.min(XZarr[:,1])

# Calculate knots (t) needed for LSQBivariateSpline
t_XZ=np.linspace(np.min(XZarr[:,0])+0.1*x_hst_lngth,np.max(XZarr[:,0])-0.1*x_hst_lngth,11)
t_YZ=np.linspace(np.min(YZarr[:,0])+0.1*y_hst_lngth,np.max(YZarr[:,0])-0.1*y_hst_lngth,11)
t_ZZ=np.linspace(np.min(XZarr[:,1])+0.1*z_hst_lngth,np.max(XZarr[:,1])-0.1*z_hst_lngth,11)


# Interpolate using LSQBivariateSpline, hash if want to use Interp2D
f_Dens_XZ=interpolate.LSQBivariateSpline(XZarr[:,0].ravel(), XZarr[:,1].ravel(),XZarr[:,2].ravel(),t_XZ,t_ZZ)
f_Dens_YZ=interpolate.LSQBivariateSpline(YZarr[:,0].ravel(), YZarr[:,1].ravel(),YZarr[:,2].ravel(),t_YZ,t_ZZ)

# Interpolate using Interp2D - unhash if want to use instead of LSQ
#f_Dens_XZ=interpolate.interp2d(XZarr[:,0].ravel(), XZarr[:,1].ravel(),XZarr[:,2].ravel())
#f_Dens_YZ=interpolate.interp2d(YZarr[:,0].ravel(), YZarr[:,1].ravel(),YZarr[:,2].ravel())

#*****Plot the density profile XZ,YZ
# Create data for plot
PLT_X=np.linspace(min(mA_X)-S_factor*size_x,max(mA_X)+S_factor*size_x,100)
PLT_Y=np.linspace(min(mA_Y)-S_factor*size_y,max(mA_Y)+S_factor*size_y,100)
PLT_Z=np.linspace(min(mA_Z)-S_factor*size_z,max(mA_Z)+S_factor*size_z,100)
plt.pcolormesh(PLT_Z, PLT_X,f_Dens_XZ(PLT_X, PLT_Z))
plt.show()
plt.pcolormesh(PLT_Z, PLT_Y,f_Dens_YZ(PLT_Y, PLT_Z))
plt.show()
#*****End of plot

# Initiate data for fitting density profile in each loop
Fit_Particle_Number=1000

New_X=np.linspace(min(mA_X)-S_factor*size_x,max(mA_X)+S_factor*size_x,Fit_Particle_Number)
New_Y=np.linspace(min(mA_Y)-S_factor*size_y,max(mA_Y)+S_factor*size_y,Fit_Particle_Number)
New_Z=np.zeros(Fit_Particle_Number)

# Initiate empty array for Z positions or particles
#density_Z=np.zeros(Num_Of_Slice_Particles)


# Calculate the min/max values for x/y along z-axis (outer shape)

minz=np.min(mA_Z)
maxz=np.max(mA_Z)
step=(maxz-minz)/NumShapeSlices
# Generate step for rescaled lenght (s_factor)
step_s_factor=((maxz+(S_factor*size_z))-(minz-(S_factor*size_z)))/NumShapeSlices
mmax_X=np.zeros(0)
mmin_X=np.zeros(0)
mmax_Y=np.zeros(0)
mmin_Y=np.zeros(0)
mm_Z=np.zeros(0)
# Create interpolated function which describes outer boundaries of initial electron beam
for i in range(0,NumShapeSlices):
#    print minz+step*(i),' ',i
#   Check if there are particles in shape slice
    ShapeParticles = mA_X[(mA_Z>=(minz+step*(i))) & (mA_Z<(minz+step*(i+1)))]
#    print len(ShapeParticles)
    if len(ShapeParticles)>1:
        mmax_X=np.append(mmax_X,np.max(mA_X[(mA_Z>=(minz+step*(i))) & (mA_Z<(minz+step*(i+1)))]))
        mmin_X=np.append(mmin_X,np.min(mA_X[(mA_Z>=(minz+step*(i))) & (mA_Z<(minz+step*(i+1)))]))
        mmax_Y=np.append(mmax_Y,np.max(mA_Y[(mA_Z>=(minz+step*(i))) & (mA_Z<(minz+step*(i+1)))]))
        mmin_Y=np.append(mmin_Y,np.min(mA_Y[(mA_Z>=(minz+step*(i))) & (mA_Z<(minz+step*(i+1)))]))
    #    mm_Z[i]=0.5*((minz+step*(i))+(minz+step*(i+1)))
    #   Project the shape of beam on rescaled length (s_factor) 
        mm_Z=np.append(mm_Z,0.5*(((minz-(S_factor*size_z))+step_s_factor*(i))+((minz-(S_factor*size_z))+step_s_factor*(i+1))))

f_mmax_X=interpolate.UnivariateSpline(mm_Z,mmax_X)
f_mmin_X=interpolate.UnivariateSpline(mm_Z,mmin_X)
f_mmax_Y=interpolate.UnivariateSpline(mm_Z,mmax_Y)
f_mmin_Y=interpolate.UnivariateSpline(mm_Z,mmin_Y)
   
#*** Procedure for placing electrons in each slice according to calculated CDF

# Create values for maximum values of X/Y of particle set
OldRange_X=max(New_X)-min(New_X)
OldRange_Y=max(New_Y)-min(New_Y)
OldMin_X=min(New_X)
OldMax_X=max(New_X)
OldMin_Y=min(New_Y)
OldMax_Y=max(New_Y)

from math import log, floor, ceil, fmod
import numpy as np
def halton(dim, nbpts):
    h = np.empty(nbpts * dim)
    h.fill(np.nan)
    p = np.empty(nbpts)
    p.fill(np.nan)
    P = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    lognbpts = log(nbpts + 1)
    for i in range(dim):
        b = P[i]
        n = int(ceil(lognbpts / log(b)))
        for t in range(n):
            p[t] = pow(b, -(t + 1) )

        for j in range(nbpts):
            d = j + 1
            sum_ = fmod(d, b) * p[0]
            for t in range(1, n):
                d = floor(d / b)
                sum_ += fmod(d, b) * p[t]

            h[j*dim + i] = sum_

    return h.reshape(nbpts, dim)

def SliceCalculate(slice_number):

# Calculate value (in meters) of current slice
    Z_Slice_Value=(slice_number*Step_Size)+(np.min(m_Z)-S_factor*size_z)
    ElectronsPerMacroParticle=((TotalNumberOfElectrons/len(m_X))/DensityFactor/Non_Zero_Z)

    Num_Of_Slice_Particles=abs(int(f_Z(Z_Slice_Value)/NumberOfSlices/ElectronsPerMacroParticle))
    HALTON_XY=halton(2,Num_Of_Slice_Particles)
    density_X=HALTON_XY[:,0]
    density_Y=HALTON_XY[:,1]
    density_Z=np.zeros(Num_Of_Slice_Particles)
    Slice_Ne=np.zeros(Num_Of_Slice_Particles)

    density_Z[:]=Z_Slice_Value

# Scale the range of new particles to new particles set (keep the new particles within original shape of beam)
    NewRange_X=f_mmax_X(Z_Slice_Value)-f_mmin_X(Z_Slice_Value)
    NewRange_Y=f_mmax_Y(Z_Slice_Value)-f_mmin_Y(Z_Slice_Value)
    New_Xl=(((New_X-OldMin_X)*NewRange_X)/OldRange_X)+f_mmin_X(Z_Slice_Value)
    New_Yl=(((New_Y-OldMin_Y)*NewRange_Y)/OldRange_Y)+f_mmin_Y(Z_Slice_Value)

# Interpolate density curev for current slice and remove values below 0
      
    Dens_XZ=f_Dens_XZ(New_X,Z_Slice_Value)
    Dens_YZ=f_Dens_YZ(New_Y,Z_Slice_Value)
    Dens_XZ=Dens_XZ.clip(min=0)
    Dens_YZ=Dens_YZ.clip(min=0)   

# Check if the density is above 0, we don't want to create slices
# with electrons charge = 0
# Note that due to used algorithm macroparticle might have charge less than
# charge of single electron - Puffin accepts this type of data
    if (np.sum(Dens_XZ) > 0 and np.sum(Dens_YZ) > 0): 

# Scale density profile so that the max value is 1
        Dens_XZ=Dens_XZ/np.sum(Dens_XZ)
        Dens_YZ=Dens_YZ/np.sum(Dens_YZ)
   
# Calculate CDF
        cumulative_XZ=np.cumsum(Dens_XZ)
        cumulative_YZ=np.cumsum(Dens_YZ)
    
# Sort the CDF and remove duplicates
        cumulative_nq_XZ=sorted(set(cumulative_XZ))
        cumulative_nq_YZ=sorted(set(cumulative_YZ))  

# Create linear array for interpolation of CDF onto new X,Y,Z
        xx_0_XZ=np.linspace(np.min(New_Xl),np.max(New_Xl),len(cumulative_nq_XZ))
        xx_0_YZ=np.linspace(np.min(New_Yl),np.max(New_Yl),len(cumulative_nq_YZ))

# Create CDF interpolation function using  UnivariateSpline 
        ff_XZ = interpolate.interp1d(cumulative_nq_XZ, xx_0_XZ)
        ff_YZ = interpolate.interp1d(cumulative_nq_YZ, xx_0_YZ)
        
# Calculate the charge for current slice taking into account that there were some
# slices with charge 0 added by using S_factor        
        Slice_Ne[:]=ElectronsPerMacroParticle*Non_Zero_Z

# If the charge of slice is > 0 then apply calculated above CDF and position the electrons
# according to CDF      
        if (np.sum(Slice_Ne))>0:        
            Full_Xl=ff_XZ(density_X)
            Full_Yl=ff_YZ(density_Y)
            Full_Zl=density_Z       
            Full_Nel=Slice_Ne
# Return the values to main program            
        return Full_Xl,Full_Yl,Full_Zl,Full_Nel

#***End of procedure

#***Start multiprocessing of 'SliceCalculate' procedure
import multiprocessing
pool = multiprocessing.Pool()
if __name__ == '__main__':
    result = []

    for slice_number in range(0,NumberOfSlices):
        ProgressValue=100*float(slice_number)/float(NumberOfSlices)
        print 'Completed = ',ProgressValue,' [%]\r',
        pool.apply_async(SliceCalculate, (slice_number,), callback=result.append)
    pool.close()
    pool.join()
#***End of multiprocessing
    
print ' \r'  

#Rearrange data from multiprocessing into more convenient arrays

ParticlesCounter=0
for i in range (0,len(result)):
    for j in range (0,len(list(result[i][0]))):
        ParticlesCounter=ParticlesCounter+1
    ProgressValue=100*float(i)/float(len(result))
    print 'Completed = ',ProgressValue,' [%]\r',
print 'Total number of macroparticles = ',ParticlesCounter


#Inititate empty arrays for positions, momentum and charge of macroparticles 
Full_X=np.zeros(ParticlesCounter)
Full_PX=np.zeros(0)
Full_Y=np.zeros(ParticlesCounter)
Full_PY=np.zeros(0)
Full_Z=np.zeros(ParticlesCounter)
Full_PZ=np.zeros(0)
Full_Ne=np.zeros(ParticlesCounter)

# Loop that rearranges array 'results' into separate X,Y,Z,Ne arrays
print 'Starting rearranging array...'


counter=0
for i in range (0,len(result)):
    SLC_X=list(result[i][0])
    SLC_Y=list(result[i][1])
    SLC_Z=list(result[i][2])
    SLC_Ne=list(result[i][3])
    for j in range (0,len(SLC_X)):
        Full_X[counter]=SLC_X[j]
        Full_Y[counter]=SLC_Y[j]
        Full_Z[counter]=SLC_Z[j]
        Full_Ne[counter]=SLC_Ne[j]
        counter=counter+1
    del SLC_X,SLC_Y,SLC_Z,SLC_Ne
    ProgressValue=100*float(i)/float(len(result))
    print 'Completed = ',ProgressValue,' [%]\r',




# Add noise
print 'Adding noise...'
Rand_Z=(Step_Size*(np.random.random(len(Full_Z)) - 0.50))/np.sqrt(Full_Ne)
Full_Z=Full_Z+Rand_Z

   
# Interpolate momentum onto new macroparticles using the momentum map from initial data   
print 'Starting to interpolate momentum data... - takes time'


def Calculate_PX(): 
    Full_PX = interpolate.griddata((mA_X.ravel(), mA_Y.ravel(), mA_Z.ravel()),mA_PX.ravel(),(Full_X, Full_Y, Full_Z), method='nearest',rescale=True)
    return Full_PX
def Calculate_PY():
    Full_PY = interpolate.griddata((mA_X.ravel(), mA_Y.ravel(), mA_Z.ravel()),mA_PY.ravel(),(Full_X, Full_Y, Full_Z), method='nearest',rescale=True)
    return Full_PY
def Calculate_PZ():
    Full_PZ = interpolate.griddata((mA_X.ravel(), mA_Y.ravel(), mA_Z.ravel()),mA_PZ.ravel(),(Full_X, Full_Y, Full_Z), method='linear',rescale=True)
    return Full_PZ
    
# Start three simultaneous processes for griddata interpolation
from multiprocessing.pool import ThreadPool
pool = ThreadPool(processes=3)
async_result_PX = pool.apply_async(Calculate_PX, ()) 
async_result_PY = pool.apply_async(Calculate_PY, ())
async_result_PZ = pool.apply_async(Calculate_PZ, ())
Full_PX = async_result_PX.get()
Full_PY = async_result_PY.get()
Full_PZ = async_result_PZ.get()
# End of interpolation


# Open output file 
output_file=tables.open_file(file_name_base+'_CDFEQ.h5','w')


# Scale units from SI to p = p/mc

Full_PX=Full_PX/(m*c)
Full_PY=Full_PY/(m*c)
Full_PZ=Full_PZ/(m*c)

# Merge all data into one array
x_px_y_py_z_pz_NE = np.vstack([Full_X.flat,Full_PX.flat,Full_Y.flat,Full_PY.flat,Full_Z.flat,Full_PZ.flat,Full_Ne.flat]).T

# Check for NaN calues in data - happens when using 'linear' option in momentum interpolations (griddata parameter)
NaN_Mask=~np.any(np.isnan(x_px_y_py_z_pz_NE), axis=1)
x_px_y_py_z_pz_NE=x_px_y_py_z_pz_NE[NaN_Mask]

ChargeFactor=InitialParticleCharge/np.sum(x_px_y_py_z_pz_NE[:,6]*e_ch)
print 'Charge scaling factor = ',ChargeFactor
x_px_y_py_z_pz_NE[:,6]=x_px_y_py_z_pz_NE[:,6]*ChargeFactor

print 'Final charge of particles = ',np.sum(x_px_y_py_z_pz_NE[:,6]*e_ch)  
print 'Saving the output to files...'

# Create Group in HDF5 and add metadata for Visit
ParticleGroup=output_file.create_array('/','Particles',x_px_y_py_z_pz_NE)
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
ParticleGroup._v_attrs.FXFELConversionTime=now.strftime("%Y-%m-%d %H:%M:%S")
ParticleGroup._v_attrs.FXFELSourceFileName=file_name_in
ParticleGroup._v_attrs.FXFELDensityFactor=DensityFactor
#Close the file
output_file.close()