# Copyright (c) 2012-2018, University of Strathclyde
# Authors: Piotr Traczykowski and Lawrence T. Campbell
# License: BSD-3-Clause

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

def getDens(sx, sy, sz, binnumber=False):
    """
    Function to calculate the peak number density of the bunch.
    """

    if (not binnumber):
        binnumber = 20

    sx = Particles[:,0]
    sy = Particles[:,2]    
    sz = Particles[:,4]
    wghts = Particles[:,6]

    # The below section calculate some initial data - 4*Pi*Rho is the one mose desired
    xyz = np.vstack([sx, sy, sz]).T
    lenx = max(sx) - min(sx)
    leny = max(sy) - min(sy)
    lenz = max(sz) - min(sz)
    # This is number of bins just to calculate initial data - don't change if not sure

    cube_volume=(lenx*leny*lenz)/float(binnumber**3)
    H, edges = np.histogramdd(xyz, bins = (binnumber,binnumber,binnumber),normed=False,weights=wghts.flat)
    np = float(np.amax(H))/cube_volume

    return np



def upsample2DDens(x, y, wghts, bnx0, bny0, S_factor, plotbasename = False):
    # Merge XZ and YZ into common array to allow creation of 2D histogram data for XZ and YZ planes          
    xy=np.vstack((x.flat, y.flat)).T

    lenx = np.amax(x) - np.amin(x)
    leny = np.amax(y) - np.amin(y)
#HxHz,edges_XZ = np.histogramdd(m_Xm_Z, bins = (binnumber_X,binnumber_Z),range=((min(mA_X)-S_factor*size_x,max(mA_X)+S_factor*size_x),(min(mA_Z)-S_factor*size_z,max(mA_Z)+S_factor*size_z)),normed=False,weights=m_WGHT)
    # Crate histogram for XZ and YZ planes and stretch it using S_factor variable
    HxHy, edges_XY = np.histogramdd(xy, bins = (bnx0, bny0),range=((min(x)-S_factor*lenx,max(x)+S_factor*lenx),(min(y)-S_factor*leny,max(y)+S_factor*leny)),normed=False,weights=wghts)

    # HxHz = ndimage.gaussian_filter(HxHz,1.50)

    if (not plotbasename):
    	qplots = False
    else:
    	qplots = True

    if (qplots):
        XL, YL = np.meshgrid(edges_XY[0],edges_XY[1])

        zyplot = plt.subplot(111)
        zyplot.pcolormesh(XL, YL, HxHy)
        #plt.show()
        plt.savefig(plotbasename + '-ORIG.png')
        plt.clf()




    lexz = len(edges_XY[0])-1
    lexz2 = len(edges_XY[1])-1

    XYarr=np.zeros(( lexz*lexz2, 3))
    

    # Convert XZ/YZ density histograms to XZ_Density/YZ_Density arrays (move the histogram to array like: Value,X,Y,Z)
    for zz in range(1,len(edges_XY[1])):
      for xx in range(1,len(edges_XY[0])):
        XYarr[(xx-1)+(zz-1)*lexz,0]=(edges_XY[0][xx]+edges_XY[0][xx-1])*0.5
        XYarr[(xx-1)+(zz-1)*lexz,1]=(edges_XY[1][zz]+edges_XY[1][zz-1])*0.5
        XYarr[(xx-1)+(zz-1)*lexz,2]=HxHy[xx-1,zz-1]

    #*** INTERPOLATE XZ AND YZ PLANES USING 2D FUNCTION
    # Calculate the length of X,Y,Z histogram for fitting
    x_hst_lngth=np.max(XYarr[:,0])-np.min(XYarr[:,0])
    y_hst_lngth=np.max(XYarr[:,1])-np.min(XYarr[:,1])

    # Calculate knots (t) needed for LSQBivariateSpline
    t_XY=np.linspace(np.min(XYarr[:,0])+0.1*x_hst_lngth,np.max(XYarr[:,0])-0.1*x_hst_lngth,25)
    t_YY=np.linspace(np.min(XYarr[:,1])+0.1*y_hst_lngth,np.max(XYarr[:,1])-0.1*y_hst_lngth,25)


    # Interpolate using LSQBivariateSpline, hash if want to use Interp2D
    f_Dens_XY=interpolate.LSQBivariateSpline(XYarr[:,0].ravel(), XYarr[:,1].ravel(),XYarr[:,2].ravel(),t_XY,t_YY)

    if (qplots):

        PLT_X=np.linspace(min(x)-S_factor*lenx,max(x)+S_factor*lenx,2500)
        PLT_Y=np.linspace(min(y)-S_factor*leny,max(y)+S_factor*leny,2500)

        zxplot = plt.subplot(111)

        zxplot.pcolormesh(PLT_Y, PLT_X,f_Dens_XY(PLT_X, PLT_Y))
        #plt.show()
        plt.savefig(plotbasename + '-UPSAMPLED.png')
        plt.clf()

    return f_Dens_XY



def currentSmoothing(z, wghts, zst, zed, bnz, S_factor, plotname = False):

    from statsmodels.nonparametric.smoothers_lowess import lowess

    if (not plotname):
    	qplots = False
    else:
    	qplots = True

    # Create histogram for Z direction and stretch it using S_factor
    Hz, edges_Z = np.histogramdd(z, bins = bnz,range=((zst,zed),(zst,zed)),normed=False,weights=wghts)

    # Create equispaced points 
    x0_Z = np.linspace(0.5*(edges_Z[0][0]+edges_Z[0][1]),0.5*(edges_Z[0][bnz]+edges_Z[0][bnz-1]),bnz)
    y0_Z = Hz / (edges_Z[0][1] - edges_Z[0][0])

    # If user want to use LSQ interpolation then unhash next three lines and comment the line where RBF is used 
    #z_hst_lngth=np.max(x0_Z)-np.min(x0_Z)
    #t_knots_z=np.linspace(np.min(x0_Z)+0.15*z_hst_lngth,np.max(x0_Z)-0.15*z_hst_lngth,13)
    #f_Z = interpolate.LSQUnivariateSpline(x0_Z, y0_Z,t_knots_z)


    #Smoothen the data
    
    smoothing_factor=0.01
    y0_Z_smth = lowess(y0_Z, x0_Z, frac=smoothing_factor)
    #print y0_Z_smth

    # Use RBF interpolation for Z-axis, hash next lines and unhash 3 lines for LSQ interpolation above  
    #f_Z = interpolate.Rbf(y0_Z_smth[:,0], y0_Z_smth[:,1])
    f_Z = interpolate.UnivariateSpline(x0_Z, y0_Z_smth[:,1])
    

    if (qplots):
    	m_Z_plt=np.linspace(zst,zed,100)
        Iplot = plt.subplot(111)
        Iplot.plot(m_Z_plt,f_Z(m_Z_plt)* 2.99792458e8 * 1.60217653e-19)
        #Iplot.show()
        plt.savefig(plotname)
        #plt.close(Iplot)
        plt.clf()

    return f_Z



def halton(dim, nbpts):
    h = np.empty(nbpts * dim)
    h.fill(np.nan)
    p = np.empty(nbpts)
    p.fill(np.nan)
    P = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    lognbpts = np.log(nbpts + 1)
    for i in range(dim):
        b = P[i]
        n = int(np.ceil(lognbpts / np.log(b)))
        for t in range(n):
            p[t] = pow(b, -(t + 1) )

        for j in range(nbpts):
            d = j + 1
            sum_ = np.fmod(d, b) * p[0]
            for t in range(1, n):
                d = np.floor(d / b)
                sum_ += np.fmod(d, b) * p[t]

            h[j*dim + i] = sum_

    return h.reshape(nbpts, dim)


def Weighted_Stats(values, weights):
#    Return the weighted average and standard deviation

    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2.0, weights=weights)  # Fast and numerically precise
    return (np.sqrt(variance))