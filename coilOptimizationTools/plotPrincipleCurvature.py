#!/usr/bin/env python
print "usage: plotPrincipleCurvature regcoil_out.XXX"

import sys
from regcoilScan import readVariable
import os
from scipy.io import netcdf
import math
import numpy as np
import matplotlib.pyplot as plt

def plotPrincipleCurvature(outputFilename):
  
      numContours = 30

      if not os.path.isfile(outputFilename):
          print "The name "+outputFilename+" is not a valid filename."
          sys.exit(0)
      try:
          f = netcdf.netcdf_file(outputFilename,'r',mmap=False)
      except:
          print "Unable to open "+outputFilename+" even though this file exists."
          sys.exit(0)
      try:
          dummy = f.variables["K2"][()]
      except:
          print "Unable to read "+outputFilename+" even though this file exists."
          sys.exit(0)
      if math.isnan(dummy[0,0,0]):
          print "Run has NaNs, so skipping it."
          sys.exit(0)

      compute_curvature = f.variables["compute_curvature"][()]
      if (compute_curvature != 1):
        print "REGCOIL has not been run with compute_curvature == 1!"
        sys.exit(0)

      principle_curvature_1 = f.variables["principle_curvature_1"][()]
      principle_curvature_2 = f.variables["principle_curvature_2"][()]
      gaussian_curvature = f.variables["gaussian_curvature"][()]
      norm_normal = f.variables["norm_normal_plasma"][()]
      mean_curvature = f.variables["mean_curvature"][()]
      theta_coil = f.variables["theta_coil"][()]
      zeta_coil = f.variables["zeta_coil"][()]
      principle_curvature_1 = np.transpose(principle_curvature_1)
      principle_curvature_2 = np.transpose(principle_curvature_2)
      norm_normal = np.transpose(norm_normal)
      theta_coil = np.array(theta_coil)
      zeta_coil = np.array(zeta_coil)
      principle_curvature_1 = np.array(principle_curvature_1)
      principle_curvature_2 = np.array(principle_curvature_2)

      print "max(mean_curvature) : " + str(np.max(mean_curvature))
      print "min(mean_curvature) : " + str(np.min(mean_curvature))

      plt.figure(facecolor='white')
      contours = np.linspace(np.min(norm_normal),np.max(norm_normal),50)
      plt.contourf(zeta_coil,theta_coil,norm_normal,levels=contours,edgecolor=None)
      plt.colorbar()
      plt.title('Norm Normal',fontsize=20)
      plt.xlabel(r'$\zeta$',fontsize=20)
      plt.ylabel(r'$\theta$',fontsize=20)

      plt.figure(facecolor='white')
      contours = np.linspace(np.min(principle_curvature_1),np.max(principle_curvature_1),50)
      plt.contourf(zeta_coil,theta_coil,principle_curvature_1,levels=contours,edgecolor=None)
      plt.colorbar()
      plt.title('Principle Curvature 1',fontsize=20)
      plt.xlabel(r'$\zeta$',fontsize=20)
      plt.ylabel(r'$\theta$',fontsize=20)

      plt.figure(facecolor='white')
      contours = np.linspace(np.min(principle_curvature_2),np.max(principle_curvature_2),50)
      plt.contourf(zeta_coil,theta_coil,principle_curvature_2,levels=contours,edgecolor=None)
      plt.colorbar()
      plt.title('Principle Curvature 2',fontsize=20)
      plt.xlabel(r'$\zeta$',fontsize=20)
      plt.ylabel(r'$\theta$',fontsize=20)

      plt.figure(facecolor='white')
      contours = np.linspace(np.min(mean_curvature),np.max(mean_curvature),50)
      plt.contourf(zeta_coil,theta_coil,mean_curvature,levels=contours,edgecolor=None)
      plt.colorbar()
      plt.title('Mean Curvature',fontsize=20)
      plt.xlabel(r'$\zeta$',fontsize=20)
      plt.ylabel(r'$\theta$',fontsize=20)

      plt.figure(facecolor='white')
      contours = np.linspace(np.min(gaussian_curvature),np.max(gaussian_curvature),50)
      plt.contourf(zeta_coil,theta_coil,gaussian_curvature,levels=contours,edgecolor=None)
      plt.colorbar()
      plt.title('Gaussian Curvature',fontsize=20)
      plt.xlabel(r'$\zeta$',fontsize=20)
      plt.ylabel(r'$\theta$',fontsize=20)

      plt.show()

if __name__ == "__main__":
  if len(sys.argv) != 2:
    print "Error! You must specify 1 arguments: regcoil_out.XXX."
    exit(1)
  plotPrincipleCurvature(sys.argv[1])
