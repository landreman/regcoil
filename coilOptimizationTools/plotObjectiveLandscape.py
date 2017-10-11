#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from scipy.io import netcdf
import math
from regcoilScan import readVariable

def plotObjectiveLandscape(inputFilename):
  
  alpha = readVariable("alpha","float",inputFilename,required=False)
  beta = readVariable("beta","float",inputFilename,required=False)
  gamma = readVariable("gamma","float",inputFilename,required=False)
  if (alpha == None):
    alpha = 0
  if (beta == None):
    beta = 0
  if (gamma == None):
    gamma = 0

  d_min = readVariable("d_min","float",inputFilename,required=False)
  if (d_min is None):
    d_min = 0.2

  outputFilename = "regcoil_out" + inputFilename[10::] + ".nc"

  coilFouriermmin = readVariable("coilFouriermmin","int",inputFilename,required=False)
  coilFouriernmin = readVariable("coilFouriernmin","int",inputFilename,required=False)
  coilFouriermmax = readVariable("coilFouriermmax","int",inputFilename,required=False)
  coilFouriernmax = readVariable("coilFouriernmax","int",inputFilename,required=False)
  sensitivity_option = readVariable("sensitivity_option","int",inputFilename,required=True)

  if (sensitivity_option < 2):
    print "Sensitivity_option must be > 1!"
    sys.exit(0)

  directories = filter(os.path.isdir, os.listdir("."))
  objective_functions = []
  omegas = []
  constraint_equaitons = []
  for whichFourier in ["rmnc","zmns"]:
    for whichM in range(coilFouriermmin, coilFouriermmax+1):
      for whichN in range(coilFouriernmin, coilFouriernmax+1):
        scanVariable = "omega"
        scanVariableName = whichFourier + '.' + str(whichM) + '.' + str(whichN)

        for directory in directories:

          if (directory != "baseCase"):
            thisFourier = directory[0:4]
            thisM,thisN,thisOmega = directory[5::].split(".",2)
            thisM = int(thisM)
            thisN = int(thisN)
            thisOmega = float(thisOmega)

            if (thisM == whichM and thisN == whichN and whichFourier == thisFourier):
              filename = directory+"/"+outputFilename
              currFileName = directory+"/"+inputFilename
              if not os.path.isfile(filename):
                  print "Directory "+directory+" does not have a "+outputFilename+" file (yet)."
                  continue
              try:
                  f = netcdf.netcdf_file(filename,'r',mmap=False)
              except:
                  print "Unable to open "+filename+" even though this file exists."
                  continue
              try:
                  dummy = f.variables["K2"][()]
              except:
                  print "Unable to read "+filename+" even though this file exists."
                  continue
              if math.isnan(dummy[0,0,0]):
                  print "Run in directory "+directory+" has NaNs, so skipping it."
                  continue
              print "Processing directory "+directory

              chi2_B = f.variables["chi2_B"][()][-1]
              coil_volume = f.variables["volume_coil"][()]
              coil_plasma_dist = f.variables["coil_plasma_dist_min_lse"][()]
              objective_function = chi2_B - alpha*coil_plasma_dist - beta*(coil_volume)**(1.0/3.0)
              constraint_equation = d_min - coil_plasma_dist
              constraint_equations.append(constraint_equation)
              objective_functions.append(objective_function)
              omegas.append(thisOmega)

        plt.figure(facecolor='white')
        plt.plot(omegas,objective_functions,'-.r')
        plt.xlabel(scanVariableName)
        plt.ylabel("Objective Function")

        plt.figure(facecolor='white')
        plt.plot(omegas,constraint_equations,'-.r')
        plt.xlabel(scanVariableName)
        plt.ylabel("Constraint Equation")

  plt.show()

if __name__ == "__main__":
  if len(sys.argv) != 2:
    print "Error! You must specify 1 arguments: regcoil_in.XXX."
    exit(1)
  plotObjectiveLandscape(sys.argv[1])



