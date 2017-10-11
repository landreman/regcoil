#!/usr/bin/env python
print "usage: plotConstraintEquation regcoil_in.XXX whichM whichN whichOmega"

import sys
from regcoilScan import readVariable
import os
from scipy.io import netcdf
import math
import numpy as np
import matplotlib.pyplot as plt
from nescinUtilities import nescinReadValue


def plotConstraintEquation(inputFilename,whichM, whichN, whichOmega):

  outputFilename = "regcoil_out" + inputFilename[10::] + ".nc"
  
  general_option = readVariable("general_option","int",inputFilename,required=True)
  d_min = readVariable("d_min","float",inputFilename,required=False)
  if (d_min is None):
    d_min = 0.2
  directories = filter(os.path.isdir, os.listdir("."))

  constraint_equations = []
  omegas = []

  for directory in directories:

      filename = directory+"/"+outputFilename
      this_input = directory+"/"+inputFilename
      currFileName = directory+"/"+inputFilename
      nescinFilename = directory+"/"+readVariable("nescin_filename","string",this_input,required=True)
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

      dist = f.variables["coil_plasma_dist_min_lse"][()]
      constraint_equation = dist - d_min
      constraint_equations.append(constraint_equation)
      omega = nescinReadValue(nescinFilename, whichM, whichN, whichOmega)
      omegas.append(omega)

  omegas = np.array(omegas)
  constraint_equations = np.array(constraint_equations)

  indices = np.argsort(omegas)
  omegas = omegas[indices]
  constraint_equations = constraint_equations[indices]

  plt.figure(facecolor='white')
  plt.plot(omegas,constraint_equations,linestyle='none',marker='.',color='black')
  plt.axhline(0)
  plt.xlabel(whichOmega + '_' + str(whichM) + '_' + str(whichN))
  plt.ylabel('Constraint equation (dmin)')
  plt.show()

if __name__ == "__main__":
  if len(sys.argv) != 5:
    print "Error! You must specify 4 arguments: regcoil_in.XXX whichM whichN whichOmega."
    exit(1)
  whichM = int(sys.argv[2])
  whichN = int(sys.argv[3])
  plotConstraintEquation(sys.argv[1],whichM,whichN,sys.argv[4])
