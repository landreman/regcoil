#!/usr/bin/env python
print "usage: plot1Doptimization.py regcoil_in.XXX whichM whichN whichOmega"

import sys
from regcoilScan import readVariable
import os
from scipy.io import netcdf
import math
import numpy as np
import matplotlib.pyplot as plt
from nescinUtilities import nescinReadValue

def plot1DOptimization(inputFilename,whichM,whichN,whichOmega):
  
  alpha = readVariable("alpha","float",inputFilename,required=False)
  beta = readVariable("beta","float",inputFilename,required=False)
  gamma = readVariable("gamma","float",inputFilename,required=False)
  if (alpha == None):
    alpha = 0
  if (beta == None):
    beta = 0
  if (gamma == None):
    gamma = 0

  outputFilename = "regcoil_out" + inputFilename[10::] + ".nc"

  directories = filter(os.path.isdir, os.listdir("."))

  chi2_Bs = []
  coil_volumes = []
  coil_plasma_dists = []
  omegas = []
  evals = []
  for directory in directories:

      filename = directory+"/"+outputFilename
      this_input = directory+"/"+inputFilename
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
      evals.append(float(directory[5::]))

      chi2_Bs.append(f.variables["chi2_B"][()][-1])
      coil_volumes.append(f.variables["volume_coil"][()])
      coil_plasma_dists.append(f.variables["coil_plasma_dist_min_lse"][()])
      coil_volume = f.variables["volume_coil"][()]
      nescinFilename = directory+"/"+readVariable("nescin_filename","string",this_input,required=True)
      omega = nescinReadValue(nescinFilename, whichM, whichN, whichOmega)
      omegas.append(omega)

  evals = np.array(evals)
  omegas = np.array(omegas)
  chi2_Bs = np.array(chi2_Bs)
  coil_plasma_dists = np.array(coil_plasma_dists)
  coil_volumes = np.array(coil_volumes)

  objective_functions = chi2_Bs - alpha*coil_plasma_dists - beta*coil_volumes**(1.0/3.0)

  indices = np.argsort(evals)
  evals = evals[indices]
  print evals
  objective_functions = objective_functions[indices]
  omegas = omegas[indices]

  final_omega = omegas[-1]
  final_objective = objective_functions[-1]
  init_omega = omegas[0]
  init_objective = objective_functions[0]

  indices = np.argsort(omegas)
  omegas = omegas[indices]
  objective_functions = objective_functions[indices]

  min_index = np.argmin(objective_functions)
  min_omega = omegas[min_index]
  print "Minimum obtained at: " + str(min_omega)

  plt.figure(facecolor='white')
  plt.axvline(min_omega,label='Minimum',color='blue')
  plt.plot(omegas, objective_functions,marker='.',color='black',linestyle='none')
  plt.axvline(final_omega,label='Final eval',color='green')
  plt.axvline(init_omega,label='Initial eval',color='red')
  plt.legend()
  plt.xlabel(whichOmega + '_' + str(whichM) + '_' + str(whichN))
  plt.ylabel('Objective function')
  plt.show()

if __name__ == "__main__":
  if len(sys.argv) != 5:
    print "Error! You must specify 4 arguments: regcoil_in.XXX whichM whichN whichOmega."
    exit(1)
  whichM = int(sys.argv[2])
  whichN = int(sys.argv[3])
  plot1DOptimization(sys.argv[1],whichM,whichN,sys.argv[4])
