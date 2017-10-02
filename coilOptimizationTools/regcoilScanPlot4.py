#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from scipy.io import netcdf
import math
from regcoilScan import readVariable

def regcoilScanPlot4(inputFilename):

  coilFouriermmin = readVariable("coilFouriermmin","int",inputFilename,required=False)
  coilFouriernmin = readVariable("coilFouriernmin","int",inputFilename,required=False)
  coilFouriermmax = readVariable("coilFouriermmax","int",inputFilename,required=False)
  coilFouriernmax = readVariable("coilFouriernmax","int",inputFilename,required=False)
  sensitivity_option = readVariable("sensitivity_option","int",inputFilename,required=True)

  outputFilename = "regcoil_out" + inputFilename[10::] + ".nc"
  plotType = readVariable("plotType","int",inputFilename,required=False)

  if (sensitivity_option < 2):
    print "Sensitivity_option must be > 1!"
    sys.exit(0)

  ilambda = -1
  if (plotType == 1):
    plotVariable = "chi2"
    plotVariableName = "chi2"
    sensitivityVariable = plotVariable
  if (plotType == 2):
    if (sensitivity_option < 3):
      print "Sensitivity_option must be > 2!"
      sys.exit(0)
    plotVariable = "chi2_B"
    plotVariableName = "chi2_B"
    sensitivityVariable = plotVariable
  if (plotType == 3):
    plotVariable = "chi2_K"
    plotVariableName = "chi2_K"
    sensitivityVariable = plotVariable
  if (plotType == 4):
    plotVariable = "volume_coil"
    plotVariableName = "volume_coil"
    sensitivityVariable = "volume_coil"
  if (plotType == 5):
    plotVariable = "LSE_current_density_with_area"
    sensitivityVariable = "LSE_current_density_with_area"
    plotVariableName = plotVariable
  if (plotType == 6):
    plotVariable = "lambda"
    sensitivityVariable = plotVariable
    plotVariableName = plotVariable
  if (plotType == 7):
    plotVariable = "coil_plasma_dist_min"
    sensitivityVariable = plotVariable
    plotVariableName = plotVariable
  if (plotType == 8):
    plotVariable = "coil_plasma_dist_max"
    sensitivityVariable = plotVariable
    plotVariableName = plotVariable


  directories = filter(os.path.isdir, os.listdir("."))

  for whichFourier in ["rmnc","zmns"]:
    for whichM in range(coilFouriermmin, coilFouriermmax+1):
      for whichN in range(coilFouriernmin, coilFouriernmax+1):

        scanVariable = "omega"
        scanVariableName = whichFourier + '.' + str(whichM) + '.' + str(whichN)

        scanVariables = []
        plotVariables = []
        dplotVariabledomegas = []

        for directory in directories:
          # Handle baseCase
          if (directory == "baseCase"):
            continue
            # baseCase = readVariable("omega","float",currFileName,required=False)
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
              lam = f.variables["lambda"][()]
              print "lambda: "+str(lam[ilambda])
              
              scanVariables.append(thisOmega)

              if (plotType == 1):
                dchi2domega = f.variables["dchi2domega"][()][ilambda,:]
                chi2 = f.variables["chi2_B"][()][ilambda] + lam[ilambda]*f.variables["chi2_K"][()][ilambda]
                dplotVariabledomegas.append(dchi2domega)
                plotVariables.append(chi2)
              elif (plotType == 2):
                plotVariables.append(f.variables["chi2_B"][()])
                dplotVariabledomegas.append(f.variables["dchi2Bdomega"][()])
              elif (plotType == 3):
                plotVariables.append(f.variables["chi2_K"][()])
                dplotVariabledomegas.append(f.variables["dchi2Kdomega"][()])
              elif (plotType == 4):
                plotVariables.append(f.variables["volume_coil"][()])
                dvolume_coildomega = f.variables["dvolume_coildomega"][()]
                dplotVariabledomegas.append(dvolume_coildomega)
              elif (plotType == 5):
                plotVariables.append(f.variables["LSE_current_density_with_area"][()])
                dplotVariabledomegas.append(f.variables["dLSE_current_density_with_areadOmega"][()])
              elif (plotType == 6):
                plotVariables.append(f.variables["lambda"][()][ilambda])
                dplotVariabledomegas.append(f.variables["dlambdadomega"][()][ilambda,:])
              elif (plotType == 7):
                plotVariables.append(f.variables["coil_plasma_dist_min_lse"][()])
                dplotVariabledomegas.append(f.variables["dcoil_plasma_dist_mindomega"][()][:])
              elif (plotType == 8):
                plotVariables.append(f.variables["coil_plasma_dist_max_lse"][()])
                dplotVariabledomegas.append(f.variables["dcoil_plasma_dist_maxdomega"][()][:])

              omega = f.variables["omega_coil"][()]
              nomega_coil = f.variables["nomega_coil"][()]
              nmax_sensitivity = f.variables["nmax_sensitivity"][()]
              mmax_sensitivity = f.variables["mmax_sensitivity"][()]
              xn = f.variables["xn_sensitivity"][()]
              xm = f.variables["xm_sensitivity"][()]
              sensitivity_symmetry_option = f.variables["sensitivity_symmetry_option"][()]

        dplotVariabledomegas = np.transpose(np.array(dplotVariabledomegas))

        # Sort by scanVariable
        scanVariables = np.array(scanVariables)
        indices = np.argsort(scanVariables)
        indices = np.array(indices)
        scanVariables = scanVariables[indices]
        plotVariables = np.array(plotVariables)

        if (plotType == 1): # chi2
          plotVariables = plotVariables[indices]
          dplotVariabledscanVariable = dplotVariabledomegas[:,indices]
        if (plotType == 2): # chi2_B
          plotVariables = plotVariables[ilambda,indices]
          dplotVariabledscanVariable = dplotVariabledomegas[:,ilambda,indices]
        if (plotType == 3): # chi2_K
          plotVariables = plotVariables[ilambda,indices]
          dplotVariabledscanVariable = dplotVariabledomegas[:,ilambda,indices]
        if (plotType == 4): # Volume
          plotVariables = plotVariables[indices]
          dplotVariabledscanVariable = dplotVariabledomegas[:,indices]
        if (plotType == 5): # LSE_current_potential_with_area
          plotVariables = plotVariables[ilambda,indices]
          dplotVariabledscanVariable = dplotVariabledomegas[:,ilambda,indices]
        if (plotType == 6): # lambda
          plotVariables = plotVariables[indices]
          dplotVariabledscanVariable = dplotVariabledomegas[:,indices]
        if (plotType == 7): # coil_plasma_dist_min
          plotVariables = plotVariables[indices]
          dplotVariabledscanVariable = dplotVariabledomegas[:,indices]
        if (plotType == 8): # coil_plasma_dist_max
          plotVariables = plotVariables[indices]
          dplotVariabledscanVariable = dplotVariabledomegas[:,indices]

      # Compute omega index corresponding to mode used for geometry scan
        if (whichM > 0):
          index_mode = 2*(nmax_sensitivity+1) + 2*(whichM - 1)*(2*nmax_sensitivity+1) + 2*(nmax_sensitivity + whichN)
        else:
          index_mode = (whichN)*2
        if (whichFourier == "rmnc"):
          offset = 0
        elif (whichFourier == "zmns"):
          offset = 1
        else:
          print "Error!"
          exit
        index_mode = index_mode + offset

        dplotVariabledscanVariable = dplotVariabledscanVariable[index_mode,:]
        if (omega[index_mode] != offset + 1):
          print "Error! Incorrect index"
          exit
        if (xm[index_mode] != whichM):
          print "Error! Incorrect index"
          exit
        if (xn[index_mode] != whichN):
          print "Error! Incorrect index"
          exit

        # Perform finite differencing
        finitediffdplotVariabledscanVariable = np.zeros(scanVariables.shape)
        for i in range(0,len(scanVariables)):
          if (i==0):
            if (scanVariables[i+1] != scanVariables[i]):
              finitediffdplotVariabledscanVariable[i] = (plotVariables[i+1]-plotVariables[i])/(scanVariables[i+1]-scanVariables[i])
          elif (i==len(scanVariables)-1):
            if (scanVariables[i-1] != scanVariables[i]):
              finitediffdplotVariabledscanVariable[i] = (plotVariables[i]-plotVariables[i-1])/(scanVariables[i]-scanVariables[i-1])
          else:
            if (scanVariables[i-1] != scanVariables[i+1]):
              finitediffdplotVariabledscanVariable[i] = (plotVariables[i+1]-plotVariables[i-1])/(scanVariables[i+1]-scanVariables[i-1])

        plt.figure(facecolor='white')
        plt.plot(scanVariables, finitediffdplotVariabledscanVariable,'-.r',label='Finite differencing')
        plt.plot(scanVariables, dplotVariabledscanVariable,'-.b',label='Gradient Computation')
        plt.xlabel(scanVariableName)
        plt.legend()
        plt.ylabel('d' + plotVariableName + 'd' + scanVariableName)

  plt.show()


