#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from scipy.io import netcdf
import math
from regcoilScan import readVariable
from regcoilScanPlot4 import regcoilScanPlot4

def regcoilScanPlot(inputFilename):

  inputFilename = sys.argv[1]
  outputFilename = "regcoil_out" + inputFilename[10::] + ".nc"
  scanType = readVariable("scanType","int",inputFilename,required=False)
  plotType = readVariable("plotType","int",inputFilename,required=False)
  if (scanType == 4):
    regcoilScanPlot4(inputFilename)



  ilambda = -1
  if (plotType == 1):
    plotVariable = "chi2"
    plotVariableName = "chi2"
    sensitivityVariable = plotVariable
  if (plotType == 6):
    plotVariable = "chi2_B"
    plotVariableName = "chi2_B"
    sensitivityVariable = plotVariable
  if (plotType == 7):
    plotVariable = "chi2_K"
    plotVariableName = "chi2_K"
    sensitivityVariable = plotVariable
  if (plotType == 16):
    plotVariable = "volume_coil"
    plotVariableName = "volume_coil"
    sensitivityVariable = "volume_coil"
  if (plotType == 18):
    plotVariable = "LSE_current_density_with_area"
    sensitivityVariable = "LSE_current_density_with_area"
    plotVariableName = plotVariable
  if (plotType == 20):
    plotVariable = "area_coil"
    sensitivityVariable = plotVariable
    plotVariableName = plotVariable
  if (plotType == 21):
    plotVariable = "lambda"
    sensitivityVariable = plotVariable
    plotVariableName = plotVariable
  if (plotType == 22):
    plotVariable = "coil_plasma_dist"
    sensitivityVariable = plotVariable
    plotVariableName = plotVariable

  plotVariables = []
  dplotVariabledomegas = []
  directories = filter(os.path.isdir, os.listdir("."))

  if (scanType == 1):
    scanVariable = "a_coil"
    scanVariableName = scanVariable
  if (scanType == 2):
    scanVariable = "R0_coil"
    scanVariableName = scanVariable
  if (scanType == 3):
    scanVariable = "R0_plasma"
    scanVariableName = scanVariable
  if (scanType == 4):
    scanVariable = "omega"
    scanVariableName = coilWhichFourier + '.' + str(coilFourierm) + '.' + str(coilFouriern)

  scanVariables = []

  for directory in directories:
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
      if (scanType < 4):
        scanVariables.append(f.variables[scanVariable][()])
        if (directory == "baseCase"): baseCase = f.variables[scanVariable][()]
      else:
        scanVariables.append(readVariable("omega","float",currFileName,required=False))
        if (directory == "baseCase"): baseCase = readVariable("omega","float",currFileName,required=False)
      if (plotType == 0):
        plotVariables.append(f.variables["chi2_B"][()] + lam[ilambda]*f.variables["chi2_K"][()])

      elif (plotType == 1):
        dchi2domega = f.variables["dchi2domega"][()][ilambda,:]
        chi2 = f.variables["chi2_B"][()][ilambda] + lam[ilambda]*f.variables["chi2_K"][()][ilambda]
        dplotVariabledomegas.append(dchi2domega)
        plotVariables.append(chi2)
      elif (plotType == 6):
        plotVariables.append(f.variables["chi2_B"][()])
        dplotVariabledomegas.append(f.variables["dchi2Bdomega"][()])
      elif (plotType == 7):
        plotVariables.append(f.variables["chi2_K"][()])
        dplotVariabledomegas.append(f.variables["dchi2Kdomega"][()])
      elif (plotType == 5):
        if (component == 0):
          dplotVariabledomegas.append(f.variables["dnormxdomega"][()])
          plotVariables.append(f.variables["normal_coil"][()][:,:,0])
        elif (component == 1):
          dplotVariabledomegas.append(f.variables["dnormydomega"][()])
          plotVariables.append(f.variables["normal_coil"][()][:,:,1])
        elif (component == 2):
          dplotVariabledomegas.append(f.variables["dnormzdomega"][()])
          plotVariables.append(f.variables["normal_coil"][()][:,:,2])
      elif (plotType == 10):
        dplotVariabledomegas.append(f.variables["dddomega"][()][:,:,component])
        if (component == 0):
          plotVariables.append(f.variables["d_x"][()])
        elif (component == 1):
          plotVariables.append(f.variables["d_y"][()])
        elif (component == 2):
          plotVariables.append(f.variables["d_z"][()])
      elif (plotType == 11):
        plotVariables.append(f.variables["chi2_B"][()] + lam[ilambda]*f.variables["chi2_K"][()])
        dchi2Bdomega = f.variables["dchi2Bdomega"][()]
        dchi2Kdomega = f.variables["dchi2Kdomega"][()]
        dchi2domega = np.zeros(dchi2Bdomega.shape)
        for ilambda in range(0,len(lam)):
          dchi2domega[ilambda,:] = dchi2Bdomega[ilambda,:] + lam[ilambda]*dchi2Kdomega[ilambda,:]
        dplotVariabledomegas.append(dchi2domega)
      elif (plotType == 16):
        plotVariables.append(f.variables["volume_coil"][()])
        dvolume_coildomega = f.variables["dvolume_coildomega"][()]
        dplotVariabledomegas.append(dvolume_coildomega)
      elif (plotType == 17):
        plotVariables.append(f.variables["single_valued_current_potential_mn"][()])
        dSolutiondOmega = f.variables["dSolutiondOmega"][()]
        dplotVariabledomegas.append(dSolutiondOmega)
      elif (plotType == 18):
        plotVariables.append(f.variables["LSE_current_density_with_area"][()])
        dplotVariabledomegas.append(f.variables["dLSE_current_density_with_areadOmega"][()])
      elif (plotType == 19):
        chi2_K = f.variables["chi2_K"][()]
        area_coil = f.variables["area_coil"][()]
        plotVariables.append(np.sqrt(chi2_K/area_coil))
        dplotVariabledomegas.append(f.variables["dRMSKdomega"][()])
      elif (plotType == 20):
        plotVariables.append(f.variables["area_coil"][()])
        dplotVariabledomegas.append(f.variables["darea_coildomega"][()])
      elif (plotType == 21):
        plotVariables.append(f.variables["lambda"][()][ilambda])
        print "f.variables[dlambdadomega][()].shape: " + str(f.variables["dlambdadomega"][()].shape)
        dplotVariabledomegas.append(f.variables["dlambdadomega"][()][ilambda,:])
      elif (plotType == 22):
        plotVariables.append(f.variables["coil_plasma_dist"][()])
        dplotVariabledomegas.append(f.variables["dcoil_plasma_distdomega"][()][:])
      else:
        plotVariables.append(f.variables[plotVariable][()])
        dplotVariabledomegas.append(f.variables["d" + sensitivityVariable + "domega"][()])


  if (plotType != 0):
    omega = f.variables["omega_coil"][()]
    nomega_coil = f.variables["nomega_coil"][()]
    nmax_sensitivity = f.variables["nmax_sensitivity"][()]
    mmax_sensitivity = f.variables["mmax_sensitivity"][()]
    xn = f.variables["xn_sensitivity"][()]
    xm = f.variables["xm_sensitivity"][()]
    sensitivity_symmetry_option = f.variables["sensitivity_symmetry_option"][()]
    
  dplotVariabledomegas = np.transpose(np.array(dplotVariabledomegas))

  ntheta_plasma = f.variables["ntheta_plasma"][()]
  ntheta_coil = f.variables["ntheta_coil"][()]
  nzeta_coil = f.variables["nzeta_coil"][()]
  scanVariables = np.transpose(scanVariables)
  plotVariables = np.transpose(np.array(plotVariables))

  if (plotType == 2 or plotType == 3 or plotType == 8):
    index_plasma = (izeta_plasma)*ntheta_plasma + (itheta_plasma+1)
  if (plotType == 3 or plotType == 9 or plotType == 10):
    index_coil = (izeta_coil)*ntheta_coil + (itheta_coil+1)

  print "dplotVariabledomegas.shape: " + str(dplotVariabledomegas.shape)
  print "plotVariables.shape: " + str(plotVariables.shape)

  # Sort by scanVariable
  indices = np.argsort(scanVariables)
  print "indices.shape: " + str(indices.shape)
  scanVariables = scanVariables[indices]

  # Take elements of plotVariable and dplotVariabledscanVariable arrays for plotting
  if (plotType == 0):
    plotVariables = plotVariables[ilambda,indices]
  if (plotType == 1): # chi2
    plotVariables = plotVariables[indices]
    dplotVariabledscanVariable = dplotVariabledomegas[:,indices]
  if (plotType == 2): # g
    plotVariables = plotVariables[index_plasma,i_basis_function,indices]
    dplotVariabledscanVariable = dplotVariabledomegas[:,index_plasma,i_basis_function,indices]
  if (plotType == 3): # inductance
    plotVariables = plotVariables[index_plasma,index_coil,indices]
    dplotVariabledscanVariable = dplotVariabledomegas[:,index_plasma,index_coil,indices]
  if (plotType == 4): # norm_normal_coil
    plotVariables = plotVariables[itheta_coil,izeta_coil,indices]
    dplotVariabledscanVariable = dplotVariabledomegas[:,itheta_coil,izeta_coil,indices]
  if (plotType == 5): # normal_coil
    plotVariables = plotVariables[itheta_coil,izeta_coil,indices]
    dplotVariabledscanVariable = dplotVariabledomegas[:,itheta_coil,izeta_coil,indices]
  if (plotType == 6): # chi2_B
    plotVariables = plotVariables[ilambda,indices]
    dplotVariabledscanVariable = dplotVariabledomegas[:,ilambda,indices]
  if (plotType == 7): # chi2_K
    plotVariables = plotVariables[ilambda,indices]
    dplotVariabledscanVariable = dplotVariabledomegas[:,ilambda,indices]
  if (plotType == 8): # h
    plotVariables = plotVariables[index_plasma,indices]
    dplotVariabledscanVariable = dplotVariabledomegas[:,index_plasma,indices]
  if (plotType == 9): # f
    plotVariables = plotVariables[index_coil,i_basis_function,indices]
    dplotVariabledscanVariable = dplotVariabledomegas[:,index_coil,i_basis_function,indices]
  if (plotType == 10): # d
    plotVariables = plotVariables[index_coil,indices]
    dplotVariabledscanVariable = dplotVariabledomegas[:,index_coil]
  if (plotType == 11): # chi2 (added with adjoint)
    plotVariables = plotVariables[ilambda,indices]
    dplotVariabledscanVariable = dplotVariabledomegas[:,ilambda,indices]
  if (plotType == 12): # RHS_K
    plotVariables = plotVariables[i_basis_function,indices]
    dplotVariabledscanVariable = dplotVariabledomegas[:,i_basis_function,indices]
  if (plotType == 13): # RHS_B
    plotVariables = plotVariables[i_basis_function,indices]
    dplotVariabledscanVariable = dplotVariabledomegas[:,i_basis_function,indices]
  if (plotType == 14): # matrix_K
    plotVariables = plotVariables[i_basis_function,j_basis_function,indices]
    dplotVariabledscanVariable = dplotVariabledomegas[:,i_basis_function,j_basis_function,indices]
  if (plotType == 15): # matrix_B
    plotVariables = plotVariables[i_basis_function,j_basis_function,indices]
    dplotVariabledscanVariable = dplotVariabledomegas[:,i_basis_function,j_basis_function,indices]
  if (plotType == 16): # Volume
    plotVariables = plotVariables[indices]
    dplotVariabledscanVariable = dplotVariabledomegas[:,indices]
  if (plotType == 17): # solution
    plotVariables = plotVariables[i_basis_function, ilambda, indices]
    dplotVariabledscanVariable = np.transpose(dplotVariabledomegas[ilambda,i_basis_function,:,indices])
  if (plotType == 18): # LSE_current_potential_with_area
    plotVariables = plotVariables[ilambda,indices]
    dplotVariabledscanVariable = dplotVariabledomegas[:,ilambda,indices]
  if (plotType == 19): # RMSK
    plotVariables = plotVariables[ilambda,indices]
    dplotVariabledscanVariable = dplotVariabledomegas[:,ilambda,indices]
  if (plotType == 20): # area_coil
    plotVariables = plotVariables[indices]
    dplotVariabledscanVariable = dplotVariabledomegas[:,indices]
  if (plotType == 21): # lambda
    plotVariables = plotVariables[indices]
    dplotVariabledscanVariable = dplotVariabledomegas[:,indices]
  if (plotType == 22): # coil_plasma_dist
    plotVariables = plotVariables[indices]
    dplotVariabledscanVariable = dplotVariabledomegas[:,indices]
  print dplotVariabledscanVariable.shape
  print plotVariables.shape

  # Compute omega index corresponding to mode used for geometry scan
  if (scanType == 1):
    if (sensitivity_symmetry_option == 3):
      index_mode = 4*((nmax_sensitivity+1) + nmax_sensitivity)
    else:
      index_mode = 2*((nmax_sensitivity+1) + nmax_sensitivity)
    dplotVariabledscanVariable = dplotVariabledscanVariable[index_mode,:] \
      + dplotVariabledscanVariable[index_mode+1,:]
  if (scanType == 2):
    index_mode = 0
    if (plotType != 0):
      dplotVariabledscanVariable = dplotVariabledscanVariable[index_mode,:]
  if (scanType == 4):
    if (sensitivity_symmetry_option == 3):
      if (coilFourierm > 0):
        index_mode = 4*(nmax_sensitivity+1) + 4*(coilFourierm - 1)*(2*nmax_sensitivity+1) + 4*(nmax_sensitivity + coilFouriern)
      else:
        index_mode = (coilFouriern)*4
    else:
      if (coilFourierm > 0):
        index_mode = 2*(nmax_sensitivity+1) + 2*(coilFourierm - 1)*(2*nmax_sensitivity+1) + 2*(nmax_sensitivity + coilFouriern)
      else:
        index_mode = (coilFouriern)*2
    if (coilWhichFourier == "rmnc"):
      offset = 0
    elif (coilWhichFourier == "zmns"):
      offset = 1
    elif (coilWhichFourier == "rmns"):
      offset = 2
    elif (coilWhichFourier == "zmnc"):
      offset = 3
    else:
      print "Error!"
      exit
    index_mode = index_mode + offset
    print("omega: %d" % (omega[index_mode]))
    print("m: %d" % (xm[index_mode]))
    print("n: %d" % (xn[index_mode]))

    dplotVariabledscanVariable = dplotVariabledscanVariable[index_mode,:]
    if (omega[index_mode] != offset + 1):
      print "Error! Incorrect index"
      exit
    if (xm[index_mode] != coilFourierm):
      print "Error! Incorrect index"
      exit
    if (xn[index_mode] != coilFouriern):
      print "Error! Incorrect index"
      exit

  plotVariables = np.array(plotVariables)
  scanVariables = np.array(scanVariables)
  print "scanVariables.shape: " + str(scanVariables.shape)
  print "plotVariables.shape: " + str(plotVariables.shape)

  plt.figure(facecolor='white')
  plt.plot(scanVariables,plotVariables, linestyle="none",marker="o",label=plotVariableName)
  plt.legend()
  plt.axvline(baseCase,color='black',linestyle='dashed')
  plt.ylabel(plotVariableName)
  plt.xlabel(scanVariableName)

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

  print "finitediffdplotVariabledscanVariable.shape: " + str(finitediffdplotVariabledscanVariable.shape)
  print "dplotVariabledscanVariable.shape: " + str(dplotVariabledscanVariable.shape)

  if (plotType != 0):
    residual = np.abs((dplotVariabledscanVariable-finitediffdplotVariabledscanVariable)/dplotVariabledscanVariable)
    
    plt.figure(facecolor='white')
    plt.title('REGCOIL Output ')
    plt.plot(scanVariables, dplotVariabledscanVariable,linestyle='none',marker='o')
    plt.xlabel(scanVariableName)
    plt.axvline(baseCase,color='black',linestyle='dashed')
    plt.ylabel('d' + plotVariableName + 'd' + scanVariableName)

    plt.figure(facecolor='white')
    plt.plot(scanVariables, finitediffdplotVariabledscanVariable,linestyle='none',marker='o',label='finite diff')
    plt.plot(scanVariables, dplotVariabledscanVariable,linestyle='none',marker='o',label='regcoil')
    plt.xlabel(scanVariableName)
    plt.legend()
    plt.ylabel('d' + plotVariableName + 'd' + scanVariableName)

    plt.figure(facecolor='white')
    plt.title('Residual')
    plt.plot(scanVariables,residual,linestyle='none',marker='o')
    plt.xlabel(scanVariableName)
    plt.axvline(baseCase,color='black',linestyle='dashed')
    plt.ylabel(plotVariableName + ' residual')

  plt.figure(facecolor='white')
  plt.title('Finite Difference')
  plt.plot(scanVariables, finitediffdplotVariabledscanVariable,linestyle='none',marker='o')
  plt.xlabel(scanVariableName)
  plt.axvline(baseCase,color='black',linestyle='dashed')
  plt.ylabel('d' + plotVariableName + 'd' + scanVariableName)


  plt.show()

if __name__ == "__main__":
  if len(sys.argv) != 2:
    print "Error! You must specify 1 arguments: regcoil_in.XXX."
    exit(1)
  regcoilScanPlot(sys.argv[1])
