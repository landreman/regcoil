#!/usr/bin/env python
import os, inspect, math, subprocess
from regcoilScan import readVariable, logspace_odd, namelistLineContains
from nescinUtilities import nescinReadValue, nescinWriteGeomScan
from regcoilGeomScan4 import regcoilGeomScan4
import sys
import numpy as np
from shutil import copyfile

# Input filename
filename = sys.argv[1]

outputFileName = "regcoil_out" + filename[10::]
scanType = readVariable("scanType","int",filename,required=False)
if (scanType < 4):
  a_coil = readVariable("a_coil","float",filename,required=False)
  R0_coil = readVariable("R0_coil","float",filename,required=False)
  R0_plasma = readVariable("R0_plasma","float",filename,required=False)
if (scanType == 1):
  a_coilMinFactor = readVariable("a_coilMinFactor","float",filename,required=False)
  a_coilMaxFactor = readVariable("a_coilMaxFactor","float",filename,required=False)
  a_coilNumRuns = readVariable("a_coilNumRuns","int",filename,required=False)
  regcoilGeomScan1(filename,a_coilMinFactor,a_coilMaxFactor,a_coilNumRuns)
if (scanType == 2):
  R0_coilMinFactor = readVariable("R0_coilMinFactor","float",filename,required=False)
  R0_coilMaxFactor = readVariable("R0_coilMaxFactor","float",filename,required=False)
  R0_coilNumRuns = readVariable("R0_coilNumRuns","float",filename,required=False)
  regcoilGeomScan2(filename,R0_coilMinFactor,R0_coilMaxFactor,R0_coilNumRuns)
if (scanType == 3):
  R0_plasmaMinFactor = readVariable("R0_plasmaMinFactor","float",filename,required=False)
  R0_plasmaMaxFactor = readVariable("R0_plasmaMaxFactor","float",filename,required=False)
  R0_plasmaNumRuns = readVariable("R0_plasmaNumRuns","float",filename,required=False)
  regcoilGeomScan3(filename,R0_plasmaMinFactor,R0_plasmaMaxFactor,R0_plasmaNumRuns)
if (scanType == 4):
  regcoilGeomScan4(filename)

#if (scanType == 1):
#  a_coils = np.linspace(a_coilMinFactor*a_coil,a_coilMaxFactor*a_coil,a_coilNumRuns)
#else:
#  a_coils = np.zeros((0,0))
#if (scanType == 2):
#  R0_coils = np.linspace(R0_coilMinFactor*R0_coil,R0_coilMaxFactor*R0_coil,R0_coilNumRuns)
#else:
#  R0_coils = np.zeros((0,0))
#if (scanType == 3):
#  R0_plasmas = np.linspace(R0_plasmaMinFactor*R0_plasma,R0_plasmaMaxFactor*R0_plasma,R0_plasmaNumRuns)
#else:
#  R0_plasmas = np.zeros((0,0))
#if (scanType == 4):
#  omegas = np.linspace(coilFourierMinFactor*omega,coilFourierMaxFactor*omega,coilFourierNumRuns)
#else:
#  omegas = np.zeros((0,0))
#
#numRunsInScan = 1+len(a_coils)+len(R0_plasmas)+len(R0_coils)+len(omegas)
#if (scanType < 4):
#  baseCase = [a_coil,R0_plasma,R0_coil]
#else:
#  baseCase = [omega]
#
#parametersForScan = []
#for i in range(numRunsInScan):
#    parametersForScan.append(list(baseCase))
#
#currentIndex = 1
#descriptions = ["baseCase"]
#
#if (scanType < 4):
#  for i in range(len(a_coils)):
#      parametersForScan[currentIndex][0] = a_coils[i]
#      descriptions.append("a_coil" + str(a_coils[i]))
#      currentIndex += 1
#  for i in range(len(R0_plasmas)):
#      parametersForScan[currentIndex][1] = R0_plasmas[i]
#      descriptions.append("R0_plasma" + str(R0_plasmas[i]))
#      currentIndex += 1
#  for i in range(len(R0_coils)):
#      parametersForScan[currentIndex][2] = R0_coils[i]
#      descriptions.append("R0_coil" + str(R0_coils[i]))
#      currentIndex += 1
#else:
#  for i in range(len(omegas)):
#      parametersForScan[currentIndex][0] = omegas[i]
#      descriptions.append(coilWhichFourier + "." + str(coilFourierm) + "." + str(coilFouriern) + "." + str(omegas[i]))
#      currentIndex += 1
#if currentIndex != numRunsInScan:
#    print "Error! Something went wrong."
#    exit(1)
#if len(parametersForScan) != len(descriptions):
#    print "Error! Something went wrong."
#    exit(1)
#
#runNum = 0
#while runNum < numRunsInScan:
#    directory = descriptions[runNum]
#    if os.path.exists(directory):
#        print "Warning: directory "+directory+" already exists, so skipping this run."
#        numRunsInScan -= 1
#        descriptions.pop(runNum)
#        parametersForScan.pop(runNum)
#        runNum -= 1
#    runNum += 1
#
#print
#print "Here are the parameters for the "+str(numRunsInScan)+" runs we will launch:"
#print "[a_coil, R0_plasma, R0_coil, omega]"
#print "-----------------------------------------------------------------------"
#for line in parametersForScan:
#    print line
#print
#print "Here are the directories that will be created:"
#print descriptions
#
#while True:
#    proceed=raw_input("Should I go ahead and launch these "+str(numRunsInScan)+" jobs? [y/n] ")
#    if proceed=="y" or proceed=="n":
#        break
#    print "You must enter either y or n."
#if proceed=="n":
#    exit(0)
#
#print "launching jobs..."
#
#with open(filename, 'r') as f:
#    inputFile = f.readlines()
#
#for mi in range(coilFouriermmin,coilFouriermmax):
#  for ni in range(coilFouriernmin,coilFouriernmax):
#  for runNum in range(numRunsInScan):
#      directory = descriptions[runNum]
#      print "Beginning to handle job "+str(runNum+1)+" of "+str(numRunsInScan)+": "+directory
#      os.makedirs(directory)
#      os.chdir(directory)
#      f = open(filename,"w")
#      for line in inputFile:
#          if (scanType < 4):
#            if namelistLineContains(line,"a_coil"):
#                line = "a_coil = "+str(parametersForScan[runNum][0])+" ! Set by regcoilScan.\n"
#            if namelistLineContains(line,"R0_plasma"):
#                line = "R0_plasma = "+str(parametersForScan[runNum][1])+" ! Set by regcoilScan.\n"
#            if namelistLineContains(line,"R0_coil"):
#                line = "R0_coil = "+str(parametersForScan[runNum][2])+" ! Set by regcoilScan.\n"
#          if (scanType == 4):
#            if namelistLineContains(line,"omega"):
#              line = "omega = "+str(parametersForScan[runNum][0])+" ! Set by regcoilScan.\n"
#            if namelistLineContains(line,"nescin_filename"):
#              line = "nescin_filename = '" + nescin_filename + "." + coilWhichFourier + str(coilFourierm) + str(coilFouriern) + "_" + str(parametersForScan[runNum][0]) + "' ! Set by regcoilScan.\n"
#          if namelistLineContains(line,"wout_filename"):
#              new_wout_filename = "../" + wout_filename
#              line = "wout_filename = " + new_wout_filename + " ! Set by regcoilScan.\n"
#          f.write(line)
#      f.close()
#      if (scanType == 4):
#        # Copy nescin file to subdirectory
#        src = "../" + nescin_filename
#        dst = nescin_filename
#        copyfile(src,dst)
#        # Edit nescin file
#        nescinWriteGeomScan(dst,coilFourierm,coilFouriern,coilWhichFourier,parametersForScan[runNum][0])
#        nescinWriteGeomScan(dst,coilFou)
#      submitCommand = "/Users/elizabethpaul/Documents/Research/regcoil_sensitivity/regcoil/regcoil " + filename
#      g = open(outputFileName,"w")
#      try:
#          submissionResult = subprocess.call(submitCommand.split(" "),stdout=g)
#      except: 
#          print "ERROR: Unable to submit run "+directory+" for some reason."
#          raise
#      else:
#          if submissionResult==0:
#              print "No errors submitting job "+directory
#          else:
#              print "Nonzero exit code returned when trying to submit job "+directory
#      g.close()
#      os.chdir('..')

