#!/usr/bin/env python
import os, inspect, math, subprocess
from regcoilScan import readVariable, logspace_odd, namelistLineContains
from nescinUtilities import nescinReadValue, nescinWriteGeomScan
import sys
import numpy as np
from shutil import copyfile

def regcoilGeomScan4(filename):
  outputFileName = "regcoil_out" + filename[10::]
  coilFouriermmin = readVariable("coilFouriermmin","int",filename,required=False)
  coilFouriermmax = readVariable("coilFouriermmax","int",filename,required=False)
  coilFouriernmin = readVariable("coilFouriernmin","int",filename,required=False)
  coilFouriernmax = readVariable("coilFouriernmax","int",filename,required=False)
  coilFourierMinFactor = readVariable("coilFourierMinFactor","float",filename,required=False)
  coilFourierMaxFactor = readVariable("coilFourierMaxFactor","float",filename,required=False)
  coilFourierNumRuns = readVariable("coilFourierNumRuns","int",filename,required=False)
  nescin_filename = readVariable("nescin_filename","string",filename,required=False)
  wout_filename = readVariable("wout_filename","string",filename,required=False)
  new_wout_filename = "'../" + wout_filename + "'"

  for whichFourier in ("rmnc","zmns"):
    for whichM in range(coilFouriermmin, coilFouriermmax+1):
      for whichN in range(coilFouriernmin, coilFouriernmax+1):
        omega = nescinReadValue(nescin_filename,whichM,whichN,whichFourier)
        print omega
        omegas = np.linspace(coilFourierMinFactor*omega,coilFourierMaxFactor*omega,coilFourierNumRuns)
        descriptions = []
        for i in range(len(omegas)):
          descriptions.append(whichFourier + "." + str(whichM) + "." + str(whichN) + "." + str(omegas[i]))
        descriptions.append("baseCase")
        # Append basecase
        omegas = np.append(omegas,omega)
        parametersForScan = list(omegas)
        numRunsInScan = coilFourierNumRuns + 1
        print len(descriptions)
        print len(omegas)
        print descriptions
        print omegas
        
        runNum = 0
        
        while runNum < numRunsInScan:
            directory = descriptions[runNum]
            if os.path.exists(directory):
                print "Warning: directory "+directory+" already exists, so skipping this run."
                numRunsInScan -= 1
                descriptions.pop(runNum)
                parametersForScan.pop(runNum)
                runNum -= 1
            runNum += 1
                              
        print
        print "Performing scan for m = "+str(whichM)+", n = "+str(whichN)+", whichFourier = "+whichFourier
        print "Here are the parameters for the "+str(numRunsInScan)+" runs we will launch:"
        print "[omega]"
        print "-----------------------------------------------------------------------"
        for line in parametersForScan:
            print line
        print
        print "Here are the directories that will be created:"
        print descriptions

        while True:
            proceed=raw_input("Should I go ahead and launch these "+str(numRunsInScan)+" jobs? [y/n] ")
            if proceed=="y" or proceed=="n":
                break
            print "You must enter either y or n."
        if proceed=="n":
            exit(0)

        print "launching jobs..."
                              
        with open(filename, 'r') as f:
          inputFile = f.readlines()
                              
        for runNum in range(numRunsInScan):
          directory = descriptions[runNum]
          print "Beginning to handle job "+str(runNum+1)+" of "+str(numRunsInScan)+": "+directory
          os.makedirs(directory)
          os.chdir(directory)
          f = open(filename,"w")
          for line in inputFile:
            if namelistLineContains(line,"omega"):
              line = "omega = "+str(parametersForScan[runNum])+" ! Set by regcoilScan.\n"
            if namelistLineContains(line,"nescin_filename"):
              line = "nescin_filename = '" + nescin_filename + "." + whichFourier + str(whichM) + str(whichN) + "_" + str(parametersForScan[runNum]) + "' ! Set by regcoilScan.\n"
            if namelistLineContains(line,"wout_filename"):
              line = "wout_filename = " + new_wout_filename + " ! Set by regcoilScan.\n"
            f.write(line)
          f.close()
          # Copy nescin file to subdirectory
          src = "../" + nescin_filename
          dst = nescin_filename
          copyfile(src,dst)
          # Edit nescin file
          nescinWriteGeomScan(dst,whichM,whichN,whichFourier,parametersForScan[runNum])
          submitCommand = "/Users/elizabethpaul/Documents/Research/regcoil_sensitivity/regcoil/regcoil " + filename
          g = open(outputFileName,"w")
          try:
              submissionResult = subprocess.call(submitCommand.split(" "),stdout=g)
          except: 
              print "ERROR: Unable to submit run "+directory+" for some reason."
              raise
          else:
              if submissionResult==0:
                  print "No errors submitting job "+directory
              else:
                  print "Nonzero exit code returned when trying to submit job "+directory
          g.close()
          os.chdir('..')

if __name__ == "__main__":
  if len(sys.argv) != 1:
    print "Error! You must specify 1 arguments: regcoil_in.XXX."
    exit(1)
  regcoilGeomScan4(sys.argv[1])



