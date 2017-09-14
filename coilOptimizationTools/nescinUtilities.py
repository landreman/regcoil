import re
import sys
import numpy as np
import os

def nescinReadValue(inputFile,m,n,geom):

  if re.match(geom, "rmnc"):
    iOmega = 0
  elif re.match(geom, "zmns"):
    iOmega = 1
  elif re.match(geom, "rmns"):
    iOmega = 2
  elif re.match(geom, "zmnc"):
    iOmega = 3
  else:
    print "Invalid geom parameter"
    exit

  file = open(inputFile,"r")
  inCurrentGeometry = 0
  # m,n,crc2,czs2,crs2,czc2
  ms = []
  ns = []
  rmncs = []
  zmnss = []
  rmnss = []
  zmncs = []
  for line in file:
      if (inCurrentGeometry):
          list = line.split()
          mCurr = int(list[0])
          nCurr = int(list[1])
          if (mCurr == m and nCurr == n):
            return float(list[iOmega+2])
      if re.match("------ Current Surface",line):
          inCurrentGeometry = 1
          next(file)
          NumModes = int(next(file))
          next(file)
          next(file)
  file.close()
  print "Error in reading desired nescin value."

def nescinWriteGeomScan(inputFile,m,n,geom,omega):
  #print os.system('pwd')
  #print os.system('ls')
  file = open(inputFile,"r")
  inCurrentGeometry = 0
  # m,n,crc2,czs2,crs2,czc2
  ms = []
  ns = []
  rmncs = []
  zmnss = []
  rmnss = []
  zmncs = []
  for line in file:
    if (inCurrentGeometry):
      list = line.split()
      ms.append(int(list[0]))
      ns.append(int(list[1]))
      rmncs.append(float(list[2]))
      zmnss.append(float(list[3]))
      rmnss.append(float(list[4]))
      zmncs.append(float(list[5]))
    if re.match("------ Current Surface",line):
      inCurrentGeometry = 1
      next(file)
      NumModes = int(next(file))
      next(file)
      next(file)
  file.close()
    
  #geom = raw_input("Scan over rmnc, rmns, zmnc, or zmns: ")
  # m = int(raw_input("m to use for scan: "))
  if (m < min(ms) or m > max(ms)):
    print "Choice of m out of range"
    exit
  # n = int(raw_input("n to use for scan: "))
  if (abs(n) < min(ns) or abs(n) > max(ns)):
    print "Choice of n out of range"
    exit
    # minFactor = float(raw_input("minFactor: "))
    # maxFactor = float(raw_input("maxFactor: "))
    # NumRuns = int(raw_input("NumRuns: "))

    #for i in range(len(ms)):
    #     if (ms[i] == m and ns[i] == n):
    #         modeToScan = i

  if re.match(geom, "rmnc"):
    #base = rmncs[i]
    whichOmega = 0
    fileName = inputFile + ".rmnc" + str(m) + str(n) + "_" + str(omega)
  if re.match(geom, "zmns"):
    # base = zmnss[i]
    whichOmega = 1
    fileName = inputFile + ".zmns" + str(m) + str(n) + "_" + str(omega)
  if re.match(geom, "rmns"):
    # base = rmnss[i]
    whichOmega = 2
    fileName = inputFile + ".rmns" + str(m) + str(n) + "_" + str(omega)
  if re.match(geom, "zmnc"):
    # base = zmncs[i]
    whichOmega = 3
    fileName = inputFile + ".zmnc" + str(m) + str(n) + "_" + str(omega)
    
    
  #omegas = np.linspace(base*minFactor, base*maxFactor, NumRuns)
  currentMode = -1
  # for omega in omegas:
  rmnc_write = rmncs
  rmns_write = rmnss
  zmnc_write = zmncs
  zmns_write = zmnss

  newFile = open(fileName, "w")
  oldFile = open(inputFile,"r")
  for line in oldFile:
    if (currentMode >= 0):
      if (ms[currentMode] == m and ns[currentMode] == n):
        if (whichOmega == 0):
          rmnc_write[currentMode] = omega
        if (whichOmega == 1):
          zmns_write[currentMode] = omega
        if (whichOmega == 2):
          rmns_write[currentMode] = omega
        if (whichOmega == 3):
          zmnc_write[currentMode] = omega
        lineToWrite = "\t" + str(ms[currentMode]) + "\t" + str(ns[currentMode]) + "\t" + str(rmnc_write[currentMode]) \
          + "\t" + str(zmns_write[currentMode]) + "\t" + str(rmns_write[currentMode]) + "\t" \
          + str(zmnc_write[currentMode]) + "\n"
        newFile.write(lineToWrite)
        currentMode = currentMode + 1
      else:
        currentMode = currentMode + 1
        newFile.write(line)
    elif re.match("------ Current Surface",line):
      newFile.write(line)
      lineToWrite = next(oldFile)
      newFile.write(lineToWrite)
      lineToWrite = next(oldFile)
      newFile.write(lineToWrite)
      lineToWrite = next(oldFile)
      newFile.write(lineToWrite)
      lineToWrite = next(oldFile)
      newFile.write(lineToWrite)
      currentMode = 0
    else:
      newFile.write(line)

  newFile.close()
  oldFile.close()
  os.remove(inputFile)
