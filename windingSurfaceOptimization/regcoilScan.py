import os, inspect, math, subprocess
from sys import argv
import string

# Copied from sfincsScan
def namelistLineContains(line,varName):
    line2 = line.strip().lower()
    varName = varName.lower()
# We need enough characters for the varName, =, and value:
    if len(line2)<len(varName)+2:
        return False
    if line2[0]=="!":
        return False
    nextChar = line2[len(varName)]
    if line2[:len(varName)]==varName and (nextChar==" " or nextChar=="="):
        return True
    else:
        return False

def readVariable(varName, intOrFloatOrString, inputFilename, required=True):
# This function reads normal fortran variables from the input.namelist file.
# It is assumed that the input.namelist file has been loaded into the variable "inputFile".
    with open(inputFilename, 'r') as f:
        inputFile = f.readlines()
    if (intOrFloatOrString != "int") and (intOrFloatOrString != "float") and (intOrFloatOrString != "string"):
	    print("intOrFloatOrString must be int, float, or string.")
	    exit(1)

    originalVarName = varName
    returnValue = None
    numValidLines = 0
    for line in inputFile:
        line3 = line.strip()
        if len(line3)<1:
	        continue
        if line3[0]=="!":
	        continue
        if len(line3) < len(varName)+2:
            continue
        if not line3[:len(varName)].lower()==varName.lower():
            continue
        line4 = line3[len(varName):].strip()
        if not line4[0] =="=":
            continue
        line5 = line4[1:].strip();
        if intOrFloatOrString != "string":
            line5 = line5.replace('d','e').replace('D','e')
    # Remove any comments:
        if "!" in line5:
            line5 = line5[:string.find(line5,"!")]
        if intOrFloatOrString=="int":
            try:
                returnValue = int(line5)
                numValidLines += 1
            except:
                print("Warning! I found a definition for the variable "+originalVarName+" in "+inputFilename+" but I was unable to parse the line to get an integer." #Added by AM 2015-12)
                print("Here is the line in question:")
                print(line)
        elif intOrFloatOrString=="float":
            try:
                returnValue = float(line5)
                numValidLines += 1
            except:
                print("Warning! I found a definition for the variable "+originalVarName+" in "+inputFilename+" but I was unable to parse the line to get a float." #Added by AM 2015-12)
                print("Here is the line in question:")
                print(line)
        elif intOrFloatOrString=="string":
          # Strip quotes
          if (line5.startswith('"') and line5.endswith('"')) or (line5.startswith("'") and line5.endswith("'")):
            line5 = line5[1:-1]
          returnValue = line5
          numValidLines += 1

    if required and returnValue==None:
        print("Error! Unable to find a valid setting for the variable "+originalVarName+" in "+inputFilename+".")
        exit(1) #Added by AM 2015-12

    if numValidLines > 1:
        print("Warning! More than 1 valid definition was found for the variable "+originalVarName+". The last one will be used.")

#print("Read "+originalVarName+" = "+str(returnValue))
    return returnValue

def logspace_odd(min,max,nn):
    temp = map(int,logspace(min,max,nn))
    temp2 = []
    for x in temp:
        if (x % 2 == 0):
            temp2.append(x+1)
        else:
            temp2.append(x)
    return uniq(temp2)

def namelistLineContains(line,varName):
    line2 = line.strip().lower()
    varName = varName.lower()
# We need enough characters for the varName, =, and value: 
    if len(line2)<len(varName)+2:
        return False
    if line2[0]=="!":
        return False
    nextChar = line2[len(varName)]
    if line2[:len(varName)]==varName and (nextChar==" " or nextChar=="="):
        return True
    else:
        return False
