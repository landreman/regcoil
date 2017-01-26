def readOutputFile():
    import os

    head, dirname = os.path.split(os.getcwd())
    outputFilename = "regcoil_out."+dirname+".nc"
    
    if not os.path.isfile(outputFilename):
        print "Error! The output file "+outputFilename+" has not been created."
        exit(1)
        
    from scipy.io import netcdf
    try:
        f = netcdf.netcdf_file(outputFilename,'r',mmap=False)
    except:
        print "ERROR! Unable to read netCDF output file "+outputFilename
        raise

    print "Reading output file "+outputFilename
    return f

def shouldBe(latestValue, trueValue, relativeTolerance, absoluteTolerance):
    difference = abs(latestValue-trueValue)
    if abs(trueValue) > 0:
        relativeDifference = abs(difference / trueValue)
        relativeTest = (relativeDifference <= relativeTolerance)
    else:
        relativeTest = False
    absoluteTest = (difference <= absoluteTolerance)
    string = "Variable "+variableName+" should be close to "+str(trueValue)+", and it was "+str(latestValue)
    if relativeTest:
        if absoluteTest:
            print "    Test passed. "+string+". Both abs and rel tol met."
            return 0
        else:
            print "    Test passed. "+string+". Rel tol met. Abs tol not met."
            return 0
    else:
        if absoluteTest:
            print "    Test passed. "+string+". Abs tol met. Rel tol not met."
            return 0
        else:
            print "*** TEST FAILED! "+string+". Neither rel nor abs tol met."
            return 1


def arrayShouldBe(latestValues, trueValues, relativeTolerance, absoluteTolerance, requireSameLength = True):
    # These next few lines are a hack so this function can be called on scalars without an exception
    try:
        temp = len(latestValues)
    except:
        print "arrayifying latestValues"
        latestValues = [latestValues]

    try:
        temp = len(trueValues)
    except:
        print "arrayifying trueValues"
        trueValues = [trueValues]

    if requireSameLength and (len(latestValues) != len(trueValues)):
        print "*** TEST FAILED!! Variable "+variableName+" should have length "+str(len(trueValues))+" but it instead has length "+str(len(latestValues))
        return 1

    if len(latestValues) < len(trueValues):
        print "*** TEST FAILED!! Variable "+variableName+" should have length at least "+str(len(trueValues))+" but it instead has length "+str(len(latestValues))
        return 1

    numArrayErrors = 0
    for i in range(len(trueValues)):
        numArrayErrors += shouldBe(latestValues[i],trueValues[i],relativeTolerance, absoluteTolerance)

    return numArrayErrors
