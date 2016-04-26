def readOutputFile():
    import os

    head, dirname = os.path.split(os.getcwd())
    outputFilename = "bdistrib_out."+dirname+".nc"
    
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

def shouldBe(latestValue, trueValue, relativeTolerance):
    relativeDifference = abs((latestValue - trueValue) / trueValue)
    if relativeDifference > relativeTolerance:
        print "*** TEST FAILED!!  Variable "+variableName+" should be close to "+str(trueValue)+", but it is instead "+str(latestValue)
        print "Actual / correct = ",latestValue/trueValue
        return 1
    else:
        print "    Test passed:   Variable "+variableName+" should be close to "+str(trueValue)+", and it came out to be "+str(latestValue)+", which is within tolerance."
        return 0

def arrayShouldBe(latestValues, trueValues, relativeTolerance, requireSameLength = True):
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
        numArrayErrors += shouldBe(latestValues[i],trueValues[i],relativeTolerance)

    return numArrayErrors
