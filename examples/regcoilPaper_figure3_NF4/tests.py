#!/usr/bin/env python

# This python script checks the output file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main REGCOIL directory.

execfile('../testsCommon.py')
absoluteTolerance=0

numFailures = 0

f = readOutputFile()

variableName = 'chi2_B'
data = f.variables[variableName][()]
relativeTolerance = 1e-4
numFailures += arrayShouldBe(data, [0.0106199394066055],relativeTolerance,absoluteTolerance)

variableName = 'chi2_K'
data = f.variables[variableName][()]
relativeTolerance = 1e-4
numFailures += arrayShouldBe(data, [180027864932680], relativeTolerance,absoluteTolerance)

variableName = 'max_Bnormal'
data = f.variables[variableName][()]
relativeTolerance = 1e-4
numFailures += arrayShouldBe(data, [0.0961657657267785], relativeTolerance,absoluteTolerance)

variableName = 'max_K'
data = f.variables[variableName][()]
relativeTolerance = 1e-4
numFailures += arrayShouldBe(data, [7676633.79432514], relativeTolerance,absoluteTolerance)

variableName = 'single_valued_current_potential_mn'
data = f.variables[variableName][()]
#print data.shape
relativeTolerance = 1e-4
numFailures += arrayShouldBe(data[0,:], [447.38366099137, 4639.9189494066, -1606.04570821836, -32567.0650920621,
    -40269.6076578081, 25730.5534047766, 23442.3312795834, -6320.65393830862,
    128400.032071504, 81682.2138951768, -2315.3992728423, 78024.8373333916,
    -35098.0765328628, -21023.3130109801, -38880.904286096,
    -4861.87493846693, 22920.4278347224, 111900.735160646, -349436.493977956,
    -5960.42395421097, -47528.199176313, -17603.7318829644, 31972.1106230657,
    11612.7118640209, -18751.6948378737, 13296.4775858065, -98156.5616520021,
    79971.5649827413, 104706.918633579, -105032.959996618, 36601.6555217408,
    926.735923820139, 17551.2837580923, 36130.0572868801, -14477.5458743423,
    61195.2509509964, -5825.70220383692, -48093.682021939, 46725.9456220868,
    14567.7808069598   ],relativeTolerance,absoluteTolerance,requireSameLength=False)


f.close()
print "numFailures:",numFailures
exit(numFailures > 0)
