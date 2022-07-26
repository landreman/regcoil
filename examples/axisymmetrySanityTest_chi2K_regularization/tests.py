#!/usr/bin/env python

# This python script checks the output file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main REGCOIL directory.

# In this example, the plasma and coil surfaces are both axisymmetric,
# so the single valued part of the current potential and B_normal should vanish to
# machine precision (even though the plasma and coil surfaces have different major radius.)

exec(open('../testsCommon.py').read())
absoluteTolerance = 1e-10
relativeTolerance = 1e-100 # The relative tolerance is irrelevant since the true values are 0.

numFailures = 0

f = readOutputFile()

variableName = 'chi2_B'
data = f.variables[variableName][()]
numFailures += arrayShouldBe(data, [0,0,0], relativeTolerance,absoluteTolerance)

variableName = 'max_Bnormal'
data = f.variables[variableName][()]
numFailures += arrayShouldBe(data, [0,0,0], relativeTolerance,absoluteTolerance)

variableName = 'single_valued_current_potential_mn'
data = f.variables[variableName][()]
#print data.shape
numFailures += arrayShouldBe(data[0,:], [0]*97, relativeTolerance,absoluteTolerance)
numFailures += arrayShouldBe(data[1,:], [0]*97, relativeTolerance,absoluteTolerance)
numFailures += arrayShouldBe(data[2,:], [0]*97, relativeTolerance,absoluteTolerance)


del data
f.close()
print("numFailures:",numFailures)
exit(numFailures > 0)
