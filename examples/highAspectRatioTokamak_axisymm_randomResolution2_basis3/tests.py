#!/usr/bin/env python

# This python script checks the output file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main BDISTRIB directory.

execfile('../testsCommon.py')

desiredTolerance = 0.001

numFailures = 0

f = readOutputFile()

variableName = 'svd_s_transferMatrix'
data = f.variables[variableName][()]

# Actual values returned by the code:
#numFailures += arrayShouldBe(data[0,:], [0.500107470804363, 0.250061426291945, 0.125043886863589, \
#    0.0625244758101717, 0.0312623395374004, 0.0156307958627311, \
#    0.00781502680267965, 0.00390723848185482, 0.00195343541441372, \
#    0.000976600826142138, 0.000488228275983158, 0.000244070496262005, \
#    0.000122009230538231, 6.09892799306919e-05, 3.04856349812489e-05, \
#    1.52182666772794e-05], desiredTolerance)


#numFailures += shouldBe(data[0,0], 0.500, desiredTolerance)
#numFailures += arrayShouldBe(data[0,0], 0.500, desiredTolerance)
#numFailures += arrayShouldBe(data[0,:], 0.500, desiredTolerance, requireSameLength=False)

# Compare to analytically expected values:
desiredTolerance = 0.003
analyticalResults = [0.5 ** (i+1) for i in range(17)]
numFailures += arrayShouldBe(data[0,:], analyticalResults, desiredTolerance, requireSameLength=False)

f.close()
exit(numFailures > 0)
