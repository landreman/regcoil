#!/usr/bin/env python

# This python script checks the output file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main REGCOIL directory.

execfile('../testsCommon.py')
absoluteTolerance = 1e-100

numFailures = 0

f = readOutputFile()



variableName = 'exit_code'
data = f.variables[variableName][()]
relativeTolerance = 1e-12
numFailures += shouldBe(data, -2,relativeTolerance,absoluteTolerance)


variableName = 'lambda'
data = f.variables[variableName][()]
relativeTolerance = 1e-12
numFailures += arrayShouldBe(data, [1e+200],relativeTolerance,absoluteTolerance)



f.close()
print "numFailures:",numFailures
exit(numFailures > 0)
