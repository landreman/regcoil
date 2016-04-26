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

numFailures += arrayShouldBe(data[0,:], [  0.44336924002968, 0.316270310159807, 0.269141744576492, 0.181741322974342, \
    0.175749554360671, 0.148499915381192, 0.135140685503816,    \
    0.126144947487883, 0.114859716729841, 0.102127446816665,    \
    0.0878704132354457, 0.0767381041440996, 0.072921065706652,  \
    0.0687813764094384, 0.0667612687441194, 0.0609444216218977, \
    0.0542356116843527, 0.0512526801576164, 0.0490485431598584, \
    0.0463251749503744, 0.0421386677214716, 0.0412174829985826], desiredTolerance, requireSameLength=False)


#numFailures += shouldBe(data[0,0], 0.500, desiredTolerance)
#numFailures += arrayShouldBe(data[0,0], 0.500, desiredTolerance)
#numFailures += arrayShouldBe(data[0,:], 0.500, desiredTolerance, requireSameLength=False)

f.close()
exit(numFailures > 0)
