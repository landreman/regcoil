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

numFailures += arrayShouldBe(data[0,:], [0.597902809226463, 0.565092889833134, 0.547199036539407, 0.470147029233092, \
    0.446141257518072, 0.391169410090854, 0.357065110516946, \
    0.351757082247666, 0.339196355927458, 0.335177228382406, \
    0.331322772121299, 0.304230491334018, 0.293431178465185, \
    0.292095300908778, 0.291076938451873, 0.282782114936251, \
    0.263149531031626, 0.248938060273492, 0.241762937876781], desiredTolerance, requireSameLength=False)


#numFailures += shouldBe(data[0,0], 0.500, desiredTolerance)
#numFailures += arrayShouldBe(data[0,0], 0.500, desiredTolerance)
#numFailures += arrayShouldBe(data[0,:], 0.500, desiredTolerance, requireSameLength=False)

f.close()
exit(numFailures > 0)
