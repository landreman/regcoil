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

numFailures += arrayShouldBe(data[0,:], [0.500107470804363, 0.498701861460618, 0.498570256474812, 0.494757120411698, \
    0.494526925087883, 0.250061426291945, 0.249914764213506, \
    0.24974631585096, 0.248896865997166, 0.248805348743705, \
    0.246437256346953, 0.242361254290523, 0.125043886863589, \
    0.124945329157845, 0.12494529756196, 0.124650368836285, \
    0.124650367413053, 0.0625244758101713, 0.0624889695036331, \
    0.0624889694671746, 0.0623826230023029, 0.0623826229993291, \
    0.0312623395373137, 0.0312485243982378, 0.0312485243982345, \
    0.0312071283860475, 0.0312071283860459, 0.0156307956886087, \
    0.015625154225399, 0.015625154225399, 0.0156082456386442, \
    0.0156082456386441, 0.00781500687216377, 0.00781262576668433, \
    0.00781262576668427, 0.00780548789863795, 0.00780548789863794, \
    0.0039028921106597, 0.00390186714317446, 0.00390186714317442, \
    0.00389879420871117, 0.00389879420871117], desiredTolerance)

#numFailures += shouldBe(data[0,0], 0.500, desiredTolerance)
#numFailures += arrayShouldBe(data[0,0], 0.500, desiredTolerance)
#numFailures += arrayShouldBe(data[0,:], 0.500, desiredTolerance, requireSameLength=False)


f.close()
print "numFailures: ",numFailures
exit(numFailures > 0)
