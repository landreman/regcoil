#!/usr/bin/env python

# This python script checks the output file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main BDISTRIB directory.

execfile('../testsCommon.py')

numFailures = 0

f = readOutputFile()

variableName = 'svd_s_transferMatrix'
data = f.variables[variableName][()]
a_plasma = f.variables['a_plasma'][()]
a_middle = f.variables['a_middle'][()]
R = f.variables['R0_middle'][()]

from scipy.special import iv
# iv is the modified Bessel function I_v

svd_s_transferMatrix = data
max_n = len(svd_s_transferMatrix[0,:])
analytical = []
for n in range(1,max_n+1):
    s = (a_plasma/a_middle)*(iv(-1,n*a_plasma/R)+iv(1,n*a_plasma/R))/(iv(-1,n*a_middle/R)+iv(1,n*a_middle/R))
    analytical.append(s)

# Compare to analytically expected values:
desiredTolerance = 0.015
numFailures += arrayShouldBe(data[0,:], analytical, desiredTolerance)

f.close()
exit(numFailures > 0)
