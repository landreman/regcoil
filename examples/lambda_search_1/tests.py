#!/usr/bin/env python

# This python script checks the output file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main REGCOIL directory.

execfile('../testsCommon.py')
absoluteTolerance = 1e-100

numFailures = 0

f = readOutputFile()

variableName = 'chi2_B_target'
data = f.variables[variableName][()]
relativeTolerance = 1e-4
numFailures += shouldBe(data, 0.33784814179529,relativeTolerance,absoluteTolerance)


variableName = 'exit_code'
data = f.variables[variableName][()]
relativeTolerance = 1e-12
numFailures += shouldBe(data, 0,relativeTolerance,absoluteTolerance)


variableName = 'lambda'
data = f.variables[variableName][()]
relativeTolerance = 1e-10
numFailures += arrayShouldBe(data, [1e+200, 0, 1.86103543888473e-17, 1.86103543888473e-15, \
    1.67447016215951e-15, 1.72798458404329e-15, 1.72841562509787e-15, \
    1.72840698304124e-15],relativeTolerance,absoluteTolerance)


variableName = 'chi2_B'
data = f.variables[variableName][()]
relativeTolerance = 0.001
numFailures += arrayShouldBe(data, [13.149476484944, 7.69174650558559e-06, 0.00423827116265938, \
    0.362892488869994, 0.327618094383163, 0.337768127590254, \
    0.337849778826274, 0.33784814179529],relativeTolerance,absoluteTolerance)

variableName = 'chi2_K'
data = f.variables[variableName][()]
relativeTolerance = 0.001
numFailures += arrayShouldBe(data, [1.21295356557164e+15, 1.9912738045989e+16, 2.57210417536396e+15, \
    1.67310511685225e+15, 1.69307959624817e+15, 1.68711276457985e+15, \
    1.68706551817808e+15, 1.68706646530854e+15], relativeTolerance,absoluteTolerance)

variableName = 'max_Bnormal'
data = f.variables[variableName][()]
relativeTolerance = 0.001
numFailures += arrayShouldBe(data, [0.800294392804433, 0.00131081772888197, 0.0344837973451904, 
    0.170106680639416, 0.161770814559884, 0.164214867896084, 
    0.164234374507914, 0.164233983442413], relativeTolerance,absoluteTolerance)

variableName = 'max_K'
data = f.variables[variableName][()]
relativeTolerance = 0.001
numFailures += arrayShouldBe(data, [2713725.09530369, 125251911.209827, 15467287.3730328, 
    7877128.17061528, 8052687.63164847, 8000414.46145111, 7999999.98035358, 
    8000008.28938568], relativeTolerance,absoluteTolerance)

f.close()
print "numFailures:",numFailures
exit(numFailures > 0)
