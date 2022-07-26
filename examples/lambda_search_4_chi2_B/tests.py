#!/usr/bin/env python

# This python script checks the output file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main REGCOIL directory.

exec(open('../testsCommon.py').read())
absoluteTolerance = 1e-100

numFailures = 0

f = readOutputFile()

variableName = 'chi2_B_target'
data = f.variables[variableName][()]
relativeTolerance = 1e-4
numFailures += shouldBe(data, 0.0999995031486614 ,relativeTolerance,absoluteTolerance)


variableName = 'exit_code'
data = f.variables[variableName][()]
relativeTolerance = 1e-12
numFailures += shouldBe(data, 0,relativeTolerance,absoluteTolerance)


variableName = 'lambda'
data = f.variables[variableName][()]
relativeTolerance = 1e-10
numFailures += arrayShouldBe(data, [1e+200, 0, 1.86098014413313e-17, 1.86098018572928e-15, \
    4.90046557609783e-16, 5.05252075961572e-16, 5.05061841603372e-16, \
    5.05059316300441e-16],relativeTolerance,absoluteTolerance)


variableName = 'chi2_B'
data = f.variables[variableName][()]
relativeTolerance = 0.001
numFailures += arrayShouldBe(data, [13.1506816246147, 7.71670836033572e-06, 0.00424834281462129, \
                                        0.362769364497339, 0.0970250758584328, 0.100037690832713, \
                                        0.10000000343785, 0.0999995031486614\
                                        ],relativeTolerance,absoluteTolerance)

variableName = 'chi2_K'
data = f.variables[variableName][()]
relativeTolerance = 0.001
numFailures += arrayShouldBe(data, [1.21300287142288e+15, 1.99560519855456e+16, 2.5719665173851e+15, \
    1.67324829737856e+15, 1.93334471048967e+15, 1.9272905464157e+15, \
    1.92736515173323e+15, 1.92736614228606e+15], relativeTolerance,absoluteTolerance)

variableName = 'max_Bnormal'
data = f.variables[variableName][()]
relativeTolerance = 0.001
numFailures += arrayShouldBe(data, [0.800268516799401, 0.00131646260910012, 0.0345167565909495, \
    0.170136300712312, 0.0917700178409628, 0.0923005209753598, \
    0.0922940672862682, 0.0922939815848665], relativeTolerance,absoluteTolerance)

variableName = 'max_K'
data = f.variables[variableName][()]
relativeTolerance = 0.001
numFailures += arrayShouldBe(data, [2754022.92888632, 126095244.116546, 15504329.594164, 
    7877052.99425298, 10066575.8807533, 10017193.5931966, 10017802.3175868, 
    10017810.399735], relativeTolerance,absoluteTolerance)

del data
f.close()
print("numFailures:",numFailures)
exit(numFailures > 0)
