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
numFailures += shouldBe(data, 0.337712583344148 ,relativeTolerance,absoluteTolerance)


variableName = 'exit_code'
data = f.variables[variableName][()]
relativeTolerance = 1e-12
numFailures += shouldBe(data, 0,relativeTolerance,absoluteTolerance)


variableName = 'lambda'
data = f.variables[variableName][()]
relativeTolerance = 1e-9
numFailures += arrayShouldBe(data, [1e+200, 0, 1.86098014413313e-17, 1.86098014413313e-15, \
    1.67493916804206e-15, 1.7278563450837e-15, 1.72828277387722e-15, \
    1.72827413248484e-15],relativeTolerance,absoluteTolerance)


variableName = 'chi2_B'
data = f.variables[variableName][()]
relativeTolerance = 0.001
numFailures += arrayShouldBe(data, [13.1506816246147, 7.71670836033572e-06, 0.00424834281462129, \
    0.362769356668118, 0.327597924797486, 0.337633450768883, \
    0.337714220082776, 0.337712583344148],relativeTolerance,absoluteTolerance)

variableName = 'chi2_K'
data = f.variables[variableName][()]
relativeTolerance = 0.001
numFailures += arrayShouldBe(data, [1.21300287142288e+15, 1.99560519855456e+16, 2.5719665173851e+15, \
    1.6732483015856e+15, 1.69316202257536e+15, 1.68726310825612e+15, \
    1.68721636863624e+15, 1.68721731567034e+15], relativeTolerance,absoluteTolerance)

variableName = 'max_Bnormal'
data = f.variables[variableName][()]
relativeTolerance = 0.001
numFailures += arrayShouldBe(data, [0.800268516799401, 0.00131646260910012, 0.0345167565909495, \
    0.17013629890868, 0.161820866475293, 0.16423854399308, 0.164257851045158, \
    0.164257459823674], relativeTolerance,absoluteTolerance)

variableName = 'max_K'
data = f.variables[variableName][()]
relativeTolerance = 0.001
numFailures += arrayShouldBe(data, [2754022.92888632, 126095244.116546, 15504329.594164, 
    7877053.03140448, 8052093.2255605, 8000410.04678498, 7999999.98076259, 
    8000008.28957208], relativeTolerance,absoluteTolerance)

del data
f.close()
print("numFailures:",numFailures)
exit(numFailures > 0)
