#!/usr/bin/env python

# These tests are mostly copied from the example regcoilPaper_figure10d_originalAngle_loRes.

# This python script checks the output file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main REGCOIL directory.

exec(open('../testsCommon.py').read())
absoluteTolerance = 1e-100

numFailures = 0

f = readOutputFile()

variableName = 'lambda'
data = f.variables[variableName][()]
relativeTolerance = 1e-12
numFailures += arrayShouldBe(data, [0, 1e-15, 1.33352143216332e-15, 1.77827941003892e-15,\
                                    2.37137370566166e-15, 3.16227766016838e-15, 4.21696503428582e-15,\
                                    5.62341325190349e-15, 7.49894209332456e-15, 1e-14],relativeTolerance,absoluteTolerance)


variableName = 'chi2_B'
data = f.variables[variableName][()]
relativeTolerance = 0.03
# Skip the lambda=0 case, which is pathological, and make sure the other ones are within a few % of the high-res results:
numFailures += arrayShouldBe(data[1:], [0.174519878306313, 0.233148551834137, \
    0.310446213587783, 0.411565852326815, 0.542616084841207, \
    0.710568315783369, 0.922890554110523, 1.18673864681126, 1.50770781013225],relativeTolerance,absoluteTolerance)

variableName = 'chi2_K'
data = f.variables[variableName][()]
relativeTolerance = 0.01
numFailures += arrayShouldBe(data[1:], [1.74957088182873e+15, 1.6989658516062e+15,\
    1.64892520707378e+15, 1.59982527129183e+15, 1.55209554042516e+15,\
    1.50621123410416e+15, 1.46269693947801e+15, 1.42212837760326e+15,\
    1.38509938973785e+15], relativeTolerance,absoluteTolerance)

variableName = 'max_Bnormal'
data = f.variables[variableName][()]
relativeTolerance = 0.03
numFailures += arrayShouldBe(data[1:], [0.117555544291034, 0.134561224876993,\
    0.153923792979921, 0.176347144447289, 0.201603899107085,\
    0.229977435728317, 0.261839981656698, 0.296625161466877, 0.333988791941713], relativeTolerance,absoluteTolerance)

variableName = 'max_K'
data = f.variables[variableName][()]
relativeTolerance = 0.06
numFailures += arrayShouldBe(data[1:], [ 8824481.60102707, 8355450.23653953,\
    7884212.39105115, 7412969.17310808, 6945281.57578389, 6485916.75531679,\
    6040435.92861007, 5614638.75966927, 5214018.48872153], relativeTolerance,absoluteTolerance)

variableName = 'single_valued_current_potential_mn'
data = f.variables[variableName][()]
#print data.shape
# We cannot exactly match data[0,:] with the example regcoilPaper_figure10d_originalAngle_loRes,
# since this data corresponds to no regularization and so is hyper-sensitive to differences in the offsetting algorithm.
# But we can match data[1,:], which corresponds to a sensible amount of regularization.
# Due to slight differences in the offset surfaces, we need pretty generous tolerances.
relativeTolerance = 0.03
absoluteTolerance = 800
#numFailures += arrayShouldBe(data[1,:], [ ],relativeTolerance,absoluteTolerance,requireSameLength=False)


del data
f.close()
print("numFailures:",numFailures)
exit(numFailures > 0)
