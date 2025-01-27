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
relativeTolerance = 1.0e-4
absoluteTolerance = 1.0e-2
numFailures += arrayShouldBe(data[1,:], [-636966.461168403,\
          -92094.356379444,\
          -12789.127051605,\
           -7010.275031314,\
         -1350.80696646224,\
         -359.820926459283,\
         -96.6165551832622,\
         -10.5401888362916,\
         -5.46898639120508,\
       -0.0365596619502098,\
        -0.108600468312161,\
      -0.00672381099875835,\
        0.0176356493852388,\
      -0.00671459645231758,\
       0.00310725404141035,\
     -0.000729744596224883,\
      8.75820409873942e-05,\
      2.99951629101753e-05,\
     -1.77759065953022e-05,\
      5.95043850654764e-06,\
     -6.21377719049132e-06,\
      3.26168547297748e-06,\
     -5.97283952293703e-07,\
     -1.71663483302622e-07,\
      2.85622196744898e-07,\
     -1.50929794496497e-07,\
     -4.91352655345434e-07,\
     -1.45207530361089e-07,\
     -5.02910632008792e-08,\
     -1.32059940462981e-07,\
      1.70121464698855e-07,\
      5.27683211533008e-08,\
     -4.20923526890395e-07,\
      1.48240105458905e-07,\
      1.79689851443011e-07,\
     -1.22649622399976e-07,\
     -2.06402762551819e-08,\
     -2.63786660976946e-07,\
     -7.58994648341957e-07,\
      -1.6736912030349e-08,\
      7.27128668886128e-07,\
     -6.18374425885138e-07,\
     -1.69017150409477e-06,\
      9.77887375569532e-06,\
     -1.59652315605579e-05,\
      5.38684438031965e-06,\
      3.67620321176578e-07,\
      4.78758994642156e-05,\
      -3.2972908401531e-05,\
     -0.000331756248215612,\
       0.00460642052084195,\
       -0.0147816884991402,\
        0.0549379983658079,\
       -0.0888173143699436,\
       -0.0403191959243558,\
         0.292953173664862,\
         -8.09710635241062,\
         -10.3650260389417,\
         -129.844564976801,\
         -387.005562405731,\
         -1675.24995538471,\
         -10203.9881306778,\
         -11240.4525144657,\
         -126107.639115735,\
          198983.418171164,\
          781719.532472805,\
          222578.562729136,\
          46638.5203763358,\
          1707.47248365029,\
          3001.55717298982,\
          646.488367541835,\
          165.639873325709,\
          43.1368816889835,\
          2.81512127691239,\
          1.69415084293283,\
        -0.162630603669399,\
        0.0343922501760851,\
       -0.0261906337667682,\
      -0.00266574197723834,\
       0.00188609628031843,\
      -0.00200584075887816,\
      0.000733170159693797,\
     -0.000149147632548399,\
     -1.55906996385135e-05,\
      2.28610960282525e-05,\
     -8.46053765284227e-06,\
      1.98402801173257e-06,\
      3.77587878104481e-08,\
      -6.4355213560915e-08,\
      2.69269975824045e-07,\
      9.26918586162506e-08,\
     -1.81745676946527e-08,\
      1.34703995090144e-07,\
      2.05521548209464e-09,\
      4.16624973546163e-08,\
       1.1287830153699e-07,\
     -1.54789615816159e-07,\
     -2.54999929010141e-07,\
       4.9576930235544e-07,\
     -2.14471993330123e-07],relativeTolerance,absoluteTolerance,requireSameLength=False)


del data
f.close()
print("numFailures:",numFailures)
exit(numFailures > 0)
