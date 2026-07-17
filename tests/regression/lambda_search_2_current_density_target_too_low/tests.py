#!/usr/bin/env python

# This python script checks the output file for an example to 
# see if the results are close to expected values.  This script may be
# run directly. This is a copy under tests/regression/ for pytest;
# the original lives under examples/ and is used by "make test".

from pathlib import Path as _Path
exec((_Path(__file__).resolve().parents[1] / "testsCommon.py").read_text())
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


del data
f.close()
print("numFailures:",numFailures)
exit(numFailures > 0)
