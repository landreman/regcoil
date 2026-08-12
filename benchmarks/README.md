The tests in this directory are for timing by codtest.io.
They are run in the codspeed CI but not run as part of the standard testing CI.
To run them locally, you need the `pytest-codspeed` package installed via pip.
You can run the benchmarks from the project root directory with
`pytest benchmarks/ --codtest`