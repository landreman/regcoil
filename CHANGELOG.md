Changelog
=========

v0.1.3
------

* Much faster uniform-offset surfaces (especially for the default
  standard_toroidal_angle=False).
* Updated README to include saving and plotting.
* Added `scale()` function to scale surfaces in size and magnetic field strength.

v0.1.2
------

* plasma surface can be initialized from a simsopt vmec object by @landreman in #38
* Add major_radius, minor_radius, and aspect_ratio of surfaces by @landreman in #39
* Add codspeed benchmarks in #42 #43 #45
* prefer openblas to reference blas by @landreman in #44
* Improvements in memory allocation at initialisation and speedup of matrix assembly by @viper2642 in #40

v0.1.1
------
* Remove mentions of ADRs and phases in docstrings. by @landreman in #31
* Allow python version 3.10 by @landreman in #32
* Add max_Bnormal_over_B and avg_Bnormal_over_B by @landreman in #33
* Interface with simsopt virtual casing format by @landreman in #34
* fixed missing face on finite-build coil plot by @landreman in #35
* Better test coverage of unequal theta vs zeta resolutions by @landreman in #36
* More accurate volume by @landreman in #37

v0.1.0
------

Initial release of the python package (with compiled extension) on pypi
