#from distutils.core import setup
from setuptools import setup, find_packages

setup(name='pyREGCOIL',
	version = '1.0.0',
	description = 'Python library for interfacing with REGCOIL',
	long_description =	'This software package contains python '+ \
						'software for interfacing with the REGCOIL'+\
						'inputs and outputs.',
	author = 'Samuel A. Lazerson',
	author_email = 'lazersos@gmail.com',
	url = 'https://github.com/landreman/regcoil',
	packages=['pyREGCOIL'],
	scripts = ['regcoil_util.py'],
	install_requires=['numpy','matplotlib','netCDF4','contourpy']
	)
