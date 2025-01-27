#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Main routine
if __name__=="__main__":
	import sys
	from pyREGCOIL.regcoil import REGCOIL
	from argparse import ArgumentParser
	parser = ArgumentParser(description= 
		'''Provides a utility for working with REGCOIL outputs.''')
	parser.add_argument("-f", "--file", dest="file_name",
		help="REGCOIL output filename.", default = None)
	parser.add_argument("-p", "--plot", dest="lplots", action='store_true',
		help="Produce some plots.", default = False)
	args = parser.parse_args()
	data=REGCOIL()
	data.read_regcoil(args.file_name)
	if args.lplots:
		data.plot_chisq()
		for n in range(data.nlambda):
			data.plot_current_potential(nlambda=n)
		for n in range(data.nlambda):
			data.plot_total_current_potential(nlambda=n)
		for n in range(data.nlambda):
			data.plot_current_density(nlambda=n)
		data.plot_bnormal_plasma()
		data.plot_bnormal_coil()
		for n in range(data.nlambda):
			data.plot_bnormal_total(nlambda=n)
	sys.exit(0)
