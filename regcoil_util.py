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
	parser.add_argument("-c", "--cutcoil", dest="lcutcoil", action='store_true',
		help="Cut coils from potential.", default = False)
	parser.add_argument("-l", "--lambda", dest="lam_value", type=int,
		help="Select a specific lambda value.", default = 0)
	args = parser.parse_args()
	data=REGCOIL()
	data.read_regcoil(args.file_name)
	lam_value = args.lam_value
	if (lam_value >= data.nlambda) or (lam_value < 0):
		print(rf"Requested lambda index: {lam_value:2d} greater than maximum lambda {data.nlambda-1:2d}")
		sys.exit(-1)
	if args.lplots and not args.lcutcoil:
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
	if args.lcutcoil:
		data.cut_coils(nlambda=lam_value,lplot=args.lplots)
	sys.exit(0)
