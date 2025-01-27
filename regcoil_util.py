#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Main routine
if __name__=="__main__":
	import sys
	from regcoil import REGCOIL
	data=REGCOIL()
	data.read_regcoil('regcoil_out.GIGA_v503_100.nc')
	#data.plot_chisq()
	#data.plot_current_potential(nlambda=0)
	#data.plot_total_current_potential(nlambda=0)
	#data.plot_current_density(nlambda=0)
	#data.plot_bnormal_plasma()
	#data.plot_bnormal_coil()
	#data.plot_bnormal_total()
	sys.exit(0)
