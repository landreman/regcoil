#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for handling and plotting REGCOIL data.
"""

# Libraries

# Constants

# RECOIL Class
class REGCOIL():
	"""Class for working with REGCOIL data
	"""
	def __init__(self):
		super().__init__()

	def read_regcoil(self,filename):
		"""Reads the REGCOIL netCDF file.

		This routine reads the REGCOIL netCDF regcoil_out file.

		Parameters
		----------
		file : str
			Path to regcoil_out file.
		"""
		import netCDF4
		import numpy as np
		nc = netCDF4.Dataset(filename,'r')
		for item in nc.variables:
			if item == 'lambda' or item == 'alpha':
				setattr(self,'lambdas',nc[item][()])
			elif item == 'nalpha':
				setattr(self,'nlambda',nc[item][()])
			elif item == 'chi2_J':
				setattr(self,'chi2_K',nc[item][()])
			elif item == 'J2':
				setattr(self,'K2',nc[item][()])
			else:
				setattr(self,item,nc[item][()])
		return
		# Sort in order of lambda
		permutation = np.argsort(self.lambdas)
		self.lambdas = self.lambdas[permutation]
		self.chi2_K = self.chi2_K[permutation]
		self.chi2_B = self.chi2_B[permutation]
		self.Bnormal_total = self.Bnormal_total[permutation]
		self.single_valued_current_potential_thetazeta = self.single_valued_current_potential_thetazeta[permutation]
		self.K2 = self.K2[permutation]
		self.current_potential = self.current_potential[permutation]
		if self.lambdas[-1] > 1.0E199: self.lambdas = np.inf

	def plot_chisq(self):
		"""Plots the REGCOIL chi^2

		This routine creates a plot showing the value of chi^2.

		Returns
		----------
		fig : Matplotlib figure object
			Figure object for plot.
		"""
		import matplotlib.pyplot as plt
		px = 1/plt.rcParams['figure.dpi']
		fig,ax = plt.subplots(2,3,figsize=(1024*px,768*px))
		ax[0,0].loglog(self.chi2_K,self.chi2_B,'.-r')
		ax[0,0].set_xlabel(rf'$\chi^2_K$ [$A^2$]')
		ax[0,0].set_ylabel(rf'$\chi^2_B$ [$T^2A^2$]')
		for j in range(self.nlambda):
			ax[0,0].plot(self.chi2_K[j], self.chi2_B[j],'ob')
		if self.nlambda > 1:
			ax[0,1].loglog(self.lambdas,self.chi2_B,'.-r')
			ax[0,1].set_xlabel(rf'$\lambda$ [$T^2m^2/A^2$]')
			ax[0,1].set_ylabel(rf'$\chi^2_B$ [$T^2A^2$]')
			for j in range(self.nlambda):
				ax[0,1].plot(self.lambdas[j], self.chi2_B[j],'ob')
		ax[0,2].semilogy(self.lambdas,self.chi2_B,'.-r')
		ax[0,2].set_xlabel(rf'$\lambda$ [$T^2m^2/A^2$]')
		ax[0,2].set_ylabel(rf'$\chi^2_B$ [$T^2A^2$]')
		for j in range(self.nlambda):
			ax[0,2].plot(self.lambdas[j], self.chi2_B[j],'ob')
		if self.nlambda > 1:
			ax[1,1].loglog(self.lambdas,self.chi2_K,'.-r')
			ax[1,1].set_xlabel(rf'$\lambda$ [$T^2m^2/A^2$]')
			ax[1,1].set_ylabel(rf'$\chi^2_K$ [$A^2$]')
			for j in range(self.nlambda):
				ax[1,1].plot(self.lambdas[j], self.chi2_K[j],'ob')
		ax[1,2].semilogy(self.lambdas,self.chi2_K,'.-r')
		ax[1,2].set_xlabel(rf'$\lambda$ [$T^2m^2/A^2$]')
		ax[1,2].set_ylabel(rf'$\chi^2_K$ [$A^2$]')
		for j in range(self.nlambda):
			ax[1,2].plot(self.lambdas[j], self.chi2_K[j],'ob')
		plt.tight_layout()
		plt.show()
		return fig

	def plot_current_potential(self,nlambda=0,numCoutours=20,ax=None):
		"""Plots the single valued current potential

		This routine plots the single valued current potential.

		Parameters
		----------
		nlambda : int
			Lambda value to plot (default: 0)
		numCoutours: int
			Number of contours to plot (default: 20)
		ax : axes (optional)
			Matplotlib axes object to plot to.
		"""
		import numpy as np
		import matplotlib.pyplot as plt
		# Handle passes ax
		lplotnow = False
		if not ax:
			ax = plt.axes()
			lplotnow = True
		# Plot potential
		hmesh=ax.contourf(self.zeta_coil,self.theta_coil,np.transpose(self.single_valued_current_potential_thetazeta[nlambda,:,:]),numCoutours)
		plt.colorbar(hmesh,label=rf'$\Phi$ [arb]',ax=ax)
		ax.set_xlabel(rf'Toroidal Angle ($\zeta$)')
		ax.set_ylabel(rf'Poloidal Angle ($\theta$)')
		ax.set_title(rf'Current Potential ($\lambda$ {self.lambdas[nlambda]:7.2e})')
		# Show if standalone
		if lplotnow: plt.show()

	def plot_total_current_potential(self,nlambda=0,numCoutours=5,ax=None):
		"""Plots the total single valued current potential

		This routine plots the total single valued current potential.

		Parameters
		----------
		nlambda : int
			Lambda value to plot (default: 0)
		numCoutours: int
			Number of contours to plot (default: 5)
		ax : axes (optional)
			Matplotlib axes object to plot to.
		"""
		import numpy as np
		import matplotlib.pyplot as plt
		# Handle passes ax
		lplotnow = False
		if not ax:
			ax = plt.axes()
			lplotnow = True
		# Setup Contours
		x = self.net_poloidal_current_Amperes/self.nfp
		if abs(x) > np.finfo(float).eps:
			contours = np.linspace(-0.5*x,1.5*x,numCoutours*4,endpoint=False)
		else:
			x = np.max(self.current_potential)
			contours = np.linspace(-0.5*x,1.5*x,numCoutours*4,endpoint=False)
		d = contours[1]-contours[0]
		contours = contours + d*0.5
		contours.sort()
		# Setup x/y axis
		zeta_coil_extended = np.concatenate((self.zeta_coil,[2*np.pi/self.nfp]))
		theta_coil_extended = np.concatenate((self.theta_coil,[2*np.pi]))
		# Setup Data
		data = np.transpose(self.current_potential[nlambda,:,:])
		data_extended = np.concatenate((data,data[0:1,:]),axis=0)
		data_extended = np.concatenate((data_extended,data_extended[:,0:1]+x),axis=1)
		# Plot potential
		hmesh=ax.contourf(zeta_coil_extended, theta_coil_extended, data_extended, contours)
		plt.colorbar(hmesh,label=rf'Total $\Phi$ [arb]',ax=ax)
		ax.set_xlabel(rf'Toroidal Angle ($\zeta$)')
		ax.set_ylabel(rf'Poloidal Angle ($\theta$)')
		ax.set_title(rf'Total Current Potential ($\lambda$ {self.lambdas[nlambda]:7.2e})')
		# Show if standalone
		if lplotnow: plt.show()

	def plot_current_density(self,nlambda=0,numCoutours=20,ax=None):
		"""Plots the current density

		This routine plots the current density.

		Parameters
		----------
		nlambda : int
			Lambda value to plot (default: 0)
		numCoutours: int
			Number of contours to plot (default: 20)
		ax : axes (optional)
			Matplotlib axes object to plot to.
		"""
		import numpy as np
		import matplotlib.pyplot as plt
		# Handle passes ax
		lplotnow = False
		if not ax:
			ax = plt.axes()
			lplotnow = True
		# Plot potential
		hmesh=ax.contourf(self.zeta_coil,self.theta_coil,np.transpose(self.K2[nlambda,:,:]),numCoutours)
		plt.colorbar(hmesh,label=rf'K [A/m]',ax=ax)
		ax.set_xlabel(rf'Toroidal Angle ($\zeta$)')
		ax.set_ylabel(rf'Poloidal Angle ($\theta$)')
		ax.set_title(rf'Current Density ($\lambda$ {self.lambdas[nlambda]:7.2e})')
		# Show if standalone
		if lplotnow: plt.show()

	def plot_bnormal_plasma(self,numCoutours=20,ax=None):
		"""Plots the plasma Bnormal

		This routine plots the normal magnetic field on the plasma
		surface generated by the plasma.

		Parameters
		----------
		numCoutours: int
			Number of contours to plot (default: 20)
		ax : axes (optional)
			Matplotlib axes object to plot to.
		"""
		import numpy as np
		import matplotlib.pyplot as plt
		# Handle passes ax
		lplotnow = False
		if not ax:
			ax = plt.axes()
			lplotnow = True
		# Plot potential
		hmesh=ax.contourf(self.zeta_plasma,self.theta_plasma,np.transpose(self.Bnormal_from_plasma_current),numCoutours)
		plt.colorbar(hmesh,label=r'$\vec{B}\cdot\hat{n}$ [T]',ax=ax)
		ax.set_xlabel(rf'Toroidal Angle ($\zeta$)')
		ax.set_ylabel(rf'Poloidal Angle ($\theta$)')
		ax.set_title(rf'Plasma Normal Field')
		# Show if standalone
		if lplotnow: plt.show()

	def plot_bnormal_coil(self,numCoutours=20,ax=None):
		"""Plots the coil Bnormal

		This routine plots the normal magnetic field on the plasma
		surface generated by the coil.

		Parameters
		----------
		numCoutours: int
			Number of contours to plot (default: 20)
		ax : axes (optional)
			Matplotlib axes object to plot to.
		"""
		import numpy as np
		import matplotlib.pyplot as plt
		# Handle passes ax
		lplotnow = False
		if not ax:
			ax = plt.axes()
			lplotnow = True
		# Plot potential
		hmesh=ax.contourf(self.zeta_plasma,self.theta_plasma,np.transpose(self.Bnormal_from_net_coil_currents),numCoutours)
		plt.colorbar(hmesh,label=r'$\vec{B}\cdot\hat{n}$ [T]',ax=ax)
		ax.set_xlabel(rf'Toroidal Angle ($\zeta$)')
		ax.set_ylabel(rf'Poloidal Angle ($\theta$)')
		ax.set_title(rf'Net Coil Normal Field')
		# Show if standalone
		if lplotnow: plt.show()

	def plot_bnormal_total(self,nlambda=0,numCoutours=20,ax=None):
		"""Plots the current density

		This routine plots the current density.

		Parameters
		----------
		nlambda : int
			Lambda value to plot (default: 0)
		numCoutours: int
			Number of contours to plot (default: 20)
		ax : axes (optional)
			Matplotlib axes object to plot to.
		"""
		import numpy as np
		import matplotlib.pyplot as plt
		# Handle passes ax
		lplotnow = False
		if not ax:
			ax = plt.axes()
			lplotnow = True
		# Plot potential
		hmesh=ax.contourf(self.zeta_plasma,self.theta_plasma,np.transpose(self.Bnormal_total[nlambda,:,:]),numCoutours)
		plt.colorbar(hmesh,label=r'$\vec{B}\cdot\hat{n}$ [T]',ax=ax)
		ax.set_xlabel(rf'Toroidal Angle ($\zeta$)')
		ax.set_ylabel(rf'Poloidal Angle ($\theta$)')
		ax.set_title(rf'Total Normal Field ($\lambda$ {self.lambdas[nlambda]:7.2e})')
		# Show if standalone
		if lplotnow: plt.show()




# RECOIL_INPUT Class
class REGCOIL_INPUT():
	"""Class for working with REGCOIL namelist input
	"""
	def __init__(self):
		super().__init__()

	def read_regcoil_nml(self,filename):
		"""Reads the REGCOIL namelist from a file.

		This routine reads and initializes the REGCOIL
		class with information from the REGCOIL_NML
		Fortran Namelist

		Parameters
		----------
		file : str
			Path to regcoil_in file.

		"""

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