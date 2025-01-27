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
		# Save filename
		self.filename = filename
		return

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
			Lambda index to plot (default: 0)
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
			Lambda index to plot (default: 0)
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
			Lambda index to plot (default: 0)
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
			Lambda index to plot (default: 0)
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

	def cut_coils(self,nlambda=0,numcoils=5,lplot=False):
		"""Cut coils from the REGCOIL potential

		This routine cuts coils from the REGCOIL potential.
		It allows the user to specify the number of coils per half 
		period.

		Parameters
		----------
		nlambda : int
			Lambda index to cut coils from (default: 0)
		numcoils : integer
			Number of coils per half period (default: 5)
		lplot : boolean (optional)
			Plot the potential and potential lines. (default: False)
		"""
		import numpy as np
		from contourpy import contour_generator, LineType
		import matplotlib.pyplot as plt
		# Compute total current
		Ipol = self.net_poloidal_current_Amperes
		# Setup x/y axis
		theta = self.theta_coil
		zeta  = self.zeta_coil # This is over full half field period
		nu = theta.shape[0]
		nv = zeta.shape[0]
		# Setup Data
		pot = np.transpose(self.current_potential[nlambda,:,:])
		# Compute Contour values
		cont_vals = np.zeros((numcoils))
		for k in range(numcoils):
			#u = 1
			u = round(nu*0.5)
			v = round((k+0.5)*nv/(numcoils*2))
			cont_vals[k] = pot[u,v]
		# Now calculate a larger potential map so coils can span periods
		nu2 = round(nu/(2*numcoils))
		nv2 = round(nv/(2*numcoils))
		theta_extend = np.concatenate((theta[nu-nu2:-1]-2*np.pi,theta,theta[0:nu2]+2*np.pi))
		zeta_extend = np.concatenate((zeta[nv-nv2:-1]-zeta[-1],zeta,zeta[0:nv2]+zeta[-1]))
		x = self.net_poloidal_current_Amperes/self.nfp
		pot_extend = np.concatenate((pot[nv-nv2:-1,:],pot,pot[0:nv2,:]),axis=0)
		pot_extend = np.concatenate((pot_extend[:,nu-nu2:-1]-x,pot_extend,pot_extend[:,0:nu2]+x),axis=1)
		# Now generate contours
		cont_gen = contour_generator(x=np.squeeze(zeta_extend),y=np.squeeze(theta_extend),z=np.squeeze(pot_extend), line_type=LineType.Separate)
		# Make plot if requested
		if lplot:
			px = 1/plt.rcParams['figure.dpi']
			fig=plt.figure(figsize=(1024*px,768*px))
			ax=fig.add_subplot(111)
			hmesh=ax.contourf(np.squeeze(zeta_extend),np.squeeze(theta_extend),np.squeeze(pot_extend),np.sort(cont_vals),extend='both',cmap='Greens')
			ax.contour(np.squeeze(zeta_extend),np.squeeze(theta_extend),np.squeeze(pot_extend),np.sort(cont_vals),colors='black')
			ax.set_xlabel('Toroidal angle [rad]')
			ax.set_ylabel('Poloidal angle [rad]')
			ax.set_title(r'REGCOIL Coil Cutting')
			plt.colorbar(hmesh,label=r'Potential $\Phi$ [arb]',ax=ax)
			plt.show()
		# Open file
		out_filename = self.filename.replace('regcoil_out.','coils.').replace('.nc','')
		f = open(out_filename+rf'_{nlambda:03d}','w')
		f.write(f"periods {self.nfp}\n")
		f.write(f"begin filament\n")
		f.write(f"mirror NIL\n")
		# Now loop over contours
		for k in range(numcoils):
			coil_name=f'MOD{k+1}'
			level = cont_gen.lines(cont_vals[k])
			th = np.array([]); ph = np.array([])
			for temp in level:
				th = np.append(th,temp[:,1])
				ph = np.append(ph,temp[:,0])
			# Fourier transform the coil
			npts = len(th)
			r = np.zeros((npts)); z = np.zeros((npts))
			for mn in range(self.mnmax_coil):
				mtheta = th*self.xm_coil[mn]
				nphi  = ph*self.xn_coil[mn]
				r  = r + np.cos(mtheta+nphi)*self.rmnc_coil[mn]
				z  = z + np.sin(mtheta+nphi)*self.zmns_coil[mn]
			# Check and adjust coil convention
			if (z[1]-z[0]) > 0:
				r = r[::-1]
				z = z[::-1]
			# Convert to XYZ and make current/group
			x = r * np.cos(ph)
			y = r * np.sin(ph)
			c = np.full((npts),Ipol/(self.nfp*numcoils*2))
			g = np.full((npts),k+1,dtype=int)
			c[-1] = 0.0
			# Create stellarator symmetric coil
			phn = 2.0*np.pi/self.nfp - ph
			xo = np.append(x,r[::-1]*np.cos(phn[::-1]))
			yo = np.append(y,r[::-1]*np.sin(phn[::-1]))
			zo = np.append(z,-z[::-1])
			co = np.append(c,c)
			go = np.append(g,g)
			x  = xo; y = yo; z = zo; c = co; g =go
			# Now make all field periods
			for mn in range(1,self.nfp):
				cop = np.cos(2*mn*np.pi/self.nfp)
				sip = np.sin(2*mn*np.pi/self.nfp)
				x = np.append(x,xo*cop - yo*sip)
				y = np.append(y,xo*sip + yo*cop)
				z = np.append(z,zo)
				c = np.append(c,co)
				g = np.append(g,go)
			# Write
			for j in range(len(x)):
				if c[j] == 0.0:
					f.write(f"{x[j]:.10E} {y[j]:.10E} {z[j]:.10E} {c[j]:.10E} {g[j]:02d} {coil_name}\n")
				else:
					f.write(f"{x[j]:.10E} {y[j]:.10E} {z[j]:.10E} {c[j]:.10E}\n")
		f.write(f"end\n")	
		f.close()





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