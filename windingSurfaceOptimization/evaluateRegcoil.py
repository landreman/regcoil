import numpy as np
import os
import re
from regcoilScan import readVariable, namelistLineContains
from shutil import copyfile
import subprocess
import sys
from scipy.io import netcdf

class coilFourier:

  def __init__(self,nmax,mmax,regcoil_input_file):
    self.regcoil_input_file = regcoil_input_file
	# Load in variables from input file
    nmax_sensitivity = readVariable("nmax_sensitivity","int",regcoil_input_file,required=True)
    mmax_sensitivity = readVariable("mmax_sensitivity","int",regcoil_input_file,required=True)
    geometry_option_coil = readVariable("geometry_option_coil","int",regcoil_input_file,required=False)
    geometry_option_plasma = readVariable("geometry_option_plasma","int",regcoil_input_file,required=False)
    load_bnorm = readVariable("load_bnorm","string",regcoil_input_file,required=False)
    general_option = readVariable("general_option","int",regcoil_input_file,required=False)
    fixed_norm_sensitivity_option = readVariable("fixed_norm_sensitivity_option","string",regcoil_input_file,required=False)
    sensitivity_option = readVariable("sensitivity_option","int",regcoil_input_file,required=False)

  # Set to default values if absent
    if (geometry_option_plasma is None):
		  geometry_option_plasma = 0
    if (geometry_option_coil is None):
		  geometry_option_coil = 0
    if (load_bnorm is None):
			load_bnorm = False
    if (load_bnorm == ".true." or load_bnorm == "T" or load_bnorm == ".TRUE."):
			load_bnorm = True
    if (general_option is None):
			general_option = 1
    if (fixed_norm_sensitivity_option is None):
			fixed_norm_sensitivity_option = False
    if (fixed_norm_sensitivity_option == ".true." or fixed_norm_sensitivity_option == "T" or fixed_norm_sensitivity_option == ".TRUE."):
			fixed_norm_sensitivity_option = True
    if (sensitivity_option is None):
			sensitivity_option = 1
	
	  # Set class variables
    self.geometry_option_coil = geometry_option_coil
    self.geometry_option_plasma = geometry_option_plasma
    self.load_bnorm = load_bnorm
    self.general_option = general_option
    self.fixed_norm_sensitivity_option = fixed_norm_sensitivity_option
    self.sensitivity_option = sensitivity_option
	
    if (general_option > 3):
      target_value = readVariable("target_value","float",regcoil_input_file,required=False)
      self.target_value = target_value
      self.target_value_init = target_value
      target_option = readVariable("target_option","string",regcoil_input_file,required=True)
      self.target_option = target_option.strip()[1:-1]
	
    spectral_norm_p = readVariable("spectral_norm_p","float",regcoil_input_file,required=False)
    spectral_norm_q = readVariable("spectral_norm_q","float",regcoil_input_file,required=False)
    self.spectral_norm_p = spectral_norm_p
    self.spectral_norm_q = spectral_norm_q
    if (self.spectral_norm_p==None):
      self.spectral_norm_p = 2
    if (self.spectral_norm_q==None):
      self.spectral_norm_q = 2

    alphaV = readVariable("alphaV","float",regcoil_input_file,required=False)
    alphaS = readVariable("alphaS","float",regcoil_input_file,required=False)
    alphaD = readVariable("alphaD","float",regcoil_input_file,required=False)
    alphaB = readVariable("alphaB","float",regcoil_input_file,required=False)
    alphaK = readVariable("alphaK","float",regcoil_input_file,required=False)
    alphaD_tanh = readVariable("alphaD_tanh","float",regcoil_input_file,required=False)
    alphaMaxK = readVariable("alphaMaxK","float",regcoil_input_file,required=False)

    scaleFactor = readVariable("scaleFactor","float",regcoil_input_file,required=False)
    d_min_target = readVariable("d_min_target","float",regcoil_input_file,required=False)
    alphaD_tanh_scale = readVariable("alphaD_tanh_scale","float",regcoil_input_file,required=False)
    
    self.alphaV = alphaV
    self.alphaS = alphaS
    self.alphaD = alphaD
    self.alphaB = alphaB
    self.alphaK = alphaK
    self.alphaMaxK = alphaMaxK
    self.scaleFactor = scaleFactor
    self.alphaD_tanh = alphaD_tanh
    self.d_min_target = d_min_target
    self.alphaD_tanh_scale = alphaD_tanh_scale
    if (alphaV == None):
      self.alphaV = 0
    if (alphaS == None):
      self.alphaS = 0
    if (alphaD == None):
      self.alphaD = 0
    if (alphaB == None):
      self.alphaB = 0
    if (alphaK == None):
      self.alphaK = 0
    if (alphaMaxK == None):
      self.alphaMaxK = 0
    if (scaleFactor == None):
      self.scaleFactor = 1.0
    if (alphaD_tanh == None):
      self.alphaD_tanh = 0
    if (alphaD_tanh_scale == None):
      self.alphaD_tanh_scale = 1.0
    if (d_min_target == None):
      self.d_min_target = 0.1
		
    if (self.alphaV == 0 and self.alphaS == 0 and self.alphaD == 0 and self.alphaB == 0 and self.alphaK == 0 and self.alphaMaxK == 0 and self.alphaD_tanh == 0):
      print("Error! One of the alpha parameters must be non-zero.")
      sys.exit(0)

		# Check for correct input
    if (geometry_option_coil != 3 and geometry_option_coil != 4):
      print("Error! This script is only compatible with geometry_option_coil=3 or 4 at the) moment."
      sys.exit(0)

    if (abs(self.alphaMaxK)>0 and (target_option=='max_K_lse') and fixed_norm_sensitivity_option):
			print("Error! K_max included in objective function but K_max is held fixed in gradient calculation.")
			sys.exit(0)

    if (abs(self.alphaB)>0 and (target_option=='chi2_B') and fixed_norm_sensitivity_option):
			print("Error! chi2_B included in objective function but chi2_B is held fixed in gradient calculation.")
			sys.exit(0)

    if ((abs(self.alphaB)>0 or abs(self.alphaK)>0 or abs(self.alphaMaxK)>0) and self.sensitivity_option < 2):
			print("Warning! chi2_B, K_rms, or K_max are included in the objective function but their derivatives are not computed in REGCOIL. Make sure that grad_option = 1 (using scipy_optmize) so finite difference derivatives are performed.")

    self.nmax_sensitivity = nmax_sensitivity
    self.mmax_sensitivity = mmax_sensitivity
    self.nmax = nmax
    self.mmax = mmax
    # Number of modes - does not include factor of 2 from (rmnc,zmns) - not # of Fourier coefficients
    nmodes_sensitivity = (nmax_sensitivity+1) + (2*nmax_sensitivity+1)*mmax_sensitivity
    self.dspectral_normdomegas = np.zeros(2*nmodes_sensitivity)
    self.spectral_norm = 0
    nmodes = (nmax+1) + (2*nmax+1)*mmax
    self.nmodes = nmodes
    self.nmodes_sensitivity = nmodes_sensitivity
    self.rmncs = np.zeros(nmodes)
    self.zmnss = np.zeros(nmodes)
    self.rmnss = np.zeros(nmodes)
    self.zmncs = np.zeros(nmodes)
    self.xn = np.zeros(nmodes)
    self.xm = np.zeros(nmodes)
    self.xn_sensitivity = np.zeros(nmodes_sensitivity)
    self.xm_sensitivity = np.zeros(nmodes_sensitivity)
    self.omegas = np.zeros(2*nmodes)
    self.omegas_sensitivity = np.zeros(2*nmodes_sensitivity)
    self.objective_function = 0
    self.dobjective_functiondomegas_sensitivity = np.zeros(2*nmodes_sensitivity)
    self.evaluated = False # has function been evaluated at current omegas_sensitivity?
    self.feval = 0 # Number of function evaluations
    self.chi2B = 0
    self.chi2K = 0
    self.rms_K = 0
    self.area_coil = 0
    self.dchi2Bdomega = 0
    self.dchi2Kdomega = 0
    self.darea_coildomega = 0
    self.coil_volume = 0
    self.dcoil_volumedomega = 0
    self.dcoil_plasma_dist_mindomega = 0
    self.drms_Kdomega = 0
    self.coil_plasma_dist_min_lse = 0
    self.coil_plasma_dist_max_lse = 0
    self.increased_target_current = False
    self.decreased_target_current = False
    self.current_factor = 0.1
    # Initialize m = 0 modes
    imode = 0
    for ni in range(0,nmax+1):
      self.xn[imode] = ni
      self.xm[imode] = 0
      imode = imode + 1
    # Initialize m > 0 modes
    for mi in range(1,mmax+1):
      for ni in range(-nmax,nmax+1):
        self.xn[imode] = ni
        self.xm[imode] = mi
        imode = imode + 1
    imode = 0
    for ni in range(0,nmax_sensitivity+1):
      self.xn_sensitivity[imode] = ni
      self.xm_sensitivity[imode] = 0
      imode = imode + 1
    for mi in range(1,mmax_sensitivity+1):
      for ni in range(-nmax_sensitivity,nmax_sensitivity+1):
        self.xn_sensitivity[imode] = ni
        self.xm_sensitivity[imode] = mi
        imode = imode + 1

	  # Open output file
    filename = "windingSurfaceOptimization.txt"
    file = open(filename,'w')
    self.file = file
    file.write('eval \t')
    file.write('objective_function \t norm(dobjective_functiondomega) \t')
    if (abs(self.alphaK) > 0):
      file.write('rms_K \t norm(drms_Kdomega) \t')
    if (abs(self.alphaB) > 0):
      file.write('chi2B \t norm(dchi2Bdomega) \t')
    if (abs(self.alphaV) > 0):
      file.write('coil_volume \t norm(dcoil_volumedomega) \t')
    if (abs(self.alphaD) or abs(self.alphaD_tanh)>0):
      file.write('coil_plasma_dist_min \t norm(dcoil_plasma_dist_mindomega) \t')
    if (abs(self.alphaS)):
      file.write('spectral_norm \t norm(dspectral_normdomegas) \t')

  def set_omegas_sensitivity(self,new_omegas_sensitivity):
    self.omegas_sensitivity = new_omegas_sensitivity
    for imn in range(0, self.nmodes_sensitivity):
      for jmn in range(0, self.nmodes):
        if (self.xn[jmn] == self.xn_sensitivity[imn]):
          if (self.xm[jmn] == self.xm_sensitivity[imn]):
            self.rmncs[jmn] = self.omegas_sensitivity[2*imn]
            self.zmnss[jmn] = self.omegas_sensitivity[2*imn+1]
            self.omegas[2*jmn] = self.omegas_sensitivity[2*imn]
            self.omegas[2*jmn+1] = self.omegas_sensitivity[2*imn+1]
  # Function has not been evaluated at new omegas
    self.evaluated = False
  
  def set_drms_Kdomega(self,drms_Kdomega):
    self.drms_Kdomega = drms_Kdomega

  def set_dchi2Bdomega(self,dchi2Bdomega):
    self.dchi2Bdomega = dchi2Bdomega
  
  def set_dchi2Kdomega(self,dchi2Kdomega):
    self.dchi2Kdomega = dchi2Kdomega
  
  def set_darea_coildomega(self,darea_coildomega):
    self.darea_coildomega = darea_coildomega

  def set_dcoil_volumedomega(self,dcoil_volumedomega):
    self.dcoil_volumedomega = dcoil_volumedomega

  def set_dcoil_plasma_dist_mindomega(self,dcoil_plasma_dist_mindomega):
    self.dcoil_plasma_dist_mindomega = dcoil_plasma_dist_mindomega
  
  def set_dobjective_functiondomegas(self,dobjective_functiondomegas_sensitivity):
    self.dobjective_functiondomegas_sensitivity = dobjective_functiondomegas_sensitivity
  
  def increment_feval(self):
    self.feval = self.feval + 1
  
  def compute_spectral_norm(self):
    self.spectral_norm = 0
    dspectral_normdomegas = np.zeros(2*self.nmodes_sensitivity)
    for imode in range(0,self.nmodes_sensitivity):
      self.spectral_norm = self.spectral_norm + self.xm_sensitivity[imode]**(self.spectral_norm_p)*(self.omegas_sensitivity[2*imode]**2 + self.omegas_sensitivity[2*imode+1]**2)
      dspectral_normdomegas[2*imode] = self.xm_sensitivity[imode]**(self.spectral_norm_p)*2*self.omegas_sensitivity[2*imode]
      dspectral_normdomegas[2*imode+1] = self.xm_sensitivity[imode]**(self.spectral_norm_p)*2*self.omegas_sensitivity[2*imode+1]
    self.set_dspectral_normdomegas(dspectral_normdomegas)
  
  def set_Fourier_from_nescin(self,nescin_file):
    file = open(nescin_file, "r")
    imode = 0
    inCurrentGeometry = 0
    for line in file:
      if (inCurrentGeometry):
        list = line.split()
        this_m = int(list[0])
        this_n = int(list[1])
        rmnc = float(list[2])
        zmns = float(list[3])
        for imn in range(0,self.nmodes):
          if (self.xm[imn]==this_m and self.xn[imn]==this_n):
            self.rmncs[imn] = rmnc
            self.zmnss[imn] = zmns
            self.omegas[2*imn] = rmnc
            self.omegas[2*imn+1] = zmns
          elif (self.xm[imn]==0 and this_m == 0 and self.xn[imn]==-this_n):
						self.rmncs[imn] = rmnc
						self.zmnss[imn] = -zmns
						self.omegas[2*imn] = rmnc
						self.omegas[2*imn+1] = -zmns
        for imn in range(0,self.nmodes_sensitivity):
          if (self.xn_sensitivity[imn]==this_n and self.xm_sensitivity[imn]==this_m):
            self.omegas_sensitivity[2*imn] = rmnc
            self.omegas_sensitivity[2*imn+1] = zmns
          elif (self.xm_sensitivity[imn] == 0 and this_m == 0 and self.xn_sensitivity[imn]==-this_n):
						self.omegas_sensitivity[2*imn] = rmnc
						self.omegas_sensitivity[2*imn+1] = -zmns
      if re.match("------ Current Surface",line):
        inCurrentGeometry = 1
        next(file)
        next(file)
        next(file)
        next(file)
    file.close()

  def set_dspectral_normdomegas(self,new_dspectral_normdomegas):
    self.dspectral_normdomegas = new_dspectral_normdomegas
  
  def evaluateObjectiveFunction(self):
    self.compute_spectral_norm()
    dobjective_functiondomegas = 0
    self.objective_function = 0
    if (abs(self.alphaV)>0):
        self.objective_function = self.objective_function - self.alphaV*self.coil_volume**(1.0/3.0)
        dobjective_functiondomegas = self.objective_function - self.alphaV*(1.0/3.0)*(self.coil_volume**(-2.0/3.0))*self.dcoil_volumedomega
    if (abs(self.alphaB)>0):
        self.objective_function = self.objective_function + self.alphaB*self.chi2B
        dobjective_functiondomegas = dobjective_functiondomegas + self.alphaB*self.dchi2Bdomega
    if (abs(self.alphaS)>0):
				self.objective_function = self.objective_function + self.alphaS*self.spectral_norm
				dobjective_functiondomegas = dobjective_functiondomegas + self.alphaS*self.dspectral_normdomegas
    if (abs(self.alphaK)>0):
				self.objective_function = self.objective_function + self.alphaK*self.rms_K
				dobjective_functiondomegas = dobjective_functiondomegas + self.alphaK*self.drms_Kdomega
    if (abs(self.alphaD)>0):
				self.objective_function = self.objective_function - self.alphaD*self.coil_plasma_dist_min_lse
				dobjective_functiondomegas = dobjective_functiondomegas - self.alphaD*self.dcoil_plasma_dist_mindomega
    if (abs(self.alphaD_tanh)>0):
        self.objective_function = self.objective_function + self.alphaD_tanh*(1+np.tanh((self.d_min_target-self.coil_plasma_dist_min_lse)/self.alphaD_tanh_scale))
        dobjective_functiondomegas = dobjective_functiondomegas - (self.alphaD_tanh/self.alphaD_tanh_scale)*(np.cosh((self.d_min_target-self.coil_plasma_dist_min_lse)/self.alphaD_tanh_scale)**(-2.0))*self.dcoil_plasma_dist_mindomega

    self.objective_function = self.scaleFactor*self.objective_function
    self.set_dobjective_functiondomegas(self.scaleFactor*dobjective_functiondomegas)

		# Write output to file
    self.file.write('\n')
    self.file.write(str(self.feval)+'\t')
    self.file.write(str(self.objective_function)+'\t')
    self.file.write(str(np.linalg.norm(self.dobjective_functiondomegas_sensitivity,2))+'\t')
    if (abs(self.alphaK) > 0):
        self.file.write(str(self.rms_K)+'\t'+str(np.linalg.norm(self.drms_Kdomega,2))+'\t')
    if (abs(self.alphaB) > 0):
        self.file.write(str(self.chi2B)+'\t'+str(np.linalg.norm(self.dchi2Bdomega,2))+'\t')
    if (abs(self.alphaV) > 0):
        self.file.write(str(self.coil_volume)+'\t'+str(np.linalg.norm(self.dcoil_volumedomega,2))+'\t')
    if (abs(self.alphaD) > 0 or abs(self.alphaD_tanh > 0)):
        self.file.write(str(self.coil_plasma_dist_min_lse)+'\t'+str(np.linalg.norm(self.dcoil_plasma_dist_mindomega,2))+'\t')
    if (abs(self.alphaS) > 0):
        self.file.write(str(self.spectral_norm)+'\t'+str(np.linalg.norm(self.dspectral_normdomegas,2)))

  # This is a script to be called within a nonlinear optimization routine in order to evaluate
  # chi2 and its gradient with respect to the Fourier coefficients
  def evaluateRegcoil(self,omegas_sensitivity_new,target_value=0):
    
    self.set_omegas_sensitivity(omegas_sensitivity_new.copy())
    
    regcoil_input_file = self.regcoil_input_file
	
    if (self.geometry_option_plasma == 2 or self.geometry_option_plasma == 3 or self.geometry_option_plasma == 4):
    	wout_filename = readVariable("wout_filename","string",regcoil_input_file,required=True)
    elif (self.geometry_option_plasma == 5):
			efit_filename = readVariable("efit_filename","string",regcoil_input_file,required=True)
    elif (self.geometry_option_plasma == 6):
			shape_filename = readVariable("shape_filename_plasma","string",regcoil_input_file,required=True)
    nescin_filename = readVariable("nescin_filename","string",regcoil_input_file,required=True)
    if (self.load_bnorm):
    	bnorm_filename = readVariable("bnorm_filename","string",regcoil_input_file,required=True)
	
    # Create new directory
    directory = "eval_" + str(self.feval)
    if (not os.path.isdir(directory)):
      os.makedirs(directory)
    os.chdir(directory)

    # Copy nescin file
    src = "../" + nescin_filename
    dst = nescin_filename
    copyfile(src,dst)

    # Copy regcoil_in file
    src = "../" + regcoil_input_file
    dst = regcoil_input_file
    copyfile(src,dst)

    # Edit new nescin file
    new_nescin = nescin_filename + "_" + str(self.feval)
    self.create_nescin(nescin_filename,new_nescin)
    # Remove old nescin file
    os.remove(nescin_filename)

    # Edit regcoil_in file
    with open(regcoil_input_file, 'r') as f:
      inputFile = f.readlines()
    f = open(regcoil_input_file,"w")
    for line in inputFile:
        if namelistLineContains(line,"nescin_filename"):
            line = 'nescin_filename = "'+new_nescin+'"\n'
        if (self.geometry_option_plasma ==2 or self.geometry_option_plasma == 3 or self.geometry_option_plasma == 4):
				if namelistLineContains(line,"wout_filename"):
					if (wout_filename[0]=='/'):
						new_wout = wout_filename
					else:
						new_wout = '../' + wout_filename
					line = 'wout_filename = "'+new_wout+'"\n'
				elif (self.geometry_option_plasma == 5):
					if namelistLineContains(line,"efit_filename"):
						new_efit = '../' + efit_filename
				elif (self.geometry_option_plasma == 6):
					if namelistLineContains(line,"shape_filename_plasma"):
						new_shape_file = '../' + shape_filename
				if (self.load_bnorm):
					if namelistLineContains(line,"bnorm_filename"):
						new_bnorm = '../' + bnorm_filename
						line = 'bnorm_filename = "'+new_bnorm+'"\n'
				if (self.general_option > 3):
					if namelistLineContains(line,"target_value"):
						line = 'target_value = '+str(target_value)+'\n'
				f.write(line)
    f.close()

    submitCommand = "regcoil " + regcoil_input_file
    outputFileName = "regcoil_out" + regcoil_input_file[10::]
    g = open(outputFileName,"w")
    try:
      submissionResult = subprocess.call(submitCommand.split(" "),stdout=g)
    except:
      print("ERROR: Unable to submit run "+directory+" for some reason.")
      raise
    else:
      if submissionResult==0:
        print("No errors submitting job "+directory)
      else:
        print("Nonzero exit code returned when trying to submit job "+directory)
        exit
    g.close()

    # Obtain objective function and its derivative
    cdfFileName = outputFileName + ".nc"
    try:
      f = netcdf.netcdf_file(cdfFileName,'r',mmap=False)
    except:
      print("Unable to open "+cdfFileName+" even though this file exists.")
      sys.exit(0)
    try:
      dummy = f.variables["K2"][()]
    except:
      print("Unable to read "+cdfFileName+" even though this file exists.")
      sys.exit(0)

    exit_code = f.variables["exit_code"][()]
    if (exit_code == 0):
			self.area_coil = f.variables["area_coil"][()]
			self.chi2B = f.variables["chi2_B"][()][-1]
			self.chi2K = f.variables["chi2_K"][()][-1]
			self.rms_K = np.sqrt(self.chi2K/self.area_coil)
			self.coil_volume = f.variables["volume_coil"][()]
		
			if (self.sensitivity_option>1):
				self.set_darea_coildomega(f.variables["darea_coildomega"][()])
				self.set_dcoil_volumedomega(f.variables["dvolume_coildomega"][()])
				self.set_dcoil_plasma_dist_mindomega(f.variables["dcoil_plasma_dist_mindomega"][()])
				self.coil_plasma_dist_min_lse = f.variables["coil_plasma_dist_min_lse"][()]
				self.coil_plasma_dist_max_lse = f.variables["coil_plasma_dist_max_lse"][()]
			if (self.sensitivity_option>2):
				self.set_dchi2Bdomega(f.variables["dchi2Bdomega"][()][-1])
				self.set_dchi2Kdomega(f.variables["dchi2Kdomega"][()][-1])
				drms_Kdomega = 0.5*self.rms_K**(-1.0)*(self.dchi2Kdomega/self.area_coil - self.chi2K*self.darea_coildomega/self.area_coil**2)
			
			self.set_drms_Kdomega(drms_Kdomega)

			self.evaluateObjectiveFunction()
    
			os.chdir('..')
			self.increment_feval()
			self.evaluated = True
			if (self.general_option > 3):
				self.increased_target_current = False
				self.decreased_target_current = False
				self.current_factor = 0.1
				self.target_value = self.target_value_init
    else:
      print("Error! Job did not complete.")
      if (exit_code == -1): # did not converge in nlambda iterations
        nlambda = f.variables["nlambda"][()]
        new_nlambda = nlambda*2
        print("Trying again with nlambda = " + str(new_nlambda))
        # Edit input file with more nlambda
        os.chdir('..')
        with open(regcoil_input_file, 'r') as f:
          inputFile = f.readlines()
          f = open(regcoil_input_file,"w")
          for line in inputFile:
            if namelistLineContains(line,"nlambda"):
              line = 'nlambda = '+str(new_nlambda)+'\n'
            f.write(line)
          f.close()
        if (self.general_option > 3):
          self.evaluateRegcoil(omegas_sensitivity_new,self.target_value)
        else:
          self.evaluateRegcoil(omegas_sensitivity_new)
      # exit_code == -2 or -3 should only happen with general_option > 3
      elif (exit_code == -2): # current density too low or chi2B too high
        if (self.target_option == 'max_K_lse' or self.target_option == 'lp_norm_K'):
          print("Current density too low.")
          # Decrease factor of increase/decrease
          if (self.decreased_target_current): # previously tried decreasing target
            self.current_factor = self.current_factor*0.5
            print("current_factor is now: " + str(self.current_factor))
          self.target_value = (1.0+self.current_factor)*self.target_value
          print("Trying again with target_value = " + str(self.target_value))
          os.chdir('..')
          self.increased_target_current = True
          self.evaluateRegcoil(omegas_sensitivity_new,self.target_value)
        elif (self.target_option == 'chi2_B'):
          print("chi2B too high.")
          if (self.increased_target_current): # previously tried increasing
            self.current_factor = self.current_factor*0.5
            print("current_factor is now: " + str(self.current_factor))
          self.target_value = (1.0-self.current_factor)*self.target_value
          print "Trying again with target_value = " + str(self.target_value)
          os.chdir('..')
          self.decreased_target_current = True
          self.evaluateRegcoil(omegas_sensitivity_new,self.target_value)
        else:
					print("Error! Incorrect target_option: "+str(self.target_option))
					sys.exit(0)
      elif (exit_code == -3): # current density too high
        if (self.target_option == "max_K_lse"):
          print("Current density too high.")
          # Target has been bracketed. Decrease interval.
          if (self.increased_target_current):
            self.current_factor = self.current_factor*0.5
            print("current_factor is now: " + str(self.current_factor))
          self.target_value = (1.0-self.current_factor)*self.target_value
          print("Trying again with target_value = " + str(self.target_value))
          os.chdir('..')
          self.decreased_target_current = True
          self.evaluateRegcoil(omegas_sensitivity_new,self.target_value)
        elif (self.target_option == 'chi2_B'):
          print("chi2_B too low.")
          # Target has been bracketed. Decrease interval.
          if (self.decreased_target_current):
            self.current_factor = self.current_factor*0.5
            print("current_factor is now: " + str(self.current_factor))
          self.target_value = (1.0+self.current_factor)*self.target_value
          print("Trying again with target_value = " + str(self.target_value)
          os.chdir('..'))
          self.increased_target_current = True
          self.evaluateRegcoil(omegas_sensitivity_new,self.target_value)
        else:
					print("Error! Incorrect target_option: "+str(self.target_option))
					sys.exit(0)
      else:
        sys.exit(0)

  # Creates new nescin file with current Fourier coefficients
  # old_nescin and new_nescin are names of nescin files
  # old_nescin is a copy of the old_nescin file in current directory
  # new_nescin file is desired name for new nescin file with
  # updated Fourier coefficients
  def create_nescin(self,old_nescin,new_nescin):
    currentMode = -1
    newFile = open(new_nescin, "w")
    oldFile = open(old_nescin, "r")

    for line in oldFile:
      if re.match("------ Current Surface",line):
        newFile.write(line)
        lineToWrite = next(oldFile)
        newFile.write(lineToWrite)
        lineToWrite = next(oldFile)
        newFile.write(str(self.nmodes) + "\n")
        lineToWrite = next(oldFile)
        newFile.write(lineToWrite)
        lineToWrite = next(oldFile)
        newFile.write(lineToWrite)
        currentMode = 0
        break
      else:
        newFile.write(line)
  
    oldFile.close()
          
    while(currentMode >= 0):
      lineToWrite = "\t" + str(int(self.xm[currentMode])) + "\t" \
        + str(int(self.xn[currentMode])) + "\t" + str(self.rmncs[currentMode]) \
        + "\t" + str(self.zmnss[currentMode]) + "\t" + str(self.rmnss[currentMode]) \
        + "\t" + str(self.zmncs[currentMode]) + "\n"
      newFile.write(lineToWrite)
      if (currentMode < self.nmodes-1):
        currentMode = currentMode + 1
      else:
        break

    newFile.close()

def evaluateFunctionRegcoil(omegas_sensitivity_new, nescinObject):
  # Check if function has already been evaluated
  if (nescinObject.evaluated == False or not np.array_equal(omegas_sensitivity_new,nescinObject.omegas_sensitivity)):
    if (nescinObject.general_option > 3):
      nescinObject.evaluateRegcoil(omegas_sensitivity_new,nescinObject.target_value)
    else:
      nescinObject.evaluateRegcoil(omegas_sensitivity_new)
  return nescinObject.objective_function

def evaluateGradientRegcoil(omegas_sensitivity_new, nescinObject):
  if (nescinObject.evaluated == False or not np.array_equal(omegas_sensitivity_new,nescinObject.omegas_sensitivity)):
    if (nescinObject.general_option > 3):
      nescinObject.evaluateRegcoil(omegas_sensitivity_new,nescinObject.target_value)
    else:
      nescinObject.evaluateRegcoil(omegas_sensitivity_new)
  return np.array(nescinObject.dobjective_functiondomegas_sensitivity)

## Testing ##
if __name__ == "__main__":
  nescinObject = coilFourier(int(sys.argv[1]),int(sys.argv[2]),sys.argv[3])
#  print(nescinObject.xn)
#  print(nescinObject.xm)
#  print(nescinObject.nmax)
#  print(nescinObject.mmax)

  file = "nescin.w7x_winding_surface_from_Drevlak_235"
  path = os.getcwd()
  filename = path + "/" + file
  nescinObject.set_Fourier_from_nescin(filename)
  print(nescinObject.xn)
  print(nescinObject.xn_sensitivity)
  omegas_old = nescinObject.omegas
  print(omegas_old)
  new_omegas = nescinObject.omegas_sensitivity
  for imn in range(0,nescinObject.nmodes_sensitivity):
    if (nescinObject.xm_sensitivity[imn] == 6 and nescinObject.xn_sensitivity[imn] == 3):
      print(new_omegas[2*imn])
      print(new_omegas[2*imn+1])
      new_omegas[2*imn] = 0
      new_omegas[2*imn+1] = 0
      print("mode obtained!"
  for imn in range(0,nescinObject.nmodes):
    if (nescinObject.xm[imn] == 6 and nescinObject.xn[imn] == 3):
      print("mode obtained in omega!")
      print(omegas_old[2*imn])
      print(omegas_old[2*imn+1])

  nescinObject.set_omegas_sensitivity(new_omegas)
  print(nescinObject.omegas)
  nescinObject.compute_spectral_norm()
  print(nescinObject.spectral_norm)
  print(nescinObject.dspectral_normdomegas)
